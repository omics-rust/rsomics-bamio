//! Validated owned and borrowing views over BAM record payloads.

use std::io::{self, BufRead, Read, Write};
use std::mem;

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, sam};
use rsomics_common::{Result, RsomicsError};

const REF_ID: usize = 0;
const POS: usize = 4;
const L_READ_NAME: usize = 8;
const MAPQ: usize = 9;
const N_CIGAR: usize = 12;
const FLAG: usize = 14;
const L_SEQ: usize = 16;
const NEXT_REF_ID: usize = 20;
const NEXT_POS: usize = 24;
const TLEN: usize = 28;
const FIXED_HEAD: usize = 32;
const CIGAR_REFERENCE_SKIP: u8 = 3;
const CIGAR_SOFT_CLIP: u8 = 4;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct RecordLayout {
    seq_start: u32,
    base_count: u32,
}

impl RecordLayout {
    fn seq_start(self) -> usize {
        self.seq_start as usize
    }

    fn base_count(self) -> usize {
        self.base_count as usize
    }

    fn quality_start(self) -> usize {
        self.seq_start() + self.base_count().div_ceil(2)
    }

    fn aux_start(self) -> usize {
        self.quality_start() + self.base_count()
    }
}

/// The 0x400 PCR/optical duplicate flag bit.
pub const FLAG_DUPLICATE: u16 = 0x400;

fn u16_at(bytes: &[u8], off: usize) -> u16 {
    u16::from_le_bytes(bytes[off..off + 2].try_into().unwrap())
}

fn i32_at(bytes: &[u8], off: usize) -> i32 {
    i32::from_le_bytes(bytes[off..off + 4].try_into().unwrap())
}

fn payload_name_len(bytes: &[u8]) -> usize {
    usize::from(bytes[L_READ_NAME])
}

fn payload_cigar_op_count(bytes: &[u8]) -> usize {
    usize::from(u16_at(bytes, N_CIGAR))
}

fn payload_base_count(bytes: &[u8]) -> usize {
    usize::try_from(u32::from_le_bytes(
        bytes[L_SEQ..L_SEQ + 4].try_into().unwrap(),
    ))
    .unwrap()
}

fn payload_name(bytes: &[u8]) -> &[u8] {
    let end = FIXED_HEAD + payload_name_len(bytes);
    let raw = &bytes[FIXED_HEAD..end];
    raw.strip_suffix(&[0]).unwrap_or(raw)
}

fn payload_cigar_ops(bytes: &[u8]) -> impl Iterator<Item = (u8, u32)> + '_ {
    let start = FIXED_HEAD + payload_name_len(bytes);
    let n = payload_cigar_op_count(bytes);
    (0..n).map(move |i| {
        let off = start + i * 4;
        let raw = u32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
        ((raw & 0xf) as u8, raw >> 4)
    })
}

fn payload_decoded_cigar(bytes: &[u8], layout: RecordLayout) -> Result<Vec<(u8, u32)>> {
    let raw = payload_cigar_ops(bytes).collect::<Vec<_>>();
    let cigar = if raw.len() == 2
        && raw[0] == (CIGAR_SOFT_CLIP, layout.base_count)
        && raw[1].0 == CIGAR_REFERENCE_SKIP
        && payload_aux_type(bytes, layout, *b"CG") == Some(b'B')
    {
        payload_long_cigar(bytes, layout, raw[1].1)?
    } else {
        raw
    };
    for &(kind, length) in &cigar {
        if kind > 8 {
            return Err(invalid_cigar(
                bytes,
                format!("unsupported CIGAR operation {kind}"),
            ));
        }
        if length == 0 {
            return Err(invalid_cigar(
                bytes,
                "CIGAR contains a zero-length operation",
            ));
        }
    }
    Ok(cigar)
}

fn payload_long_cigar(
    bytes: &[u8],
    layout: RecordLayout,
    placeholder_span: u32,
) -> Result<Vec<(u8, u32)>> {
    let value = payload_aux_value(bytes, layout, *b"CG")
        .ok_or_else(|| invalid_cigar(bytes, "CG long-CIGAR field is absent"))?;
    if value.len() < 5 || value[0] != b'I' {
        return Err(invalid_cigar(
            bytes,
            "CG long-CIGAR field is not a B,I array",
        ));
    }
    let count = usize::try_from(u32::from_le_bytes(
        value[1..5]
            .try_into()
            .expect("long-CIGAR count has four bytes"),
    ))
    .map_err(|error| invalid_cigar(bytes, error))?;
    if count <= usize::from(u16::MAX) {
        return Err(invalid_cigar(
            bytes,
            "CG long-CIGAR field does not exceed the BAM operation limit",
        ));
    }
    let expected = count
        .checked_mul(4)
        .and_then(|length| length.checked_add(5))
        .ok_or_else(|| invalid_cigar(bytes, "CG long-CIGAR length overflows"))?;
    if value.len() != expected {
        return Err(invalid_cigar(
            bytes,
            "CG long-CIGAR field has the wrong length",
        ));
    }
    let cigar = value[5..]
        .chunks_exact(4)
        .map(|chunk| {
            let packed = u32::from_le_bytes(
                chunk
                    .try_into()
                    .expect("long-CIGAR operation has four bytes"),
            );
            ((packed & 0x0f) as u8, packed >> 4)
        })
        .collect::<Vec<_>>();
    let reference_span = cigar
        .iter()
        .filter(|(kind, _)| matches!(kind, 0 | 2 | 3 | 7 | 8))
        .try_fold(0u64, |span, (_, length)| {
            span.checked_add(u64::from(*length))
        })
        .ok_or_else(|| invalid_cigar(bytes, "CG long-CIGAR reference span overflows"))?;
    if reference_span != u64::from(placeholder_span) {
        return Err(invalid_cigar(
            bytes,
            "CG long-CIGAR reference span differs from its placeholder",
        ));
    }
    Ok(cigar)
}

fn invalid_cigar(bytes: &[u8], reason: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "read {}: {reason}",
        String::from_utf8_lossy(payload_name(bytes))
    ))
}

fn payload_seq_nibble(bytes: &[u8], layout: RecordLayout, i: usize) -> u8 {
    let byte = bytes[layout.seq_start() + i / 2];
    if i.is_multiple_of(2) {
        byte >> 4
    } else {
        byte & 0x0f
    }
}

fn payload_quality_scores(bytes: &[u8], layout: RecordLayout) -> &[u8] {
    let qual = &bytes[layout.quality_start()..layout.aux_start()];
    if qual.iter().all(|&b| b == 0xff) {
        &[]
    } else {
        qual
    }
}

/// Byte range of the aux field with `tag`, including the 2-byte tag, the 1-byte
/// type code, and the value. `None` if absent.
fn payload_aux_field_range(
    bytes: &[u8],
    layout: RecordLayout,
    tag: [u8; 2],
) -> Option<std::ops::Range<usize>> {
    let mut pos = layout.aux_start();
    let end = bytes.len();
    while pos + 3 <= end {
        let field_tag = [bytes[pos], bytes[pos + 1]];
        let type_code = bytes[pos + 2];
        let value_len = aux_value_len(bytes, pos + 3, type_code)?;
        let field_end = pos.checked_add(3)?.checked_add(value_len)?;
        if field_end > end {
            return None;
        }
        if field_tag == tag {
            return Some(pos..field_end);
        }
        pos = field_end;
    }
    None
}

fn payload_aux_value(bytes: &[u8], layout: RecordLayout, tag: [u8; 2]) -> Option<&[u8]> {
    let range = payload_aux_field_range(bytes, layout, tag)?;
    Some(&bytes[range.start + 3..range.end])
}

fn payload_aux_type(bytes: &[u8], layout: RecordLayout, tag: [u8; 2]) -> Option<u8> {
    let range = payload_aux_field_range(bytes, layout, tag)?;
    Some(bytes[range.start + 2])
}

macro_rules! raw_field_accessors {
    (owned $ty:ty) => {
        raw_field_accessors!(@impl impl $ty);
    };
    (borrowed <$life:lifetime> $ty:ty) => {
        raw_field_accessors!(@impl impl<$life> $ty);
    };
    (@impl $($head:tt)+) => {
        $($head)+ {
            fn payload_bytes(&self) -> &[u8] {
                &self.bytes[..]
            }

            fn layout(&self) -> RecordLayout {
                self.layout
            }

            pub fn as_bytes(&self) -> &[u8] {
                self.payload_bytes()
            }

            pub fn flags(&self) -> u16 {
                u16_at(self.payload_bytes(), FLAG)
            }

            pub fn reference_sequence_id(&self) -> i32 {
                i32_at(self.payload_bytes(), REF_ID)
            }

            /// 0-based; -1 = unmapped.
            pub fn alignment_start(&self) -> i32 {
                i32_at(self.payload_bytes(), POS)
            }

            pub fn mapping_quality(&self) -> u8 {
                self.payload_bytes()[MAPQ]
            }

            pub fn mate_reference_sequence_id(&self) -> i32 {
                i32_at(self.payload_bytes(), NEXT_REF_ID)
            }

            pub fn mate_alignment_start(&self) -> i32 {
                i32_at(self.payload_bytes(), NEXT_POS)
            }

            pub fn template_length(&self) -> i32 {
                i32_at(self.payload_bytes(), TLEN)
            }

            /// The read name without its trailing NUL.
            pub fn name(&self) -> &[u8] {
                payload_name(self.payload_bytes())
            }

            /// CIGAR operations as BAM `(kind, length)` pairs.
            pub fn cigar_ops(&self) -> impl Iterator<Item = (u8, u32)> + '_ {
                payload_cigar_ops(self.payload_bytes())
            }

            /// The decoded BAM CIGAR, including a `CG:B,I` long CIGAR.
            pub fn decoded_cigar(&self) -> Result<Vec<(u8, u32)>> {
                payload_decoded_cigar(self.payload_bytes(), self.layout())
            }

            /// The number of query bases.
            pub fn sequence_len(&self) -> usize {
                self.layout().base_count()
            }

            /// The BAM `seq_nt16` code at query index `i`.
            pub fn seq_nibble(&self, i: usize) -> u8 {
                payload_seq_nibble(self.payload_bytes(), self.layout(), i)
            }

            /// The packed BAM sequence bytes.
            pub fn seq_bytes_packed(&self) -> &[u8] {
                let bytes = self.payload_bytes();
                let layout = self.layout();
                &bytes[layout.seq_start()..layout.quality_start()]
            }

            /// Raw Phred scores, or an empty slice for the missing-score sentinel.
            pub fn quality_scores(&self) -> &[u8] {
                payload_quality_scores(self.payload_bytes(), self.layout())
            }

            /// Auxiliary value bytes without the type code.
            pub fn aux_value(&self, tag: [u8; 2]) -> Option<&[u8]> {
                payload_aux_value(self.payload_bytes(), self.layout(), tag)
            }

            /// The BAM type code for an auxiliary field.
            pub fn aux_type(&self, tag: [u8; 2]) -> Option<u8> {
                payload_aux_type(self.payload_bytes(), self.layout(), tag)
            }
        }
    };
}

/// An owned, validated BAM record payload.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RawRecord {
    bytes: Vec<u8>,
    layout: RecordLayout,
}

raw_field_accessors!(owned RawRecord);

impl Default for RawRecord {
    fn default() -> Self {
        Self {
            bytes: vec![
                0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x02, 0xff, 0x48, 0x12, 0x00, 0x00,
                0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                0x00, 0x00, 0x00, 0x00, b'*', 0x00,
            ],
            layout: RecordLayout {
                seq_start: 34,
                base_count: 0,
            },
        }
    }
}

impl TryFrom<Vec<u8>> for RawRecord {
    type Error = RsomicsError;

    fn try_from(bytes: Vec<u8>) -> Result<Self> {
        let layout = validate_payload(&bytes)?;
        Ok(Self { bytes, layout })
    }
}

/// Encodes alignment records into validated owned BAM payloads.
pub struct RawRecordEncoder {
    writer: bam::io::Writer<Vec<u8>>,
}

impl Default for RawRecordEncoder {
    fn default() -> Self {
        Self {
            writer: bam::io::Writer::from(Vec::new()),
        }
    }
}

impl RawRecordEncoder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn encode(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> Result<RawRecord> {
        self.writer
            .write_alignment_record(header, record)
            .map_err(RsomicsError::Io)?;
        let mut framed = mem::take(self.writer.get_mut());
        let size = framed
            .get(..4)
            .and_then(|bytes| <[u8; 4]>::try_from(bytes).ok())
            .map(u32::from_le_bytes)
            .and_then(|value| usize::try_from(value).ok())
            .ok_or_else(invalid_encoded_size)?;
        if size != framed.len() - 4 {
            return Err(invalid_encoded_size());
        }
        framed.drain(..4);
        RawRecord::try_from(framed)
    }
}

pub(crate) fn from_htslib_record(record: &rust_htslib::bam::Record) -> Result<RawRecord> {
    let inner = record.inner();
    let core = &inner.core;
    let data_len = usize::try_from(inner.l_data)
        .map_err(|_| invalid_record("HTSlib record has a negative data length"))?;
    let data = if data_len == 0 {
        &[]
    } else {
        // SAFETY: HTSlib owns `data` and guarantees `l_data` initialized bytes until mutation.
        unsafe { std::slice::from_raw_parts(inner.data, data_len) }
    };
    let qname_storage = usize::from(core.l_qname);
    let extra_nuls = usize::from(core.l_extranul);
    let qname_len = qname_storage
        .checked_sub(extra_nuls)
        .filter(|&len| len > 0 && len <= u8::MAX as usize)
        .ok_or_else(|| invalid_record("HTSlib record has an invalid read name length"))?;
    let cigar_len = usize::try_from(core.n_cigar)
        .ok()
        .and_then(|len| len.checked_mul(4))
        .ok_or_else(|| invalid_record("HTSlib record CIGAR length overflows"))?;
    let cigar_end = qname_storage
        .checked_add(cigar_len)
        .filter(|&end| end <= data.len())
        .ok_or_else(|| invalid_record("HTSlib record data is truncated"))?;
    let position = i32::try_from(core.pos)
        .map_err(|_| invalid_record("HTSlib record position exceeds BAM limits"))?;
    let mate_position = i32::try_from(core.mpos)
        .map_err(|_| invalid_record("HTSlib mate position exceeds BAM limits"))?;
    let template_length = i32::try_from(core.isize_)
        .map_err(|_| invalid_record("HTSlib template length exceeds BAM limits"))?;

    let long_cigar = core.n_cigar > u16::MAX.into();
    let encoded_cigar_count = if long_cigar { 2 } else { core.n_cigar as u16 };
    let mut bytes = Vec::with_capacity(32 + data_len - extra_nuls + usize::from(long_cigar) * 16);
    bytes.extend_from_slice(&core.tid.to_le_bytes());
    bytes.extend_from_slice(&position.to_le_bytes());
    bytes.extend_from_slice(
        &((u32::from(core.bin) << 16) | (u32::from(core.qual) << 8) | qname_len as u32)
            .to_le_bytes(),
    );
    bytes.extend_from_slice(
        &((u32::from(core.flag) << 16) | u32::from(encoded_cigar_count)).to_le_bytes(),
    );
    bytes.extend_from_slice(&core.l_qseq.to_le_bytes());
    bytes.extend_from_slice(&core.mtid.to_le_bytes());
    bytes.extend_from_slice(&mate_position.to_le_bytes());
    bytes.extend_from_slice(&template_length.to_le_bytes());
    bytes.extend_from_slice(&data[..qname_len]);

    if long_cigar {
        let read_len = u32::try_from(core.l_qseq)
            .ok()
            .filter(|&len| len < 1 << 28)
            .ok_or_else(|| invalid_record("HTSlib record sequence length exceeds BAM limits"))?;
        let reference_len = data[qname_storage..cigar_end]
            .chunks_exact(4)
            .try_fold(0_u32, |sum, chunk| {
                let op = u32::from_ne_bytes(chunk.try_into().unwrap());
                if matches!(op & 0x0f, 0 | 2 | 3 | 7 | 8) {
                    sum.checked_add(op >> 4)
                } else {
                    Some(sum)
                }
            })
            .filter(|&len| len < 1 << 28)
            .ok_or_else(|| invalid_record("HTSlib record reference span exceeds BAM limits"))?;
        bytes.extend_from_slice(&((read_len << 4) | 4).to_le_bytes());
        bytes.extend_from_slice(&((reference_len << 4) | 3).to_le_bytes());
        bytes.extend_from_slice(&data[cigar_end..]);
        bytes.extend_from_slice(b"CGBI");
        bytes.extend_from_slice(&core.n_cigar.to_le_bytes());
        for chunk in data[qname_storage..cigar_end].chunks_exact(4) {
            bytes.extend_from_slice(&u32::from_ne_bytes(chunk.try_into().unwrap()).to_le_bytes());
        }
    } else {
        bytes.extend_from_slice(&data[qname_storage..]);
    }

    RawRecord::try_from(bytes)
}

fn invalid_encoded_size() -> RsomicsError {
    RsomicsError::InvalidInput("BAM encoder produced an invalid record size".to_owned())
}

impl RawRecord {
    fn set_i32_at(&mut self, off: usize, value: i32) {
        self.bytes[off..off + 4].copy_from_slice(&value.to_le_bytes());
    }

    /// Mutable per-base quality bytes, including the missing-score sentinel.
    pub fn quality_scores_mut(&mut self) -> &mut [u8] {
        &mut self.bytes[self.layout.quality_start()..self.layout.aux_start()]
    }

    /// Mutable packed BAM sequence bytes.
    pub fn seq_bytes_mut(&mut self) -> &mut [u8] {
        &mut self.bytes[self.layout.seq_start()..self.layout.quality_start()]
    }

    /// Set the given FLAG bits (offset 14). Other bits are left untouched.
    pub fn set_flag_bits(&mut self, bits: u16) {
        let new = self.flags() | bits;
        self.bytes[FLAG..FLAG + 2].copy_from_slice(&new.to_le_bytes());
    }

    /// Clear the given FLAG bits (offset 14). Other bits are left untouched.
    pub fn clear_flag_bits(&mut self, bits: u16) {
        let new = self.flags() & !bits;
        self.bytes[FLAG..FLAG + 2].copy_from_slice(&new.to_le_bytes());
    }

    /// Set the reference-sequence ID.
    pub fn set_reference_sequence_id(&mut self, value: i32) {
        self.set_i32_at(REF_ID, value);
    }

    pub fn set_alignment_start(&mut self, value: i32) {
        self.set_i32_at(POS, value);
    }

    pub fn set_mate_reference_sequence_id(&mut self, value: i32) {
        self.set_i32_at(NEXT_REF_ID, value);
    }

    pub fn set_mate_alignment_start(&mut self, value: i32) {
        self.set_i32_at(NEXT_POS, value);
    }

    pub fn set_template_length(&mut self, value: i32) {
        self.set_i32_at(TLEN, value);
    }

    /// Append a complete aux field: 2-byte `tag`, 1-byte `type_code`, then
    /// `value`. The caller is responsible for the value's on-disk encoding
    /// matching `type_code`. The field is added at the end of the aux tail.
    pub fn append_aux(&mut self, tag: [u8; 2], type_code: u8, value: &[u8]) -> Result<()> {
        validate_aux_value(type_code, value)?;
        self.bytes.push(tag[0]);
        self.bytes.push(tag[1]);
        self.bytes.push(type_code);
        self.bytes.extend_from_slice(value);
        Ok(())
    }

    /// Remove the aux field with `tag`. No-op if absent. Returns whether a field
    /// was removed.
    pub fn remove_aux(&mut self, tag: [u8; 2]) -> bool {
        match payload_aux_field_range(&self.bytes, self.layout, tag) {
            Some(range) => {
                self.bytes.drain(range);
                true
            }
            None => false,
        }
    }

    /// Replace the aux field with `tag`, or append it if absent. Equivalent to
    /// [`remove_aux`](Self::remove_aux) followed by
    /// [`append_aux`](Self::append_aux).
    pub fn set_aux(&mut self, tag: [u8; 2], type_code: u8, value: &[u8]) -> Result<()> {
        validate_aux_value(type_code, value)?;
        self.remove_aux(tag);
        self.append_aux(tag, type_code, value)
    }

    /// Copy a record while retaining this allocation.
    #[inline]
    pub fn clone_from_raw(&mut self, other: &RawRecord) {
        self.bytes.clear();
        self.bytes.extend_from_slice(&other.bytes);
        self.layout = other.layout;
    }
}

fn aux_value_len(bytes: &[u8], pos: usize, type_code: u8) -> Option<usize> {
    match type_code {
        b'A' | b'c' | b'C' => Some(1),
        b's' | b'S' => Some(2),
        b'i' | b'I' | b'f' => Some(4),
        b'Z' | b'H' => {
            let nul = bytes[pos..].iter().position(|&b| b == 0)?;
            Some(nul + 1)
        }
        b'B' => {
            let sub = *bytes.get(pos)?;
            let count = u32::from_le_bytes(bytes.get(pos + 1..pos + 5)?.try_into().ok()?) as usize;
            let elem = match sub {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            count.checked_mul(elem)?.checked_add(5)
        }
        _ => None,
    }
}

fn validate_aux_value(type_code: u8, value: &[u8]) -> Result<()> {
    if aux_value_len(value, 0, type_code) == Some(value.len()) {
        Ok(())
    } else {
        Err(RsomicsError::InvalidInput(format!(
            "invalid BAM auxiliary value for type {}",
            char::from(type_code)
        )))
    }
}

fn validate_payload(bytes: &[u8]) -> Result<RecordLayout> {
    if bytes.len() < FIXED_HEAD {
        return Err(invalid_record("payload is shorter than the fixed fields"));
    }

    let name_len = payload_name_len(bytes);
    if name_len == 0 {
        return Err(invalid_record("read name length is zero"));
    }

    let name_end = FIXED_HEAD
        .checked_add(name_len)
        .ok_or_else(|| invalid_record("record layout overflows"))?;
    if name_end > bytes.len() || bytes[name_end - 1] != 0 {
        return Err(invalid_record(
            "read name is truncated or not NUL-terminated",
        ));
    }

    let cigar_len = payload_cigar_op_count(bytes)
        .checked_mul(4)
        .ok_or_else(|| invalid_record("record layout overflows"))?;
    let base_count = payload_base_count(bytes);
    let seq_start = name_end
        .checked_add(cigar_len)
        .ok_or_else(|| invalid_record("record layout overflows"))?;
    let quality_start = seq_start
        .checked_add(base_count.div_ceil(2))
        .ok_or_else(|| invalid_record("record layout overflows"))?;
    let aux_start = quality_start
        .checked_add(base_count)
        .ok_or_else(|| invalid_record("record layout overflows"))?;
    if aux_start > bytes.len() {
        return Err(invalid_record("variable-length fields are truncated"));
    }

    let layout = RecordLayout {
        seq_start: u32::try_from(seq_start)
            .map_err(|_| invalid_record("record layout overflows"))?,
        base_count: u32::try_from(base_count)
            .map_err(|_| invalid_record("record layout overflows"))?,
    };

    let mut pos = aux_start;
    while pos < bytes.len() {
        if bytes.len() - pos < 3 {
            return Err(invalid_record("auxiliary field header is truncated"));
        }
        let value_len = aux_value_len(bytes, pos + 3, bytes[pos + 2])
            .ok_or_else(|| invalid_record("auxiliary field is malformed"))?;
        pos = pos
            .checked_add(3)
            .and_then(|end| end.checked_add(value_len))
            .filter(|&end| end <= bytes.len())
            .ok_or_else(|| invalid_record("auxiliary field is truncated"))?;
    }

    Ok(layout)
}

fn invalid_record(message: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("invalid BAM record: {message}"))
}

/// Read and validate one BAM record payload, returning zero at EOF.
pub fn read_record<R: Read>(reader: &mut R, dst: &mut RawRecord) -> Result<usize> {
    let mut size_buf = [0u8; 4];
    if !read_exact_or_eof(reader, &mut size_buf).map_err(RsomicsError::Io)? {
        return Ok(0);
    }
    let block_size = u32::from_le_bytes(size_buf) as usize;
    if block_size == 0 {
        return Ok(0);
    }
    dst.bytes.resize(block_size, 0);
    reader
        .read_exact(&mut dst.bytes)
        .map_err(RsomicsError::Io)?;
    dst.layout = validate_payload(&dst.bytes)?;
    Ok(block_size)
}

/// A borrowed, validated BAM record payload.
#[derive(Clone, Copy, Debug)]
pub struct RecordRef<'a> {
    bytes: &'a [u8],
    layout: RecordLayout,
}

raw_field_accessors!(borrowed <'a> RecordRef<'a>);

impl<'a> RecordRef<'a> {
    /// Construct a borrowing view over a raw BAM record payload slice.
    pub fn from_bytes(bytes: &'a [u8]) -> Result<Self> {
        let layout = validate_payload(bytes)?;
        Ok(Self { bytes, layout })
    }

    /// The validated payload with the record's full borrow lifetime.
    pub fn payload(&self) -> &'a [u8] {
        self.bytes
    }
}

/// A borrowing record scanner over a buffered BAM stream.
pub struct RecordReader<'r, R: BufRead> {
    reader: &'r mut R,
    scratch: Vec<u8>,
    pending_consume: usize,
}

impl<'r, R: BufRead> RecordReader<'r, R> {
    /// Wrap a BGZF [`BufRead`] for borrowing record scans.
    pub fn new(reader: &'r mut R) -> Self {
        RecordReader {
            reader,
            scratch: Vec::new(),
            pending_consume: 0,
        }
    }

    /// Read and validate the next record.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<Option<RecordRef<'_>>> {
        if self.pending_consume != 0 {
            self.reader.consume(self.pending_consume);
            self.pending_consume = 0;
        }

        let mut size_buf = [0u8; 4];
        if !read_exact_or_eof(self.reader, &mut size_buf).map_err(RsomicsError::Io)? {
            return Ok(None);
        }
        let block_size = u32::from_le_bytes(size_buf) as usize;
        if block_size == 0 {
            return Ok(None);
        }

        let buffered = self.reader.fill_buf().map_err(RsomicsError::Io)?.len();
        if buffered >= block_size {
            // Consumption is deferred because the returned record borrows this block.
            self.pending_consume = block_size;
            let bytes = &self.reader.fill_buf().map_err(RsomicsError::Io)?[..block_size];
            let layout = validate_payload(bytes)?;
            Ok(Some(RecordRef { bytes, layout }))
        } else {
            self.scratch.resize(block_size, 0);
            self.reader
                .read_exact(&mut self.scratch)
                .map_err(RsomicsError::Io)?;
            let layout = validate_payload(&self.scratch)?;
            Ok(Some(RecordRef {
                bytes: &self.scratch,
                layout,
            }))
        }
    }
}

/// Write an owned raw BAM record with its block-size prefix.
pub fn write_record<W: Write>(writer: &mut W, record: &RawRecord) -> Result<()> {
    write_payload(writer, record.as_bytes())
}

/// Write a borrowed raw BAM record, including its block-size prefix.
pub fn write_record_ref<W: Write>(writer: &mut W, record: &RecordRef<'_>) -> Result<()> {
    write_payload(writer, record.payload())
}

fn write_payload<W: Write>(writer: &mut W, payload: &[u8]) -> Result<()> {
    let block_size = u32::try_from(payload.len())
        .map_err(|e| RsomicsError::InvalidInput(format!("record too large: {e}")))?;
    writer
        .write_all(&block_size.to_le_bytes())
        .map_err(RsomicsError::Io)?;
    writer.write_all(payload).map_err(RsomicsError::Io)?;
    Ok(())
}

fn read_exact_or_eof<R: Read>(reader: &mut R, buf: &mut [u8]) -> io::Result<bool> {
    let mut filled = 0;
    while filled < buf.len() {
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    if filled == 0 {
        Ok(false)
    } else if filled == buf.len() {
        Ok(true)
    } else {
        Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "truncated BAM record block size",
        ))
    }
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::record::{Cigar, CigarString};

    use super::*;

    #[test]
    fn converts_htslib_long_cigar_to_bam_encoding() {
        let cigar = CigarString(vec![Cigar::Match(1); usize::from(u16::MAX) + 1]);
        let sequence = vec![b'A'; cigar.len()];
        let qualities = vec![30; cigar.len()];
        let mut record = rust_htslib::bam::Record::new();
        record.set(b"long", Some(&cigar), &sequence, &qualities);
        record.set_tid(0);
        record.set_pos(0);

        let raw = from_htslib_record(&record).unwrap();
        assert_eq!(
            raw.cigar_ops().collect::<Vec<_>>(),
            [(4, 65_536), (3, 65_536)]
        );
        let decoded = raw.decoded_cigar().unwrap();
        assert_eq!(decoded.len(), 65_536);
        assert!(decoded.iter().all(|&op| op == (0, 1)));
    }
}
