//! Raw-record BAM edit path.
//!
//! noodles' `bam::Record` holds the on-disk record bytes but never exposes them
//! (`Record(pub(crate) Vec<u8>)`, no getter, no `From<Vec<u8>>`), and its writer
//! re-encodes every record through the BAM codec — decoding and re-emitting
//! seq/qual/cigar even when a tool only flips a fixed-offset field. For tools
//! that edit a few header fields (markdup sets the 0x400 duplicate bit; fixmate
//! rewrites mate fields + MC/MQ), that round-trip is the bottleneck.
//!
//! [`RawRecord`] owns a record's bytes (the buffer after the 4-byte
//! `block_size`, exactly as the BGZF stream stores it) and edits fixed-offset
//! fields and the aux tail in place. seq/qual/cigar/name are never decoded —
//! they are passed through byte-for-byte. [`read_record`] / [`write_record`]
//! move records over the BGZF stream directly (block_size + payload), reusing
//! noodles' `bam::io::Reader`/`Writer` only for the header.
//!
//! Field offsets follow the BAM spec record layout (SAMv1 §4.2), measured from
//! the start of the payload (after `block_size`):
//!
//! ```text
//! refID@0 pos@4 l_read_name@8 mapq@9 bin@10 n_cigar@12 flag@14 l_seq@16
//! next_refID@20 next_pos@24 tlen@28
//! read_name(l_read_name) cigar(4*n_cigar) seq((l_seq+1)/2) qual(l_seq) aux
//! ```

use std::io::{self, Read, Write};

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

/// The 0x400 PCR/optical duplicate flag bit.
pub const FLAG_DUPLICATE: u16 = 0x400;

/// An owned BAM record's raw payload bytes (everything after `block_size`).
///
/// Edits operate on fixed-offset fields and the aux tail; the variable-length
/// region (name/cigar/seq/qual) is never decoded. Construct via [`read_record`]
/// or [`RawRecord::default`]; emit via [`write_record`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RawRecord {
    bytes: Vec<u8>,
}

impl RawRecord {
    /// The raw payload bytes (after `block_size`), as stored on disk.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    fn u16_at(&self, off: usize) -> u16 {
        u16::from_le_bytes(self.bytes[off..off + 2].try_into().unwrap())
    }

    fn i32_at(&self, off: usize) -> i32 {
        i32::from_le_bytes(self.bytes[off..off + 4].try_into().unwrap())
    }

    fn set_i32_at(&mut self, off: usize, value: i32) {
        self.bytes[off..off + 4].copy_from_slice(&value.to_le_bytes());
    }

    /// FLAG (offset 14, u16 LE).
    pub fn flags(&self) -> u16 {
        self.u16_at(FLAG)
    }

    /// refID (offset 0, i32 LE; -1 = unmapped/none).
    pub fn reference_sequence_id(&self) -> i32 {
        self.i32_at(REF_ID)
    }

    /// pos (offset 4, i32 LE, 0-based; -1 = none).
    pub fn alignment_start(&self) -> i32 {
        self.i32_at(POS)
    }

    /// MAPQ (offset 9, u8).
    pub fn mapping_quality(&self) -> u8 {
        self.bytes[MAPQ]
    }

    /// next_refID (offset 20, i32 LE).
    pub fn mate_reference_sequence_id(&self) -> i32 {
        self.i32_at(NEXT_REF_ID)
    }

    /// next_pos (offset 24, i32 LE, 0-based).
    pub fn mate_alignment_start(&self) -> i32 {
        self.i32_at(NEXT_POS)
    }

    /// tlen (offset 28, i32 LE).
    pub fn template_length(&self) -> i32 {
        self.i32_at(TLEN)
    }

    fn name_len(&self) -> usize {
        usize::from(self.bytes[L_READ_NAME])
    }

    fn cigar_op_count(&self) -> usize {
        usize::from(self.u16_at(N_CIGAR))
    }

    fn base_count(&self) -> usize {
        usize::try_from(u32::from_le_bytes(
            self.bytes[L_SEQ..L_SEQ + 4].try_into().unwrap(),
        ))
        .unwrap()
    }

    /// read_name without the trailing NUL. `*` (the BAM "no name" sentinel) is
    /// returned as-is — callers that need to distinguish handle it.
    pub fn name(&self) -> &[u8] {
        let end = FIXED_HEAD + self.name_len();
        let raw = &self.bytes[FIXED_HEAD..end];
        raw.strip_suffix(&[0]).unwrap_or(raw)
    }

    /// CIGAR operations as `(kind, len)` where `kind` is the 0..=8 op code
    /// (0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X). The packed 4-byte encoding is read
    /// directly, never materialised into noodles types.
    pub fn cigar_ops(&self) -> impl Iterator<Item = (u8, u32)> + '_ {
        let start = FIXED_HEAD + self.name_len();
        let n = self.cigar_op_count();
        (0..n).map(move |i| {
            let off = start + i * 4;
            let raw = u32::from_le_bytes(self.bytes[off..off + 4].try_into().unwrap());
            ((raw & 0xf) as u8, raw >> 4)
        })
    }

    fn quality_start(&self) -> usize {
        FIXED_HEAD + self.name_len() + self.cigar_op_count() * 4 + self.base_count().div_ceil(2)
    }

    /// Per-base quality scores (Phred, not ASCII-offset). A fully-`0xff` block
    /// (the BAM "missing qualities" sentinel) yields an empty slice.
    pub fn quality_scores(&self) -> &[u8] {
        let base_count = self.base_count();
        let start = self.quality_start();
        let qual = &self.bytes[start..start + base_count];
        if qual.iter().all(|&b| b == 0xff) {
            &[]
        } else {
            qual
        }
    }

    fn aux_start(&self) -> usize {
        let base_count = self.base_count();
        FIXED_HEAD
            + self.name_len()
            + self.cigar_op_count() * 4
            + base_count.div_ceil(2)
            + base_count
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

    /// Set refID (offset 0). The unmapped-read salvage step copies a mapped
    /// mate's reference onto its unmapped partner so a coordinate sort keeps the
    /// pair adjacent.
    pub fn set_reference_sequence_id(&mut self, value: i32) {
        self.set_i32_at(REF_ID, value);
    }

    /// Set pos (offset 4, 0-based). Paired with [`set_reference_sequence_id`]
    /// for the unmapped-read salvage step.
    pub fn set_alignment_start(&mut self, value: i32) {
        self.set_i32_at(POS, value);
    }

    /// Set next_refID (offset 20).
    pub fn set_mate_reference_sequence_id(&mut self, value: i32) {
        self.set_i32_at(NEXT_REF_ID, value);
    }

    /// Set next_pos (offset 24).
    pub fn set_mate_alignment_start(&mut self, value: i32) {
        self.set_i32_at(NEXT_POS, value);
    }

    /// Set tlen (offset 28).
    pub fn set_template_length(&mut self, value: i32) {
        self.set_i32_at(TLEN, value);
    }

    /// Byte range of the aux field with `tag`, including the 2-byte tag, the
    /// 1-byte type code, and the value. `None` if absent.
    fn aux_field_range(&self, tag: [u8; 2]) -> Option<std::ops::Range<usize>> {
        let mut pos = self.aux_start();
        let end = self.bytes.len();
        while pos + 3 <= end {
            let field_tag = [self.bytes[pos], self.bytes[pos + 1]];
            let type_code = self.bytes[pos + 2];
            let value_len = aux_value_len(&self.bytes, pos + 3, type_code)?;
            let field_end = pos + 3 + value_len;
            if field_tag == tag {
                return Some(pos..field_end);
            }
            pos = field_end;
        }
        None
    }

    /// Aux field value bytes (type code excluded) for `tag`, or `None` if absent.
    pub fn aux_value(&self, tag: [u8; 2]) -> Option<&[u8]> {
        let range = self.aux_field_range(tag)?;
        Some(&self.bytes[range.start + 3..range.end])
    }

    /// The aux field's BAM type code (`A`/`c`/`C`/`s`/`S`/`i`/`I`/`f`/`Z`/`H`/`B`)
    /// for `tag`, or `None` if absent.
    pub fn aux_type(&self, tag: [u8; 2]) -> Option<u8> {
        let range = self.aux_field_range(tag)?;
        Some(self.bytes[range.start + 2])
    }

    /// Append a complete aux field: 2-byte `tag`, 1-byte `type_code`, then
    /// `value`. The caller is responsible for the value's on-disk encoding
    /// matching `type_code`. The field is added at the end of the aux tail.
    pub fn append_aux(&mut self, tag: [u8; 2], type_code: u8, value: &[u8]) {
        self.bytes.push(tag[0]);
        self.bytes.push(tag[1]);
        self.bytes.push(type_code);
        self.bytes.extend_from_slice(value);
    }

    /// Remove the aux field with `tag`. No-op if absent. Returns whether a field
    /// was removed.
    pub fn remove_aux(&mut self, tag: [u8; 2]) -> bool {
        match self.aux_field_range(tag) {
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
    pub fn set_aux(&mut self, tag: [u8; 2], type_code: u8, value: &[u8]) {
        self.remove_aux(tag);
        self.append_aux(tag, type_code, value);
    }
}

/// Length in bytes of an aux value (excluding the type code) starting at `pos`.
/// `None` on a malformed/truncated field — callers treat that as a hard error.
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
            Some(1 + 4 + count * elem)
        }
        _ => None,
    }
}

/// Read one raw record's payload from a BGZF `Read` stream into `dst`, returning
/// the payload length. Returns `0` at end of records (a 0 `block_size`, which
/// the BGZF EOF block presents). The reader must be positioned after the header
/// (call `bam::io::Reader::read_header` first, then pass `reader.get_mut()`).
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
    Ok(block_size)
}

/// Write one raw record's payload to a BGZF `Write` stream, prefixed with its
/// `block_size` (u32 LE). Pass `writer.get_mut()` of a `bam::io::Writer` whose
/// header has already been written.
pub fn write_record<W: Write>(writer: &mut W, record: &RawRecord) -> Result<()> {
    let block_size = u32::try_from(record.bytes.len())
        .map_err(|e| RsomicsError::InvalidInput(format!("record too large: {e}")))?;
    writer
        .write_all(&block_size.to_le_bytes())
        .map_err(RsomicsError::Io)?;
    writer.write_all(&record.bytes).map_err(RsomicsError::Io)?;
    Ok(())
}

/// Fill `buf` from `reader`, returning `false` on a clean EOF before any byte is
/// read and `true` once `buf` is full. A partial read (EOF mid-buffer) is a hard
/// error.
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
