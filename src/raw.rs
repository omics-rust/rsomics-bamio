//! Raw-record BAM read/edit path.
//!
//! noodles' `bam::Record` holds the on-disk record bytes but never exposes them
//! (`Record(pub(crate) Vec<u8>)`, no getter, no `From<Vec<u8>>`), and its writer
//! re-encodes every record through the BAM codec — decoding and re-emitting
//! seq/qual/cigar even when a tool only flips a fixed-offset field. For tools
//! that edit a few header fields (markdup sets the 0x400 duplicate bit; fixmate
//! rewrites mate fields + MC/MQ), that round-trip is the bottleneck.
//!
//! Two record views share one field-decoding core (the `payload_*` free
//! functions over `&[u8]`):
//!
//! - [`RawRecord`] **owns** a record's bytes (the buffer after the 4-byte
//!   `block_size`, exactly as the BGZF stream stores it). It supports in-place
//!   edits of fixed-offset fields and the aux tail. seq/qual/cigar/name are
//!   never decoded — they are passed through byte-for-byte. Use it when the
//!   record must outlive the next read, or is mutated.
//! - [`RecordRef`] **borrows** a record's bytes. For read-only scans
//!   ([`RecordReader`]) the record is read directly out of the BGZF
//!   reader's already-decompressed block buffer with no per-record allocation
//!   or copy whenever it lies within one block; only a block-straddling record
//!   spills into a caller-owned scratch buffer. It exposes the same accessors
//!   as `RawRecord`.
//!
//! [`read_record`] / [`write_record`] move owned records over the BGZF stream
//! directly (block_size + payload), reusing noodles' `bam::io::Reader`/`Writer`
//! only for the header.
//!
//! Field offsets follow the BAM spec record layout (SAMv1 §4.2), measured from
//! the start of the payload (after `block_size`):
//!
//! ```text
//! refID@0 pos@4 l_read_name@8 mapq@9 bin@10 n_cigar@12 flag@14 l_seq@16
//! next_refID@20 next_pos@24 tlen@28
//! read_name(l_read_name) cigar(4*n_cigar) seq((l_seq+1)/2) qual(l_seq) aux
//! ```

use std::io::{self, BufRead, Read, Write};

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

fn payload_quality_start(bytes: &[u8]) -> usize {
    FIXED_HEAD
        + payload_name_len(bytes)
        + payload_cigar_op_count(bytes) * 4
        + payload_base_count(bytes).div_ceil(2)
}

fn payload_seq_nibble(bytes: &[u8], i: usize) -> u8 {
    let seq_start = FIXED_HEAD + payload_name_len(bytes) + payload_cigar_op_count(bytes) * 4;
    let byte = bytes[seq_start + i / 2];
    if i.is_multiple_of(2) {
        byte >> 4
    } else {
        byte & 0x0f
    }
}

fn payload_quality_scores(bytes: &[u8]) -> &[u8] {
    let base_count = payload_base_count(bytes);
    let start = payload_quality_start(bytes);
    let qual = &bytes[start..start + base_count];
    if qual.iter().all(|&b| b == 0xff) {
        &[]
    } else {
        qual
    }
}

fn payload_aux_start(bytes: &[u8]) -> usize {
    let base_count = payload_base_count(bytes);
    FIXED_HEAD
        + payload_name_len(bytes)
        + payload_cigar_op_count(bytes) * 4
        + base_count.div_ceil(2)
        + base_count
}

/// Byte range of the aux field with `tag`, including the 2-byte tag, the 1-byte
/// type code, and the value. `None` if absent.
fn payload_aux_field_range(bytes: &[u8], tag: [u8; 2]) -> Option<std::ops::Range<usize>> {
    let mut pos = payload_aux_start(bytes);
    let end = bytes.len();
    while pos + 3 <= end {
        let field_tag = [bytes[pos], bytes[pos + 1]];
        let type_code = bytes[pos + 2];
        let value_len = aux_value_len(bytes, pos + 3, type_code)?;
        let field_end = pos + 3 + value_len;
        if field_tag == tag {
            return Some(pos..field_end);
        }
        pos = field_end;
    }
    None
}

fn payload_aux_value(bytes: &[u8], tag: [u8; 2]) -> Option<&[u8]> {
    let range = payload_aux_field_range(bytes, tag)?;
    Some(&bytes[range.start + 3..range.end])
}

fn payload_aux_type(bytes: &[u8], tag: [u8; 2]) -> Option<u8> {
    let range = payload_aux_field_range(bytes, tag)?;
    Some(bytes[range.start + 2])
}

/// Generate the read-only field accessors as inherent methods on both
/// [`RawRecord`] (owned) and [`RecordRef`] (borrowed). Each method delegates to
/// the `payload_*` free functions over `self.bytes`, keeping the decoding core
/// in one place while preserving the inherent-method API the `rsomics-bam-*`
/// tools already call (no trait import required by callers).
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

            /// The raw payload bytes (after `block_size`), as stored on disk.
            pub fn as_bytes(&self) -> &[u8] {
                self.payload_bytes()
            }

            /// FLAG (offset 14, u16 LE).
            pub fn flags(&self) -> u16 {
                u16_at(self.payload_bytes(), FLAG)
            }

            /// refID (offset 0, i32 LE; -1 = unmapped/none).
            pub fn reference_sequence_id(&self) -> i32 {
                i32_at(self.payload_bytes(), REF_ID)
            }

            /// pos (offset 4, i32 LE, 0-based; -1 = none).
            pub fn alignment_start(&self) -> i32 {
                i32_at(self.payload_bytes(), POS)
            }

            /// MAPQ (offset 9, u8).
            pub fn mapping_quality(&self) -> u8 {
                self.payload_bytes()[MAPQ]
            }

            /// next_refID (offset 20, i32 LE).
            pub fn mate_reference_sequence_id(&self) -> i32 {
                i32_at(self.payload_bytes(), NEXT_REF_ID)
            }

            /// next_pos (offset 24, i32 LE, 0-based).
            pub fn mate_alignment_start(&self) -> i32 {
                i32_at(self.payload_bytes(), NEXT_POS)
            }

            /// tlen (offset 28, i32 LE).
            pub fn template_length(&self) -> i32 {
                i32_at(self.payload_bytes(), TLEN)
            }

            /// read_name without the trailing NUL. `*` (the BAM "no name"
            /// sentinel) is returned as-is — callers that need to distinguish
            /// handle it.
            pub fn name(&self) -> &[u8] {
                payload_name(self.payload_bytes())
            }

            /// CIGAR operations as `(kind, len)` where `kind` is the 0..=8 op
            /// code (0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X). The packed 4-byte
            /// encoding is read directly, never materialised into noodles types.
            pub fn cigar_ops(&self) -> impl Iterator<Item = (u8, u32)> + '_ {
                payload_cigar_ops(self.payload_bytes())
            }

            /// Number of query bases (`l_seq`, offset 16). This is the SEQ/QUAL
            /// length, not the byte length of either packed field.
            pub fn sequence_len(&self) -> usize {
                payload_base_count(self.payload_bytes())
            }

            /// The 4-bit `seq_nt16` code (0..=15) of the query base at 0-based
            /// index `i`, read straight from the packed SEQ nibbles (high nibble
            /// = even index). The pileup base encoding maps this through
            /// `seq_nt16_str` (`=ACMGRSVTWYHKDBN`), so the engine carries the raw
            /// code rather than a decoded byte.
            pub fn seq_nibble(&self, i: usize) -> u8 {
                payload_seq_nibble(self.payload_bytes(), i)
            }

            /// The raw packed-nibble SEQ bytes as stored in BAM (2 bases per byte,
            /// high nibble = base at even query position). Length = `(seq_len+1)/2`.
            ///
            /// Prefer this over repeated `seq_nibble` calls in inner loops: it
            /// computes the sequence start offset once rather than per-call.
            pub fn seq_bytes_packed(&self) -> &[u8] {
                let bytes = self.payload_bytes();
                let seq_start = FIXED_HEAD
                    + payload_name_len(bytes)
                    + payload_cigar_op_count(bytes) * 4;
                let seq_len = payload_base_count(bytes);
                &bytes[seq_start..seq_start + seq_len.div_ceil(2)]
            }

            /// Per-base quality scores (Phred, not ASCII-offset). A fully-`0xff`
            /// block (the BAM "missing qualities" sentinel) yields an empty slice.
            pub fn quality_scores(&self) -> &[u8] {
                payload_quality_scores(self.payload_bytes())
            }

            /// Aux field value bytes (type code excluded) for `tag`, or `None` if
            /// absent.
            pub fn aux_value(&self, tag: [u8; 2]) -> Option<&[u8]> {
                payload_aux_value(self.payload_bytes(), tag)
            }

            /// The aux field's BAM type code
            /// (`A`/`c`/`C`/`s`/`S`/`i`/`I`/`f`/`Z`/`H`/`B`) for `tag`, or `None`
            /// if absent.
            pub fn aux_type(&self, tag: [u8; 2]) -> Option<u8> {
                payload_aux_type(self.payload_bytes(), tag)
            }
        }
    };
}

/// An owned BAM record's raw payload bytes (everything after `block_size`).
///
/// Edits operate on fixed-offset fields and the aux tail; the variable-length
/// region (name/cigar/seq/qual) is never decoded. Construct via [`read_record`]
/// or [`RawRecord::default`]; emit via [`write_record`]. For allocation-free
/// read-only scans use [`RecordRef`] via [`RecordReader`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RawRecord {
    bytes: Vec<u8>,
}

raw_field_accessors!(owned RawRecord);

impl From<Vec<u8>> for RawRecord {
    /// Construct a [`RawRecord`] from a pre-validated payload byte vector.
    ///
    /// `bytes` must be the payload of exactly one BAM record (everything after
    /// the 4-byte `block_size` field). No validation is performed.
    fn from(bytes: Vec<u8>) -> Self {
        RawRecord { bytes }
    }
}

impl RawRecord {
    fn set_i32_at(&mut self, off: usize, value: i32) {
        self.bytes[off..off + 4].copy_from_slice(&value.to_le_bytes());
    }

    /// Mutable per-base quality scores. The pileup overlap-removal step zeroes /
    /// scales qualities of overlapping mate bases in place (htslib
    /// `tweak_overlap_quality`), so it needs to write back into the raw payload.
    /// Unlike [`quality_scores`](Self::quality_scores) this does not mask the
    /// `0xff` "missing" sentinel — callers check `sequence_len` and the sentinel
    /// themselves before mutating.
    pub fn quality_scores_mut(&mut self) -> &mut [u8] {
        let base_count = payload_base_count(&self.bytes);
        let start = payload_quality_start(&self.bytes);
        &mut self.bytes[start..start + base_count]
    }

    /// Mutable packed SEQ nibbles. The calmd `use_equal` path rewrites matched
    /// base nibbles to 0 (the `=` code in `seq_nt16`) in place. The packed
    /// representation stores two bases per byte: high nibble = even index, low
    /// nibble = odd index. The returned slice is `(seq_len + 1) / 2` bytes.
    pub fn seq_bytes_mut(&mut self) -> &mut [u8] {
        let base_count = payload_base_count(&self.bytes);
        let start =
            FIXED_HEAD + payload_name_len(&self.bytes) + payload_cigar_op_count(&self.bytes) * 4;
        let end = start + base_count.div_ceil(2);
        &mut self.bytes[start..end]
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

    /// Set pos (offset 4, 0-based). Paired with
    /// [`set_reference_sequence_id`](Self::set_reference_sequence_id) for the
    /// unmapped-read salvage step.
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
        match payload_aux_field_range(&self.bytes, tag) {
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

    /// Copy `other`'s bytes into `self`, reusing `self`'s allocation when possible.
    ///
    /// Unlike `self.clone_from(other)` (which does the same via the derived
    /// `Clone` impl), this method is explicitly named to signal intent: the caller
    /// is recycling a pooled buffer to avoid per-record malloc overhead.
    #[inline]
    pub fn clone_from_raw(&mut self, other: &RawRecord) {
        self.bytes.clear();
        self.bytes.extend_from_slice(&other.bytes);
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

/// A borrowed BAM record's raw payload bytes, valid only until the next read.
///
/// Produced by [`RecordReader`] for allocation-free read-only scans: the bytes
/// are borrowed straight out of the BGZF reader's already-inflated block buffer
/// when the record lies within one block, or out of a reused scratch buffer when
/// it straddles a block boundary. Exposes the same read accessors as
/// [`RawRecord`]; there are no setters (the borrow is read-only).
#[derive(Clone, Copy, Debug)]
pub struct RecordRef<'a> {
    bytes: &'a [u8],
}

raw_field_accessors!(borrowed <'a> RecordRef<'a>);

impl<'a> RecordRef<'a> {
    /// Construct a borrowing view over a raw BAM record payload slice.
    ///
    /// `bytes` must be the payload of exactly one BAM record (everything after
    /// the 4-byte `block_size` field). No validation is performed; callers that
    /// read from a slab guarantee this invariant by construction.
    pub fn from_bytes(bytes: &'a [u8]) -> Self {
        RecordRef { bytes }
    }

    /// The raw payload bytes with the record's full borrow lifetime `'a`, so the
    /// slice (and slices derived from it) outlive a by-value `RecordRef`. Prefer
    /// this over [`as_bytes`](Self::as_bytes) — whose result is tied to `&self` —
    /// when building a borrowing view (e.g. an aux-field iterator) over a record
    /// passed by value.
    pub fn payload(&self) -> &'a [u8] {
        self.bytes
    }
}

/// A borrowing, allocation-free record scanner over a BGZF [`BufRead`] stream.
///
/// [`next`](RecordReader::next) hands out a [`RecordRef`] borrowed straight out
/// of the reader's inflated block buffer whenever the whole record sits within
/// the currently-buffered bytes; a record that straddles a BGZF block boundary
/// is copied once into a reused scratch buffer. The borrow ties each `RecordRef`
/// to `&mut self`, so it must be dropped before the next `next()` — the natural
/// shape for a streaming scan, and the discipline that makes the no-copy borrow
/// sound without `unsafe`.
///
/// The consume of the previous record's bytes is deferred to the start of the
/// next `next()`: BGZF's `consume` only advances the in-block cursor (it does
/// not overwrite the block buffer until the next inflate), so deferring keeps
/// the borrow live for the caller while still advancing correctly.
///
/// Construct via [`RecordReader::new`] over a reader positioned after the header
/// (open with [`open_with_workers`](crate::open_with_workers), call
/// `read_header`, pass `reader.get_mut()`).
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

    /// Read the next record, borrowing its payload. `Ok(None)` at end of records.
    ///
    /// The returned [`RecordRef`] borrows `self`, so the borrow checker forbids
    /// the next call until it is dropped — exactly the invariant the no-copy path
    /// needs.
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
            // Whole record sits in the inflated block buffer: borrow it, no copy.
            // Defer the consume so this borrow stays valid for the caller.
            self.pending_consume = block_size;
            let bytes = &self.reader.fill_buf().map_err(RsomicsError::Io)?[..block_size];
            Ok(Some(RecordRef { bytes }))
        } else {
            // Record straddles a BGZF block boundary: spill into the scratch buffer.
            self.scratch.resize(block_size, 0);
            self.reader
                .read_exact(&mut self.scratch)
                .map_err(RsomicsError::Io)?;
            Ok(Some(RecordRef {
                bytes: &self.scratch,
            }))
        }
    }
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
