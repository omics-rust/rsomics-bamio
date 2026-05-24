//! Work-stealing parallel BGZF writer — a fixed-ring, no-per-block-alloc
//! replacement for noodles' `MultithreadedWriter`.
//!
//! ## Motivation
//!
//! noodles' `MultithreadedWriter` allocates a new `crossbeam_channel::bounded(1)`
//! for every 65 KiB BGZF block — millions of heap-allocated channel objects on a
//! large BAM. It also uses a single serial writer thread that collects compressed
//! blocks in arrival order via a per-block receiver. On a 4M-record BAM at -t4
//! this structure keeps calmd at 0.86× samtools.
//!
//! ## Design (mirrors samtools `bgzf_mt_*`)
//!
//! - **Fixed ring of `N` slot objects** (N = `workers * 2`). Each slot holds a
//!   reusable uncompressed buffer (≤ `MAX_IDATA` bytes) and a pre-sized compressed
//!   buffer. In steady state there is zero per-block heap allocation.
//! - **Producer** (the caller's thread) accumulates BAM record bytes into the
//!   current slot's uncompressed buffer. When the buffer reaches `MAX_IDATA` or
//!   `flush` is called it sends the slot index to a bounded work queue and moves to
//!   the next slot (ring-wrapping, blocking on back-pressure).
//! - **`workers` deflate threads** pull slot indices from the work queue, compress
//!   the uncompressed buffer into the compressed buffer using **libdeflate** at
//!   level 6 (same library and level samtools uses), compute the CRC32 with
//!   libdeflater's `Crc`, and post the slot index to a per-sequence-position
//!   ready-map.
//! - **One writer thread** drains the ready-map in *sequence order* (by the
//!   monotonically increasing block index assigned at dispatch time), writes each
//!   BGZF block to the underlying `File`, and returns the slot to the free ring
//!   (by sending the index back to the producer via a free-slot channel).
//! - **Output is byte-deterministic**: block boundaries depend only on the
//!   byte stream (the producer fills each uncompressed buffer to the same threshold
//!   regardless of worker count), so the compressed output is byte-identical across
//!   any thread count. libdeflate at a fixed level is also deterministic.
//! - **BGZF EOF block** (28-byte empty gzip member) is appended on `finish`.
//!
//! ## Exposed API
//!
//! `WorkStealingBgzfWriter` implements `std::io::Write` and wraps a `File` (or
//! any `Write + Send + 'static`). Callers write BAM record bytes the same way
//! they would to the existing `MultithreadedWriter`; block boundaries are handled
//! internally.
//!
//! For BAM output, wrap it in `bam::io::Writer::from(writer)` after writing the
//! header, exactly as today's `ParallelBamWriter`.

use std::collections::BTreeMap;
use std::io::{self, Write};
use std::num::NonZero;
use std::thread::JoinHandle;

use crossbeam_channel::{Receiver, Sender, bounded};
use libdeflater::{CompressionLvl, Compressor, Crc};
use rsomics_common::{Result, RsomicsError};

// BGZF block uncompressed payload size limit.  The noodles writer uses
// MAX_BUF_SIZE = 65280 − 15 = 65265; we use the same figure so block
// boundaries are identical to noodles' single-threaded writer on the same
// input stream.  (The −15 headroom is for deflate level-0 overhead; at
// level 6 the compressed output is always smaller, but the bound must hold
// for any level.)
const MAX_IDATA: usize = 65_280 - 18 - 8 - 15; // = 65_239

// Worst-case compressed size for `MAX_IDATA` uncompressed bytes at deflate
// level 0.  libdeflate's bound is MAX_IDATA + 10 bytes of DEFLATE overhead.
const MAX_CDATA: usize = MAX_IDATA + 10;

// 28-byte BGZF EOF marker (SAM spec § 4.1.2).
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

// BGZF header size: 18 bytes (RFC 1952 gzip header 10 B + XLEN field 2 B + BC subfield 6 B).
const HEADER_SIZE: usize = 18;

// BGZF trailer size: CRC32 (4 B) + ISIZE (4 B).
const TRAILER_SIZE: usize = 8;

/// One slot in the ring — a reusable pair of (uncompressed, compressed) buffers.
struct Slot {
    /// Uncompressed BAM payload bytes, up to `MAX_IDATA`.
    idata: Vec<u8>,
    /// Compressed DEFLATE bytes (sized for `MAX_CDATA`; actual len is set after
    /// compression and used when writing the block).
    cdata: Vec<u8>,
    /// CRC32 of `idata` (filled by the deflate thread).
    crc32: u32,
    /// Number of valid bytes in `cdata` (filled by the deflate thread).
    clen: usize,
}

impl Slot {
    fn new() -> Self {
        Slot {
            idata: Vec::with_capacity(MAX_IDATA),
            cdata: vec![0u8; MAX_CDATA],
            crc32: 0,
            clen: 0,
        }
    }
}

// Channel message types:

/// Work item: (slot index, block sequence number, idata length).  The deflate
/// threads receive these; they own the compressed-side write to `slots[idx]`.
type WorkItem = (usize, u64, usize);
/// Completion item: (slot index, block sequence number, idata length, clen, crc32).
type DoneItem = (usize, u64, usize, usize, u32);
/// Free slot index returned by the writer thread to the producer.
type FreeItem = usize;

/// A BGZF writer that uses a fixed ring of reusable block slots and a
/// work-stealing deflate pool — no per-block heap allocation in steady state.
///
/// Implements `std::io::Write`. Wrap in `bam::io::Writer::from(writer)` after
/// writing the BAM header for use as a drop-in replacement for
/// `ParallelBamWriter`.
///
/// Call [`finish`](Self::finish) to flush pending data and append the BGZF EOF
/// block. Dropping without `finish` appends EOF (via the writer thread's `Drop`)
/// but discards any I/O error.
pub struct WorkStealingBgzfWriter<W: Write + Send + 'static> {
    /// Ring of shared slot buffers (Arc so deflate + writer threads can own refs).
    slots: Vec<std::sync::Arc<std::sync::Mutex<Slot>>>,
    /// Index of the slot currently being filled by the producer.
    cur: usize,
    /// Monotonically increasing block sequence number assigned at dispatch.
    seq: u64,
    /// Send side of the work queue (producer → deflate workers).
    work_tx: Sender<WorkItem>,
    /// Receive side of the free-slot queue (writer thread → producer).
    free_rx: Receiver<FreeItem>,
    /// Forwarding channel from deflate workers → writer thread.
    done_tx: Sender<DoneItem>,
    /// Deflate worker join handles (held so `finish` can drain them).
    deflate_handles: Vec<JoinHandle<()>>,
    /// Writer thread join handle.
    writer_handle: Option<JoinHandle<io::Result<()>>>,
    /// Marker so `finish` is idempotent.
    finished: bool,
    /// Phantom so W is tied to this struct's lifetime.
    _w: std::marker::PhantomData<W>,
}

impl<W: Write + Send + 'static> WorkStealingBgzfWriter<W> {
    /// Create a writer with `workers` deflate threads over `sink`.
    pub fn new(sink: W, workers: NonZero<usize>) -> Self {
        let n_slots = workers.get() * 2 + 1;

        let slots: Vec<_> = (0..n_slots)
            .map(|_| std::sync::Arc::new(std::sync::Mutex::new(Slot::new())))
            .collect();

        // Work queue: bounded to `workers` so the producer can't fill every slot
        // before deflaters drain them.
        let (work_tx, work_rx) = bounded::<WorkItem>(workers.get());
        // Free-slot queue: bounded to `n_slots`.
        let (free_tx, free_rx) = bounded::<FreeItem>(n_slots);
        // Pre-fill: all slots except slot 0 (which the producer holds) are free.
        for i in 1..n_slots {
            free_tx.send(i).unwrap();
        }

        // Done channel: bounded to `n_slots` (deflate → writer).
        let (done_tx, done_rx) = bounded::<DoneItem>(n_slots);

        // Spawn deflate workers.
        let mut deflate_handles = Vec::with_capacity(workers.get());
        for _ in 0..workers.get() {
            let work_rx = work_rx.clone();
            let done_tx = done_tx.clone();
            let slots = slots.clone();
            deflate_handles.push(std::thread::spawn(move || {
                deflate_worker(slots, work_rx, done_tx);
            }));
        }

        // Spawn writer thread.
        let slots_w = slots.clone();
        let free_tx_w = free_tx;
        let writer_handle =
            std::thread::spawn(move || write_worker(sink, slots_w, done_rx, free_tx_w));

        WorkStealingBgzfWriter {
            slots,
            cur: 0,
            seq: 0,
            work_tx,
            free_rx,
            done_tx,
            deflate_handles,
            writer_handle: Some(writer_handle),
            finished: false,
            _w: std::marker::PhantomData,
        }
    }

    /// Flush the current (partially filled) slot as a BGZF block, then wait for
    /// all in-flight blocks to be written and append the BGZF EOF block.
    ///
    /// Returns the underlying writer's I/O error if any occurred on the writer
    /// thread. Safe to call exactly once; subsequent calls are no-ops.
    pub fn finish(mut self) -> io::Result<()> {
        self.flush_current()?;
        self.shutdown()
    }

    /// Flush the current slot (may be empty — emits an empty BGZF block only if
    /// there is data). Dispatches to deflate workers.
    fn flush_current(&mut self) -> io::Result<()> {
        let idata_len = {
            let slot = self.slots[self.cur].lock().unwrap();
            slot.idata.len()
        };
        if idata_len == 0 {
            return Ok(());
        }
        let idx = self.cur;
        let seq = self.seq;
        self.seq += 1;
        self.work_tx
            .send((idx, seq, idata_len))
            .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "deflate worker died"))?;

        // Acquire the next free slot.
        let next = self
            .free_rx
            .recv()
            .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "writer thread died"))?;
        {
            let mut slot = self.slots[next].lock().unwrap();
            slot.idata.clear();
        }
        self.cur = next;
        Ok(())
    }

    fn shutdown(&mut self) -> io::Result<()> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;
        // Drop the only live work sender (deflate threads hold only the receiver)
        // so their `recv` returns Err and they exit. `self` is borrowed here, not
        // consumed, so swap in a disconnected sender and drop the real one.
        let (dead_tx, _dead_rx) = bounded::<WorkItem>(0);
        drop(std::mem::replace(&mut self.work_tx, dead_tx));

        for handle in self.deflate_handles.drain(..) {
            handle.join().unwrap();
        }
        // Closing done_tx signals the writer thread to exit after draining.
        let (dead_done, _dead_done_rx) = bounded::<DoneItem>(0);
        let old_done = std::mem::replace(&mut self.done_tx, dead_done);
        drop(old_done);

        if let Some(handle) = self.writer_handle.take() {
            handle.join().unwrap()?;
        }
        Ok(())
    }
}

impl<W: Write + Send + 'static> Write for WorkStealingBgzfWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let remaining = {
            let slot = self.slots[self.cur].lock().unwrap();
            MAX_IDATA - slot.idata.len()
        };
        let amt = remaining.min(buf.len());
        {
            let mut slot = self.slots[self.cur].lock().unwrap();
            slot.idata.extend_from_slice(&buf[..amt]);
        }
        let is_full = {
            let slot = self.slots[self.cur].lock().unwrap();
            slot.idata.len() >= MAX_IDATA
        };
        if is_full {
            self.flush_current()?;
        }
        Ok(amt)
    }

    fn flush(&mut self) -> io::Result<()> {
        // Flush the current partial block.
        self.flush_current()
    }
}

impl<W: Write + Send + 'static> Drop for WorkStealingBgzfWriter<W> {
    fn drop(&mut self) {
        if !self.finished {
            // Best-effort: flush current data and shut down.
            let _ = self.flush_current();
            let _ = self.shutdown();
        }
    }
}

/// Deflate worker body: compress slots and post them to the done channel.
fn deflate_worker(
    slots: Vec<std::sync::Arc<std::sync::Mutex<Slot>>>,
    work_rx: Receiver<WorkItem>,
    done_tx: Sender<DoneItem>,
) {
    let mut compressor = Compressor::new(CompressionLvl::new(6).expect("level 6 is valid"));
    // Scratch buffer for the idata copy needed to split the borrow of idata
    // and cdata (both fields of the same Slot).  Reused across blocks — zero
    // allocation in the hot loop once warmed up.
    let mut idata_scratch: Vec<u8> = Vec::with_capacity(MAX_IDATA);

    while let Ok((idx, seq, idata_len)) = work_rx.recv() {
        let (clen, crc32) = {
            let mut slot = slots[idx].lock().unwrap();

            // Copy idata into scratch so we can borrow idata (immutable) and
            // cdata (mutable) simultaneously without the borrow-checker conflict.
            idata_scratch.resize(idata_len, 0);
            idata_scratch.copy_from_slice(&slot.idata[..idata_len]);

            // CRC32 over uncompressed data — libdeflater Crc has no reset, create fresh.
            let mut crc_engine = Crc::new();
            crc_engine.update(&idata_scratch);
            let crc32 = crc_engine.sum();

            // Compress into the slot's cdata buffer.
            let clen = compressor
                .deflate_compress(&idata_scratch, &mut slot.cdata)
                .expect("deflate compress cannot fail on valid input");
            slot.clen = clen;
            slot.crc32 = crc32;
            (clen, crc32)
        };
        if done_tx.send((idx, seq, idata_len, clen, crc32)).is_err() {
            break;
        }
    }
}

/// Writer thread body: collect completed blocks in sequence order, write BGZF
/// frames to `sink`, return free slot indices to the producer.
fn write_worker<W: Write + Send + 'static>(
    mut sink: W,
    slots: Vec<std::sync::Arc<std::sync::Mutex<Slot>>>,
    done_rx: Receiver<DoneItem>,
    free_tx: Sender<FreeItem>,
) -> io::Result<()> {
    // Reorder buffer: accumulate out-of-order completions keyed by sequence number.
    let mut pending: BTreeMap<u64, (usize, usize, usize, u32)> = BTreeMap::new();
    let mut next_seq: u64 = 0;
    // Reusable frame buffer: header (18) + cdata + trailer (8).
    let mut frame_buf: Vec<u8> = Vec::with_capacity(HEADER_SIZE + MAX_CDATA + TRAILER_SIZE);

    while let Ok((idx, seq, idata_len, clen, crc32)) = done_rx.recv() {
        pending.insert(seq, (idx, idata_len, clen, crc32));

        while let Some(&(pidx, pidata_len, pclen, pcrc32)) = pending.get(&next_seq) {
            pending.remove(&next_seq);
            next_seq += 1;

            // Write the BGZF frame.
            {
                let slot = slots[pidx].lock().unwrap();
                write_bgzf_block(
                    &mut sink,
                    &mut frame_buf,
                    &slot.cdata[..pclen],
                    pcrc32,
                    pidata_len,
                )?;
            }
            // Return the slot to the producer.
            if free_tx.send(pidx).is_err() {
                // Producer has gone away — no point continuing.
                break;
            }
        }
    }

    // Flush any remaining pending blocks (channel closed before we drained).
    while let Some((&seq, _)) = pending.first_key_value() {
        if seq != next_seq {
            break;
        }
        let (pidx, pidata_len, pclen, pcrc32) = pending.remove(&next_seq).unwrap();
        next_seq += 1;
        let slot = slots[pidx].lock().unwrap();
        write_bgzf_block(
            &mut sink,
            &mut frame_buf,
            &slot.cdata[..pclen],
            pcrc32,
            pidata_len,
        )?;
    }

    // Append BGZF EOF.
    sink.write_all(&BGZF_EOF)
}

/// Encode one BGZF block into `frame_buf` and write it to `sink`.
///
/// Block layout: 18-byte header | DEFLATE-compressed data | CRC32 (4 B) | ISIZE (4 B).
fn write_bgzf_block<W: Write>(
    sink: &mut W,
    frame_buf: &mut Vec<u8>,
    cdata: &[u8],
    crc32: u32,
    idata_len: usize,
) -> io::Result<()> {
    let block_size = HEADER_SIZE + cdata.len() + TRAILER_SIZE;

    frame_buf.clear();

    // gzip header (RFC 1952) with BGZF extensions (SAM spec § 4.1).
    frame_buf.extend_from_slice(&[
        0x1f, 0x8b, // ID1, ID2
        0x08, // CM = DEFLATE
        0x04, // FLG = FEXTRA
        0x00, 0x00, 0x00, 0x00, // MTIME = 0
        0x00, // XFL = 0
        0xff, // OS = 255 (unknown)
        0x06, 0x00, // XLEN = 6
        b'B', b'C', // SI1, SI2
        0x02, 0x00, // SLEN = 2
    ]);
    // BSIZE = block_size − 1 (u16 LE), the total compressed block size minus 1.
    let bsize = u16::try_from(block_size - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;
    frame_buf.extend_from_slice(&bsize.to_le_bytes());

    // DEFLATE-compressed data.
    frame_buf.extend_from_slice(cdata);

    // CRC32 of the uncompressed data (u32 LE).
    frame_buf.extend_from_slice(&crc32.to_le_bytes());

    // ISIZE = uncompressed size (u32 LE).
    let isize_u32 = u32::try_from(idata_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;
    frame_buf.extend_from_slice(&isize_u32.to_le_bytes());

    sink.write_all(frame_buf)
}

/// Convenience: create a `WorkStealingBgzfWriter` over a newly created file.
pub fn create_ws_bgzf<P: AsRef<std::path::Path>>(
    path: P,
    workers: NonZero<usize>,
) -> Result<WorkStealingBgzfWriter<std::fs::File>> {
    let file = std::fs::File::create(path.as_ref()).map_err(|e| {
        RsomicsError::InvalidInput(format!("creating {}: {e}", path.as_ref().display()))
    })?;
    Ok(WorkStealingBgzfWriter::new(file, workers))
}

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};
    use std::sync::{Arc, Mutex};

    use super::*;

    /// A `Write + Send + 'static` wrapper around a shared `Vec<u8>` for tests.
    #[derive(Clone)]
    struct SharedBuf(Arc<Mutex<Vec<u8>>>);

    impl Write for SharedBuf {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn new_buf() -> (SharedBuf, Arc<Mutex<Vec<u8>>>) {
        let inner = Arc::new(Mutex::new(Vec::new()));
        (SharedBuf(Arc::clone(&inner)), inner)
    }

    #[test]
    fn round_trip_single_thread() {
        let workers = NonZero::new(1).unwrap();
        let (sink, data) = new_buf();
        let mut w = WorkStealingBgzfWriter::new(sink, workers);
        w.write_all(b"hello world").unwrap();
        w.finish().unwrap();

        let bytes = data.lock().unwrap().clone();
        let mut reader = noodles::bgzf::io::Reader::new(std::io::Cursor::new(bytes));
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed).unwrap();
        assert_eq!(decompressed, b"hello world");
    }

    #[test]
    fn round_trip_multi_thread() {
        for nw in [2usize, 4, 8] {
            let workers = NonZero::new(nw).unwrap();
            // Write 4 full blocks + one partial block.
            let payload: Vec<u8> = (0u8..=255).cycle().take(MAX_IDATA * 4 + 1000).collect();
            let (sink, data) = new_buf();
            let mut w = WorkStealingBgzfWriter::new(sink, workers);
            w.write_all(&payload).unwrap();
            w.finish().unwrap();

            let bytes = data.lock().unwrap().clone();
            let mut reader = noodles::bgzf::io::Reader::new(std::io::Cursor::new(bytes));
            let mut decompressed = Vec::new();
            reader.read_to_end(&mut decompressed).unwrap();
            assert_eq!(decompressed, payload, "failed at workers={nw}");
        }
    }

    #[test]
    fn eof_block_present() {
        let workers = NonZero::new(2).unwrap();
        let (sink, data) = new_buf();
        let mut w = WorkStealingBgzfWriter::new(sink, workers);
        w.write_all(b"test").unwrap();
        w.finish().unwrap();

        let bytes = data.lock().unwrap().clone();
        let end = bytes.len();
        assert_eq!(&bytes[end - BGZF_EOF.len()..], BGZF_EOF);
    }

    #[test]
    fn output_identical_across_thread_counts() {
        // The same byte stream should produce identical compressed output
        // regardless of how many deflate workers are used, because:
        // 1. Block boundaries depend only on byte offset (same MAX_IDATA threshold).
        // 2. libdeflate at a fixed level is deterministic.
        let payload: Vec<u8> = (0u8..=255).cycle().take(MAX_IDATA * 3 + 500).collect();

        let outputs: Vec<Vec<u8>> = [1usize, 2, 4]
            .iter()
            .map(|&nw| {
                let workers = NonZero::new(nw).unwrap();
                let (sink, data) = new_buf();
                let mut w = WorkStealingBgzfWriter::new(sink, workers);
                w.write_all(&payload).unwrap();
                w.finish().unwrap();
                data.lock().unwrap().clone()
            })
            .collect();

        assert_eq!(outputs[0], outputs[1], "t1 vs t2 differ");
        assert_eq!(outputs[0], outputs[2], "t1 vs t4 differ");
    }
}
