//! Batched BAM write path that lifts the multi-thread write ceiling.
//!
//! [`write_record`](crate::raw::write_record) frames one record at a time onto
//! the [`MultithreadedWriter`](noodles::bgzf::io::MultithreadedWriter): the
//! deflate + CRC + BGZF framing already run on noodles' worker pool, but the
//! caller's thread still pays for every record — two `write_all` calls each
//! going through the `MultithreadedWriter`'s per-write truncation loop, plus the
//! `block_size`-prefix + payload `memcpy` into the writer's staging buffer. On a
//! 4M-record file that is ~8M `write` calls and ~800 MB of staging-buffer copy
//! serialised on the caller's thread, which caps a `par_iter`-parallelised
//! producer (calmd at `-t4`/`-t8` lost to `samtools calmd -@N` purely here).
//!
//! [`BatchBamWriter`] moves that per-record cost off the caller's thread. It
//! owns the `MultithreadedWriter` on a dedicated framing thread and accepts an
//! owned `Vec<RawRecord>` per batch over a bounded channel. The caller's only
//! per-batch cost is moving the `Vec` into the channel (a pointer move, no copy);
//! the framing thread does the `block_size`-prefix + payload concatenation and
//! the `write_all`, after which noodles' deflate pool fans the blocks out exactly
//! as before. The on-disk byte stream is identical to writing the same records in
//! the same order through [`write_record`](crate::raw::write_record): the framing
//! thread emits the same `block_size` (u32 LE) + payload bytes into the same
//! `MultithreadedWriter`, whose BGZF block boundaries depend only on the byte
//! sequence, not on which thread fed it.
//!
//! The channel is bounded so a fast producer cannot outrun the framing thread
//! and grow memory without limit; a full channel back-pressures the producer,
//! which is the desired pipeline behaviour. [`finish`](BatchBamWriter::finish)
//! drains the channel, flushes the writer (appending the BGZF EOF block), and
//! propagates any write error that occurred on the framing thread.

use std::io::Write;
use std::thread::JoinHandle;

use crossbeam_channel::{Sender, bounded};
use rsomics_common::{Result, RsomicsError};

use crate::ParallelBamWriter;
use crate::raw::RawRecord;

/// Bounded depth of the batch channel. Two in-flight batches let the framing
/// thread work on one while the producer assembles the next, without letting an
/// unbounded queue of un-framed batches accumulate in memory.
const CHANNEL_DEPTH: usize = 2;

/// A batch of records handed to the framing thread, kept owned so the producer's
/// per-batch cost is a `Vec` move rather than a copy.
type Batch = Vec<RawRecord>;

/// A BAM writer that frames batches of raw records on a dedicated thread.
///
/// Construct from a [`ParallelBamWriter`] whose header has already been written
/// (call [`ParallelBamWriter::write_alignment_header`] /
/// `write_header` first). Hand whole batches to
/// [`write_records_batch`](Self::write_records_batch); call
/// [`finish`](Self::finish) at the end to flush and surface any error.
///
/// This is opt-in and additive: it does not change
/// [`write_record`](crate::raw::write_record), which the rest of the
/// `rsomics-bam-*` family keeps using unchanged. The two paths produce identical
/// bytes for identical record sequences.
pub struct BatchBamWriter {
    sender: Option<Sender<Batch>>,
    handle: Option<JoinHandle<Result<ParallelBamWriter>>>,
}

impl BatchBamWriter {
    /// Wrap a [`ParallelBamWriter`] (header already written) for batched writes.
    ///
    /// The writer is moved onto a framing thread that drains batches off a
    /// bounded channel, concatenates each record's `block_size` prefix + payload,
    /// and writes the concatenation into the underlying `MultithreadedWriter`.
    pub fn new(writer: ParallelBamWriter) -> Self {
        let (sender, receiver) = bounded::<Batch>(CHANNEL_DEPTH);
        let handle = std::thread::spawn(move || frame_loop(writer, receiver));
        BatchBamWriter {
            sender: Some(sender),
            handle: Some(handle),
        }
    }

    /// Hand a batch of records to the framing thread, in record order.
    ///
    /// Records are written in the order they appear in `batch`, and successive
    /// `write_records_batch` calls concatenate in call order, so the byte stream
    /// matches [`write_record`](crate::raw::write_record) over the same records.
    /// The `Vec` is moved into the channel (no copy). A full channel blocks until
    /// the framing thread drains a batch (intended back-pressure). If the framing
    /// thread has already failed, the error is surfaced from
    /// [`finish`](Self::finish), not here — this returns `Ok` for the enqueue.
    pub fn write_records_batch(&mut self, batch: Vec<RawRecord>) -> Result<()> {
        if batch.is_empty() {
            return Ok(());
        }
        let sender = self
            .sender
            .as_ref()
            .expect("write_records_batch after finish");
        // A send error means the framing thread died (and dropped the receiver);
        // its real error is recovered by finish() joining the thread.
        sender.send(batch).ok();
        Ok(())
    }

    /// Flush all queued batches, append the BGZF EOF block, and return the
    /// underlying writer. Any I/O error raised on the framing thread surfaces
    /// here. Must be called exactly once; dropping without `finish` still flushes
    /// (via the writer's own `Drop`) but discards the error.
    pub fn finish(mut self) -> Result<ParallelBamWriter> {
        // Dropping the sender closes the channel so the framing thread's recv
        // loop ends and it returns the writer (or its error).
        drop(self.sender.take());
        let handle = self.handle.take().expect("finish called twice");
        handle.join().expect("framing thread panicked")
    }
}

impl Drop for BatchBamWriter {
    fn drop(&mut self) {
        // If finish() was not called, still shut the framing thread down cleanly
        // so the writer flushes and the BGZF EOF block is appended on its Drop.
        // The error is unavailable here — callers that need it must use finish().
        drop(self.sender.take());
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

/// Framing-thread body: drain batches, concatenate each record's `block_size`
/// (u32 LE) + payload into one reused buffer, and write the buffer into the
/// `MultithreadedWriter`. Concatenating a whole batch into one `write_all`
/// collapses the per-record write-call overhead; the writer's pool still deflates
/// and frames the BGZF blocks. Returns the writer so `finish` can flush it on the
/// caller's behalf and recover any error.
fn frame_loop(
    mut writer: ParallelBamWriter,
    receiver: crossbeam_channel::Receiver<Batch>,
) -> Result<ParallelBamWriter> {
    let inner = writer.get_mut();
    let mut buf: Vec<u8> = Vec::new();
    while let Ok(batch) = receiver.recv() {
        buf.clear();
        for record in &batch {
            let payload = record.as_bytes();
            let block_size = u32::try_from(payload.len())
                .map_err(|e| RsomicsError::InvalidInput(format!("record too large: {e}")))?;
            buf.extend_from_slice(&block_size.to_le_bytes());
            buf.extend_from_slice(payload);
        }
        inner.write_all(&buf).map_err(RsomicsError::Io)?;
    }
    Ok(writer)
}
