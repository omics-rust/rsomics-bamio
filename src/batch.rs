//! Bounded batched framing for raw BAM records.

use std::io::{self, Write};
use std::thread::JoinHandle;

use crossbeam_channel::{Sender, bounded};
use rsomics_common::{Result, RsomicsError};

use crate::raw::RawRecord;
use crate::{ParallelBamWriter, RingBamWriter, finish_ring_bam_writer};

const CHANNEL_DEPTH: usize = 2;

type Batch = Vec<RawRecord>;

/// A BAM writer that frames batches of raw records on a dedicated thread.
pub struct BatchBamWriter {
    sender: Option<Sender<Batch>>,
    handle: Option<JoinHandle<Result<ParallelBamWriter>>>,
}

impl BatchBamWriter {
    /// Wrap a [`ParallelBamWriter`] after its header has been written.
    pub fn new(writer: ParallelBamWriter) -> Self {
        let (sender, receiver) = bounded::<Batch>(CHANNEL_DEPTH);
        let handle = std::thread::spawn(move || frame_loop(writer, receiver));
        BatchBamWriter {
            sender: Some(sender),
            handle: Some(handle),
        }
    }

    /// Enqueue a batch in record order.
    pub fn write_records_batch(&mut self, batch: Vec<RawRecord>) -> Result<()> {
        if batch.is_empty() {
            return Ok(());
        }
        let sender = self
            .sender
            .as_ref()
            .expect("write_records_batch after finish");
        sender.send(batch).ok();
        Ok(())
    }

    /// Flush all queued batches and return the writer.
    pub fn finish(mut self) -> Result<ParallelBamWriter> {
        drop(self.sender.take());
        let handle = self.handle.take().expect("finish called twice");
        handle
            .join()
            .map_err(|_| RsomicsError::Io(io::Error::other("BAM framing thread panicked")))?
    }
}

impl Drop for BatchBamWriter {
    fn drop(&mut self) {
        drop(self.sender.take());
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

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

/// Batched framing backed by a [`RingBamWriter`].
pub struct RingBatchBamWriter {
    sender: Option<Sender<Batch>>,
    handle: Option<JoinHandle<Result<RingBamWriter>>>,
}

impl RingBatchBamWriter {
    /// Wrap a [`RingBamWriter`] (header already written) for batched writes.
    pub fn new(writer: RingBamWriter) -> Self {
        let (sender, receiver) = bounded::<Batch>(CHANNEL_DEPTH);
        let handle = std::thread::spawn(move || ring_frame_loop(writer, receiver));
        RingBatchBamWriter {
            sender: Some(sender),
            handle: Some(handle),
        }
    }

    /// Enqueue a batch in record order.
    pub fn write_records_batch(&mut self, batch: Vec<RawRecord>) -> Result<()> {
        if batch.is_empty() {
            return Ok(());
        }
        let sender = self
            .sender
            .as_ref()
            .expect("write_records_batch after finish");
        sender.send(batch).ok();
        Ok(())
    }

    /// Flush queued batches and finish the BGZF stream.
    pub fn finish(mut self) -> Result<()> {
        drop(self.sender.take());
        let handle = self.handle.take().expect("finish called twice");
        let writer = handle
            .join()
            .map_err(|_| RsomicsError::Io(io::Error::other("BAM framing thread panicked")))??;
        finish_ring_bam_writer(writer)
    }
}

impl Drop for RingBatchBamWriter {
    fn drop(&mut self) {
        drop(self.sender.take());
        if let Some(handle) = self.handle.take()
            && let Ok(Ok(writer)) = handle.join()
        {
            let _ = finish_ring_bam_writer(writer);
        }
    }
}

fn ring_frame_loop(
    mut writer: RingBamWriter,
    receiver: crossbeam_channel::Receiver<Batch>,
) -> Result<RingBamWriter> {
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
