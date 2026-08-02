//! Validated raw BAM records and bounded parallel BGZF I/O.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::num::NonZero;
use std::path::Path;

use noodles::{bam, bgzf};
use rsomics_common::{Result, RsomicsError};

pub mod batch;
mod indexed;
pub mod raw;
pub mod ring_writer;

pub use batch::{BatchBamWriter, RingBatchBamWriter};
pub use indexed::{IndexedAlignmentReader, open_indexed_alignment, visit_raw_alignment_records};
pub use ring_writer::{RingBgzfWriter, create_ring_bgzf};

const READ_BUFFER: usize = 256 * 1024;

/// A BAM reader over a BGZF stream, forward-only.
pub type ParallelBamReader = bam::io::Reader<Box<dyn BufRead + Send>>;

/// A BAM writer whose BGZF blocks are deflated across a worker pool.
pub type ParallelBamWriter = bam::io::Writer<bgzf::io::MultithreadedWriter<File>>;

/// A BAM writer backed by reusable BGZF compression buffers.
pub type RingBamWriter = bam::io::Writer<RingBgzfWriter<File>>;

/// Open `input` with one inflate worker per available core.
pub fn open_parallel(input: &Path) -> Result<ParallelBamReader> {
    let workers = std::thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN);
    open_with_workers(input, workers)
}

/// Open `input` with an explicit inflate-worker count.
pub fn open_with_workers(input: &Path, workers: NonZero<usize>) -> Result<ParallelBamReader> {
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let inner: Box<dyn BufRead + Send> = if workers.get() == 1 {
        let buffered = BufReader::with_capacity(READ_BUFFER, file);
        Box::new(bgzf::io::Reader::new(buffered))
    } else {
        Box::new(bgzf::io::MultithreadedReader::with_worker_count(
            workers, file,
        ))
    };
    Ok(bam::io::Reader::from(inner))
}

/// Create `output` with one deflate worker per available core.
pub fn create_parallel(output: &Path) -> Result<ParallelBamWriter> {
    let workers = std::thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN);
    create_with_workers(output, workers)
}

/// Create `output` with an explicit deflate-worker count.
pub fn create_with_workers(output: &Path, workers: NonZero<usize>) -> Result<ParallelBamWriter> {
    let file = File::create(output)
        .map_err(|e| RsomicsError::InvalidInput(format!("creating {}: {e}", output.display())))?;
    Ok(bam::io::Writer::from(
        bgzf::io::MultithreadedWriter::with_worker_count(workers, file),
    ))
}

/// Create `output` with reusable BGZF buffers and `workers` deflate threads.
pub fn create_ring_with_workers(output: &Path, workers: NonZero<usize>) -> Result<RingBamWriter> {
    let file = File::create(output)
        .map_err(|e| RsomicsError::InvalidInput(format!("creating {}: {e}", output.display())))?;
    Ok(bam::io::Writer::from(RingBgzfWriter::new(file, workers)))
}

/// Finish a [`RingBamWriter`] and append the BGZF EOF marker.
pub fn finish_ring_bam_writer(writer: RingBamWriter) -> Result<()> {
    writer
        .into_inner()
        .finish()
        .map(drop)
        .map_err(RsomicsError::Io)
}
