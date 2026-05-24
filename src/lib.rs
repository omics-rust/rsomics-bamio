//! Shared BAM input/output for the `rsomics-bam-*` tools.
//!
//! `samtools` inflates BGZF on a single thread by default, and a single-threaded
//! pure-Rust reader (zlib-rs) loses to its libdeflate inner loop. The lever that
//! wins is libdeflate on every block (enabled here) plus the right threading
//! shape for each direction:
//!
//! - **Reader**: at one worker, noodles' `MultithreadedReader` is pure overhead
//!   — a crossbeam channel per block and a cross-thread hand-off at `bounded(1)`,
//!   with no main-thread work to overlap the lone inflate worker against. So
//!   `workers == 1` routes through the plain `bgzf::io::Reader`, which beats
//!   samtools single-threaded on libdeflate alone; `workers >= 2` uses the pool.
//!   Both inner readers are erased behind `Box<dyn Read + Send>`, so the public
//!   `ParallelBamReader` stays one concrete type. The dyn dispatch is one vtable
//!   hop per BGZF block (~64 KiB), not per byte — negligible against inflate.
//! - **Writer**: the worker-pool writer wins even at one worker, because the
//!   deflate thread overlaps the caller's read+decode on the main thread (BGZF
//!   compression dominates a write-out). So file output always uses the pool.
//!
//! The reader is forward-only: BAM tools that need index-driven seeks (region
//! queries, `bedcov`) build their own seekable `bam::io::Reader<File>` directly,
//! which was never this primitive's job. So `Box<dyn Read + Send>` (no `Seek`)
//! is sufficient for every consumer of this crate.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::num::NonZero;
use std::path::Path;

use noodles::{bam, bgzf};
use rsomics_common::{Result, RsomicsError};

pub mod batch;
pub mod par_writer;
pub mod raw;

pub use batch::{BatchBamWriter, WsBatchBamWriter};
pub use par_writer::{WorkStealingBgzfWriter, create_ws_bgzf};

/// Buffer the file under the single-threaded BGZF reader so the per-block frame
/// reads coalesce into one ~256 KiB refill instead of two syscalls per block.
const READ_BUFFER: usize = 256 * 1024;

/// A BAM reader over a BGZF stream, forward-only.
///
/// The inner BGZF reader is type-erased: at `workers == 1` it is the plain
/// single-threaded reader (libdeflate, no channel overhead); at `workers >= 2`
/// it inflates blocks across a worker pool. Both are boxed as
/// `Box<dyn BufRead + Send>` — `BufRead` (not just `Read`) so the zero-copy
/// [`raw::RecordReader`] can borrow record bytes straight out of the inflated
/// block buffer via `reader.get_mut()`. `BufRead: Read`, so this stays one
/// concrete type and every existing `Read`-based consumer (`read_record`,
/// `bam::io::Reader<R>`) keeps working unchanged.
pub type ParallelBamReader = bam::io::Reader<Box<dyn BufRead + Send>>;

/// A BAM writer whose BGZF blocks are deflated across a worker pool. Even at one
/// worker the pool wins, because the deflate thread overlaps the caller's
/// read+decode (see [`create_with_workers`]).
pub type ParallelBamWriter = bam::io::Writer<bgzf::io::MultithreadedWriter<File>>;

/// A BAM writer backed by the work-stealing BGZF writer. At `workers >= 2` this
/// eliminates noodles' per-block channel allocation and matches samtools'
/// fixed-ring bgzf_mt design. Use [`create_ws_bam_writer`] to construct.
pub type WsBamWriter = bam::io::Writer<WorkStealingBgzfWriter<File>>;

/// Open `input` with one inflate worker per available core.
pub fn open_parallel(input: &Path) -> Result<ParallelBamReader> {
    let workers = std::thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN);
    open_with_workers(input, workers)
}

/// Open `input` with an explicit worker count. `workers == 1` uses the plain
/// single-threaded reader (no channel overhead); `workers >= 2` inflates blocks
/// across a worker pool. The inner reader is boxed as `Box<dyn Read + Send>`, so
/// both cases share one return type.
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

/// Create `output` with one deflate worker per available core. The BGZF EOF
/// block is appended when the returned writer is dropped.
pub fn create_parallel(output: &Path) -> Result<ParallelBamWriter> {
    let workers = std::thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN);
    create_with_workers(output, workers)
}

/// Create `output` with an explicit deflate worker count. Even at `workers == 1`
/// this uses the worker-pool writer: deflate runs on the worker thread while the
/// caller's read+decode runs on the main thread, so the two pipeline-overlap.
/// (The reader's worker pool, by contrast, is pure overhead at one worker — see
/// [`open_with_workers`] — because there is no main-thread work to overlap with
/// a single inflate worker.) The BGZF EOF block is appended on drop.
pub fn create_with_workers(output: &Path, workers: NonZero<usize>) -> Result<ParallelBamWriter> {
    let file = File::create(output)
        .map_err(|e| RsomicsError::InvalidInput(format!("creating {}: {e}", output.display())))?;
    Ok(bam::io::Writer::from(
        bgzf::io::MultithreadedWriter::with_worker_count(workers, file),
    ))
}

/// Create `output` backed by the work-stealing BGZF writer with `workers`
/// deflate threads. At `workers >= 2` this eliminates noodles' per-block channel
/// allocation, matching samtools' fixed-ring `bgzf_mt` design. The BGZF EOF
/// block is appended when the writer is finished or dropped.
pub fn create_ws_with_workers(output: &Path, workers: NonZero<usize>) -> Result<WsBamWriter> {
    let file = File::create(output)
        .map_err(|e| RsomicsError::InvalidInput(format!("creating {}: {e}", output.display())))?;
    Ok(bam::io::Writer::from(WorkStealingBgzfWriter::new(
        file, workers,
    )))
}

/// Finish a [`WsBamWriter`], flushing pending BGZF blocks and appending the EOF
/// marker. This is a free function rather than a method because the `WsBamWriter`
/// type alias wraps `bam::io::Writer<WorkStealingBgzfWriter<File>>` and the
/// inner `WorkStealingBgzfWriter::finish` method consumes the inner writer.
pub fn finish_ws_bam_writer(writer: WsBamWriter) -> Result<()> {
    writer.into_inner().finish().map_err(RsomicsError::Io)
}
