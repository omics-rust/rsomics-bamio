//! The batched writer must produce byte-identical output to the per-record
//! [`write_record`](rsomics_bamio::raw::write_record) path for the same record
//! sequence — the contract that lets a tool opt into batching without changing
//! its on-disk result. Verified at multiple worker counts and over a file large
//! enough to span many BGZF blocks.

use std::fs;
use std::num::NonZero;
use std::path::Path;

use noodles::sam::alignment::io::Write as _;
use rsomics_bamio::BatchBamWriter;
use rsomics_bamio::raw::{self, RawRecord};

fn fixture() -> &'static Path {
    Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/golden/small.bam"
    ))
}

fn read_all(path: &Path) -> (Vec<RawRecord>, noodles::sam::Header) {
    let mut reader = rsomics_bamio::open_with_workers(path, NonZero::new(1).unwrap()).unwrap();
    let header = reader.read_header().unwrap();
    let mut records = Vec::new();
    let mut rec = RawRecord::default();
    while raw::read_record(reader.get_mut(), &mut rec).unwrap() != 0 {
        records.push(rec.clone());
    }
    (records, header)
}

/// Write `records` through `write_record`, one at a time, returning the file bytes.
fn write_per_record(
    records: &[RawRecord],
    header: &noodles::sam::Header,
    workers: usize,
) -> Vec<u8> {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    {
        let mut writer =
            rsomics_bamio::create_with_workers(tmp.path(), NonZero::new(workers).unwrap()).unwrap();
        writer.write_alignment_header(header).unwrap();
        for r in records {
            raw::write_record(writer.get_mut(), r).unwrap();
        }
    }
    fs::read(tmp.path()).unwrap()
}

/// Write `records` through the batch writer in `batch_size`-record chunks,
/// returning the file bytes.
fn write_batched(
    records: &[RawRecord],
    header: &noodles::sam::Header,
    workers: usize,
    batch_size: usize,
) -> Vec<u8> {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    {
        let mut writer =
            rsomics_bamio::create_with_workers(tmp.path(), NonZero::new(workers).unwrap()).unwrap();
        writer.write_alignment_header(header).unwrap();
        let mut batch_writer = BatchBamWriter::new(writer);
        for chunk in records.chunks(batch_size) {
            batch_writer.write_records_batch(chunk.to_vec()).unwrap();
        }
        // finish flushes the framing thread, appends the BGZF EOF block, and
        // surfaces any write error.
        batch_writer.finish().unwrap();
    }
    fs::read(tmp.path()).unwrap()
}

#[test]
fn batched_output_is_byte_identical_to_per_record() {
    let (records, header) = read_all(fixture());
    for workers in [1usize, 2, 4] {
        let per_record = write_per_record(&records, &header, workers);
        for batch_size in [1usize, 3, records.len()] {
            let batched = write_batched(&records, &header, workers, batch_size);
            assert_eq!(
                batched, per_record,
                "workers={workers} batch_size={batch_size}: batched output must be byte-identical"
            );
        }
    }
}

/// A multi-block file (records repeated until the encoded stream spans several
/// 64 KiB BGZF blocks) is the realistic case: the batch writer's framing thread
/// must reproduce the same block boundaries as per-record writing, since both
/// feed the same `MultithreadedWriter`.
#[test]
fn batched_output_matches_over_multi_block_file() {
    let (originals, header) = read_all(fixture());
    let mut records = Vec::new();
    for _ in 0..2000 {
        records.extend(originals.iter().cloned());
    }

    let per_record = write_per_record(&records, &header, 2);
    let batched = write_batched(&records, &header, 2, 4096);
    assert_eq!(
        batched, per_record,
        "multi-block batched output must be byte-identical to per-record"
    );

    // And the records read back are exactly the originals, in order.
    let tmp = tempfile::NamedTempFile::new().unwrap();
    fs::write(tmp.path(), &batched).unwrap();
    let (round, _) = read_all(tmp.path());
    assert_eq!(round.len(), records.len());
    for (i, r) in round.iter().enumerate() {
        assert_eq!(r.as_bytes(), records[i].as_bytes(), "record {i}");
    }
}

/// An empty batch is a no-op and finishing with no batches still produces a
/// valid header-only BAM with the EOF block.
#[test]
fn empty_batches_and_no_batches_are_valid() {
    let (_, header) = read_all(fixture());

    let header_only = {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut writer =
            rsomics_bamio::create_with_workers(tmp.path(), NonZero::new(1).unwrap()).unwrap();
        writer.write_alignment_header(&header).unwrap();
        let mut batch_writer = BatchBamWriter::new(writer);
        batch_writer.write_records_batch(Vec::new()).unwrap();
        batch_writer.finish().unwrap();
        fs::read(tmp.path()).unwrap()
    };

    let (round, _) = read_all_from_bytes(&header_only);
    assert!(round.is_empty(), "header-only BAM has no records");
}

fn read_all_from_bytes(bytes: &[u8]) -> (Vec<RawRecord>, noodles::sam::Header) {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    fs::write(tmp.path(), bytes).unwrap();
    read_all(tmp.path())
}
