# rsomics-bamio

Layer-A primitive: a parallel-BGZF BAM reader shared by the `rsomics-bam-*`
tool family.

`samtools` inflates BGZF on a single thread by default. A single-threaded
pure-Rust reader (zlib-rs) loses to htslib's libdeflate inner loop, so reading
BAM through a worker pool is what puts the rsomics BAM tools ahead of
`samtools` default invocations on multi-core hosts.

```rust
use rsomics_bamio::open_parallel;

let mut reader = open_parallel(path)?;
let header = reader.read_header()?;
for result in reader.records() {
    let record = result?;
    // ...
}
```

## Raw-record paths

For tools that work on the on-disk record bytes directly (no decode of
seq/qual/cigar), `rsomics_bamio::raw` exposes two views over one decoding core:

- `RawRecord` — **owns** its payload bytes. Edit fixed-offset fields and the aux
  tail in place (`set_flag_bits`, `set_aux`, …), then `write_record`. Read via
  `read_record(reader.get_mut(), &mut rec)`.
- `RecordReader` / `RecordRef` — a **borrowing, allocation-free** scan for
  read-only passes. `RecordReader::new(reader.get_mut())` hands out each record
  as a `RecordRef` borrowed straight out of the BGZF reader's inflated block
  buffer (only a block-straddling record is copied into a reused scratch). Same
  field accessors as `RawRecord`; one BGZF block of resident memory regardless of
  record count.

```rust
use rsomics_bamio::{open_with_workers, raw::RecordReader};

let mut reader = open_with_workers(path, workers)?;
reader.read_header()?;
let mut scan = RecordReader::new(reader.get_mut());
while let Some(rec) = scan.next()? {
    let _ = (rec.flags(), rec.alignment_start(), rec.quality_scores());
}
```

License: MIT OR Apache-2.0.
