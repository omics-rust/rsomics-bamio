# rsomics-bamio

`rsomics-bamio` provides validated alignment records, indexed alignment input,
and bounded parallel BGZF I/O for alignment-consuming rsomics products.
Current consumers are `rsomics-bam` and `rsomics-call`; the same policy-free
record contract is planned for `rsomics-methyl` and `rsomics-minimap2`.

`RawRecord` owns a validated BAM record payload for in-place edits.
`RecordReader` yields borrowing `RecordRef` values directly from an inflated
BGZF block when possible. Both expose the same fixed-field, CIGAR, sequence,
quality, and auxiliary-data accessors. `RawRecordEncoder` converts any noodles
alignment record into the same validated payload without exposing file-format
policy.

```rust
use rsomics_bamio::{open_with_workers, raw::RecordReader};

let mut reader = open_with_workers(path, workers)?;
reader.read_header()?;

let mut records = RecordReader::new(reader.get_mut());
while let Some(record) = records.next()? {
    let flags = record.flags();
}
```

Malformed record layouts return structured errors before their fields are
exposed. Product-specific filtering, command defaults, and output policy stay
in the consuming product.

`open_indexed_alignment` opens BGZF SAM, BAM, or CRAM input with the usual
appended index names and the common alternative names such as `sample.bai`.
CRAM callers can attach an indexed reference without duplicating repository
setup in each product.

License: MIT OR Apache-2.0.
