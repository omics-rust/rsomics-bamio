# rsomics-bamio

`rsomics-bamio` provides validated raw BAM records and bounded parallel BGZF
I/O for alignment-consuming rsomics products. Initial consumers are
`rsomics-bam`; the same policy-free record contract is planned for
`rsomics-methyl`, `rsomics-call`, and `rsomics-minimap2`.

`RawRecord` owns a validated BAM record payload for in-place edits.
`RecordReader` yields borrowing `RecordRef` values directly from an inflated
BGZF block when possible. Both expose the same fixed-field, CIGAR, sequence,
quality, and auxiliary-data accessors.

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

License: MIT OR Apache-2.0.
