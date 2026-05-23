use std::num::NonZero;
use std::path::Path;

use rsomics_bamio::raw::{self, FLAG_DUPLICATE, RawRecord, RecordReader};

fn fixture() -> &'static Path {
    Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/golden/small.bam"
    ))
}

/// Fields that the raw-edit path must pass through byte-for-byte. We snapshot
/// them before editing and assert they survive the read→edit→write→read cycle
/// unchanged, while the flag we touched does change.
#[derive(Clone, PartialEq, Eq, Debug)]
struct Passthrough {
    name: Vec<u8>,
    pos: i32,
    ref_id: i32,
    cigar: Vec<(u8, u32)>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

fn snapshot(r: &RawRecord) -> Passthrough {
    Passthrough {
        name: r.name().to_vec(),
        pos: r.alignment_start(),
        ref_id: r.reference_sequence_id(),
        cigar: r.cigar_ops().collect(),
        seq: seq_bytes(r),
        qual: r.quality_scores().to_vec(),
    }
}

/// Re-extract the packed seq nibbles directly from the payload so an accidental
/// shift of any field would be caught (the public API never exposes seq bytes,
/// since seq is pure pass-through — this reaches into `as_bytes()` for the test).
fn seq_bytes(r: &RawRecord) -> Vec<u8> {
    let bytes = r.as_bytes();
    let name_len = usize::from(bytes[8]);
    let n_cigar = u16::from_le_bytes(bytes[12..14].try_into().unwrap()) as usize;
    let l_seq = u32::from_le_bytes(bytes[16..20].try_into().unwrap()) as usize;
    let start = 32 + name_len + n_cigar * 4;
    bytes[start..start + l_seq.div_ceil(2)].to_vec()
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

#[test]
fn raw_edit_sets_flag_and_passes_seq_qual_through() {
    let (originals, header) = read_all(fixture());
    assert_eq!(originals.len(), 10);

    let snapshots: Vec<Passthrough> = originals.iter().map(snapshot).collect();

    // Edit a deterministic subset: every 3rd record gets 0x400 set.
    let edited: Vec<bool> = (0..originals.len()).map(|i| i % 3 == 0).collect();

    let tmp = tempfile::NamedTempFile::new().unwrap();
    {
        let mut writer =
            rsomics_bamio::create_with_workers(tmp.path(), NonZero::new(1).unwrap()).unwrap();
        use noodles::sam::alignment::io::Write as _;
        writer.write_alignment_header(&header).unwrap();
        for (i, r) in originals.iter().enumerate() {
            let mut r = r.clone();
            if edited[i] {
                r.set_flag_bits(FLAG_DUPLICATE);
            }
            raw::write_record(writer.get_mut(), &r).unwrap();
        }
        // writer drop appends the BGZF EOF block.
    }

    let (round, _) = read_all(tmp.path());
    assert_eq!(round.len(), originals.len());

    for (i, r) in round.iter().enumerate() {
        assert_eq!(
            snapshot(r),
            snapshots[i],
            "record {i}: seq/qual/name/cigar/pos must be byte-identical"
        );
        let orig_flag = originals[i].flags();
        let expected = if edited[i] {
            orig_flag | FLAG_DUPLICATE
        } else {
            orig_flag
        };
        assert_eq!(r.flags(), expected, "record {i}: flag mismatch");
        if edited[i] && orig_flag & FLAG_DUPLICATE == 0 {
            assert_ne!(r.flags(), orig_flag, "record {i}: flag should have changed");
            assert_eq!(r.flags() & FLAG_DUPLICATE, FLAG_DUPLICATE);
        }
    }
}

#[test]
fn clear_flag_bits_only_touches_target_bits() {
    let (originals, _) = read_all(fixture());
    // Find the record that already carries 0x400 in the golden fixture.
    let idx = originals
        .iter()
        .position(|r| r.flags() & FLAG_DUPLICATE != 0)
        .expect("golden fixture has a 0x400 record");
    let mut r = originals[idx].clone();
    let before = r.flags();
    r.clear_flag_bits(FLAG_DUPLICATE);
    assert_eq!(r.flags(), before & !FLAG_DUPLICATE);
    assert_eq!(r.flags() & FLAG_DUPLICATE, 0);
}

/// The borrowing `RecordReader` must decode every field byte-identically to the
/// owning `read_record` path over the same stream, at both worker counts (single
/// = plain reader, multi = worker pool). This is the correctness contract the
/// 13+ dependent tools rely on when they switch to the zero-copy scan.
#[test]
fn borrowed_reader_matches_owned_read() {
    for workers in [1usize, 4] {
        let nz = NonZero::new(workers).unwrap();

        let owned: Vec<RawRecord> = {
            let mut reader = rsomics_bamio::open_with_workers(fixture(), nz).unwrap();
            reader.read_header().unwrap();
            let mut out = Vec::new();
            let mut rec = RawRecord::default();
            while raw::read_record(reader.get_mut(), &mut rec).unwrap() != 0 {
                out.push(rec.clone());
            }
            out
        };

        let mut reader = rsomics_bamio::open_with_workers(fixture(), nz).unwrap();
        reader.read_header().unwrap();
        let mut scanner = RecordReader::new(reader.get_mut());
        let mut i = 0;
        while let Some(r) = scanner.next().unwrap() {
            let o = &owned[i];
            assert_eq!(r.as_bytes(), o.as_bytes(), "workers={workers} record {i}");
            assert_eq!(r.flags(), o.flags());
            assert_eq!(r.reference_sequence_id(), o.reference_sequence_id());
            assert_eq!(r.alignment_start(), o.alignment_start());
            assert_eq!(r.mapping_quality(), o.mapping_quality());
            assert_eq!(r.name(), o.name());
            assert_eq!(
                r.cigar_ops().collect::<Vec<_>>(),
                o.cigar_ops().collect::<Vec<_>>()
            );
            assert_eq!(r.sequence_len(), o.sequence_len());
            assert_eq!(r.quality_scores(), o.quality_scores());
            i += 1;
        }
        assert_eq!(i, owned.len(), "workers={workers}: record count");
    }
}

/// Exercise the block-straddling spill branch: a file large enough to span many
/// BGZF blocks puts some records across a block boundary, where `RecordReader`
/// must copy into its scratch buffer instead of borrowing. The borrowed scan must
/// still reproduce every owned record byte-for-byte.
#[test]
fn borrowed_reader_handles_multi_block_files() {
    let (originals, header) = read_all(fixture());

    // Repeat the golden records until the encoded file spans several 64 KiB BGZF
    // blocks, forcing at least one record to straddle a boundary.
    let tmp = tempfile::NamedTempFile::new().unwrap();
    {
        let mut writer =
            rsomics_bamio::create_with_workers(tmp.path(), NonZero::new(1).unwrap()).unwrap();
        use noodles::sam::alignment::io::Write as _;
        writer.write_alignment_header(&header).unwrap();
        for _ in 0..5000 {
            for r in &originals {
                raw::write_record(writer.get_mut(), r).unwrap();
            }
        }
    }

    let owned = {
        let (recs, _) = read_all(tmp.path());
        recs
    };
    assert_eq!(owned.len(), originals.len() * 5000);

    let mut reader =
        rsomics_bamio::open_with_workers(tmp.path(), NonZero::new(1).unwrap()).unwrap();
    reader.read_header().unwrap();
    let mut scanner = RecordReader::new(reader.get_mut());
    let mut i = 0;
    while let Some(r) = scanner.next().unwrap() {
        assert_eq!(r.as_bytes(), owned[i].as_bytes(), "record {i}");
        i += 1;
    }
    assert_eq!(i, owned.len());
}

#[test]
fn aux_append_remove_replace_roundtrips() {
    let (originals, _) = read_all(fixture());
    let mut r = originals[0].clone();
    let before_len = r.as_bytes().len();

    // Append a new MQ tag (type 'C', one byte).
    assert!(r.aux_value(*b"MQ").is_none());
    r.append_aux(*b"MQ", b'C', &[60]);
    assert_eq!(r.aux_type(*b"MQ"), Some(b'C'));
    assert_eq!(r.aux_value(*b"MQ"), Some(&[60u8][..]));

    // Replace it with a different value.
    r.set_aux(*b"MQ", b'C', &[42]);
    assert_eq!(r.aux_value(*b"MQ"), Some(&[42u8][..]));

    // Remove it; payload returns to its original length.
    assert!(r.remove_aux(*b"MQ"));
    assert!(r.aux_value(*b"MQ").is_none());
    assert_eq!(r.as_bytes().len(), before_len);
    assert!(!r.remove_aux(*b"MQ"));
}
