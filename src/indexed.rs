use std::fs::File;
use std::io::{BufReader, Read};
use std::ops::{Deref, DerefMut};
use std::path::{Path, PathBuf};

use noodles::{bam, bgzf, cram, csi, fasta, sam};
use noodles_util::alignment;
use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::Read as _;

use crate::raw::{self, RawRecord, RawRecordEncoder};

/// A format-independent indexed SAM, BAM, or CRAM reader.
pub struct IndexedAlignmentReader {
    inner: alignment::io::IndexedReader<File>,
    input: PathBuf,
    reference: Option<PathBuf>,
}

impl Deref for IndexedAlignmentReader {
    type Target = alignment::io::IndexedReader<File>;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl DerefMut for IndexedAlignmentReader {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

/// Opens an indexed alignment file and optionally attaches an indexed reference.
///
/// Both appended index names such as `sample.bam.bai` and common alternative
/// names such as `sample.bai` are accepted. Standard input and uncompressed SAM
/// or BAM input are rejected because they cannot support indexed queries.
pub fn open_indexed_alignment(
    input: &Path,
    reference: Option<&Path>,
) -> Result<IndexedAlignmentReader> {
    if input == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "region queries require a file-backed alignment input".to_owned(),
        ));
    }
    let (format, compression) = detect_source(input)?;
    if matches!(
        (format, compression),
        (Format::Sam | Format::Bam, Compression::None)
    ) {
        return Err(RsomicsError::InvalidInput(format!(
            "region queries require BGZF SAM, BAM, or CRAM input: {}",
            input.display()
        )));
    }

    let mut builder = alignment::io::indexed_reader::Builder::default();
    if let Some(reference) = reference {
        builder = builder.set_reference_sequence_repository(reference_repository(reference)?);
    }
    builder = set_alternative_index(builder, input, format)?;
    let inner = builder.build_from_path(input).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "opening indexed alignment input {}: {error}",
            input.display()
        ))
    })?;
    Ok(IndexedAlignmentReader {
        inner,
        input: input.to_path_buf(),
        reference: reference.map(Path::to_path_buf),
    })
}

/// Visits every alignment record as a validated BAM payload.
///
/// BAM input is read without decoding and re-encoding each record. SAM and CRAM
/// records are normalized into the same representation before being visited.
pub fn visit_raw_alignment_records(
    reader: &mut IndexedAlignmentReader,
    header: &sam::Header,
    mut visit: impl FnMut(RawRecord) -> Result<()>,
) -> Result<()> {
    if matches!(reader.inner, alignment::io::IndexedReader::Cram(_)) {
        visit_htslib_records(&reader.input, reader.reference.as_deref(), visit)?;
        return Ok(());
    }

    let mut scan = open_indexed_alignment(&reader.input, reader.reference.as_deref())?;
    scan.read_header().map_err(RsomicsError::Io)?;
    if let alignment::io::IndexedReader::Bam(reader) = &mut scan.inner {
        loop {
            let mut record = RawRecord::default();
            if raw::read_record(reader.get_mut(), &mut record)? == 0 {
                break;
            }
            visit(record)?;
        }
        return Ok(());
    }

    let mut encoder = RawRecordEncoder::new();
    for result in scan.records(header) {
        let record = result.map_err(RsomicsError::Io)?;
        let record = encoder.encode(header, record.as_ref())?;
        visit(record)?;
    }
    Ok(())
}

fn visit_htslib_records(
    input: &Path,
    reference: Option<&Path>,
    mut visit: impl FnMut(RawRecord) -> Result<()>,
) -> Result<()> {
    let mut reader = rust_htslib::bam::Reader::from_path(input).map_err(|error| {
        RsomicsError::InvalidInput(format!("opening CRAM input {}: {error}", input.display()))
    })?;
    if let Some(reference) = reference {
        reader.set_reference(reference).map_err(|error| {
            RsomicsError::ConfigError(format!(
                "attaching CRAM reference {}: {error}",
                reference.display()
            ))
        })?;
    }

    let mut record = rust_htslib::bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|error| {
            RsomicsError::InvalidInput(format!("reading CRAM input {}: {error}", input.display()))
        })?;
        visit(raw::from_htslib_record(&record)?)?;
    }
    Ok(())
}

#[derive(Clone, Copy)]
enum Format {
    Sam,
    Bam,
    Cram,
}

#[derive(Clone, Copy)]
enum Compression {
    None,
    Bgzf,
}

fn set_alternative_index(
    mut builder: alignment::io::indexed_reader::Builder,
    input: &Path,
    format: Format,
) -> Result<alignment::io::indexed_reader::Builder> {
    match format {
        Format::Sam => {
            let appended = append_extension(input, "csi");
            let alternative = input.with_extension("csi");
            if !index_exists(input, &appended)? && index_exists(input, &alternative)? {
                let index = csi::fs::read(&alternative)
                    .map_err(|error| index_error(input, &alternative, error))?;
                builder = builder.set_index(index);
            }
        }
        Format::Bam => {
            let appended_bai = append_extension(input, "bai");
            let appended_csi = append_extension(input, "csi");
            if !index_exists(input, &appended_bai)? && !index_exists(input, &appended_csi)? {
                let alternative_bai = input.with_extension("bai");
                let alternative_csi = input.with_extension("csi");
                if index_exists(input, &alternative_bai)? {
                    let index = bam::bai::fs::read(&alternative_bai)
                        .map_err(|error| index_error(input, &alternative_bai, error))?;
                    builder = builder.set_index(index);
                } else if index_exists(input, &alternative_csi)? {
                    let index = csi::fs::read(&alternative_csi)
                        .map_err(|error| index_error(input, &alternative_csi, error))?;
                    builder = builder.set_index(index);
                }
            }
        }
        Format::Cram => {
            let appended = append_extension(input, "crai");
            let alternative = input.with_extension("crai");
            if !index_exists(input, &appended)? && index_exists(input, &alternative)? {
                let index = cram::crai::fs::read(&alternative)
                    .map_err(|error| index_error(input, &alternative, error))?;
                builder = builder.set_index(index);
            }
        }
    }
    Ok(builder)
}

fn append_extension(path: &Path, extension: &str) -> PathBuf {
    let mut path = path.as_os_str().to_owned();
    path.push(".");
    path.push(extension);
    PathBuf::from(path)
}

fn index_exists(input: &Path, index: &Path) -> Result<bool> {
    index
        .try_exists()
        .map_err(|error| index_error(input, index, error))
}

fn detect_source(input: &Path) -> Result<(Format, Compression)> {
    let mut source = BufReader::new(File::open(input).map_err(|error| {
        RsomicsError::InvalidInput(format!("opening {}: {error}", input.display()))
    })?);
    let mut magic = [0; 4];
    source.read_exact(&mut magic).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "detecting alignment format for {}: {error}",
            input.display()
        ))
    })?;
    if magic == *b"CRAM" {
        return Ok((Format::Cram, Compression::None));
    }
    if magic == *b"BAM\x01" {
        return Ok((Format::Bam, Compression::None));
    }
    if magic[..2] != [0x1f, 0x8b] {
        return Ok((Format::Sam, Compression::None));
    }

    let file = File::open(input).map_err(|error| {
        RsomicsError::InvalidInput(format!("opening {}: {error}", input.display()))
    })?;
    let mut reader = bgzf::io::Reader::new(file);
    reader.read_exact(&mut magic).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "detecting BGZF alignment format for {}: {error}",
            input.display()
        ))
    })?;
    Ok((
        if magic == *b"BAM\x01" {
            Format::Bam
        } else {
            Format::Sam
        },
        Compression::Bgzf,
    ))
}

fn reference_repository(path: &Path) -> Result<fasta::Repository> {
    fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map(fasta::repository::adapters::IndexedReader::new)
        .map(fasta::Repository::new)
        .map_err(|error| {
            RsomicsError::ConfigError(format!(
                "opening indexed reference {}: {error}",
                path.display()
            ))
        })
}

fn index_error(input: &Path, index: &Path, error: std::io::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "reading alignment index {} for {}: {error}",
        index.display(),
        input.display()
    ))
}

#[cfg(test)]
mod tests {
    use noodles::core::Region;
    use noodles::sam;
    use noodles::sam::alignment::io::Write as _;

    use super::*;

    #[test]
    fn opens_alternative_bai_and_queries_records() {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("records.bam");
        let mut source = sam::io::Reader::new(
            b"@HD\tVN:1.6\tSO:coordinate\n\
              @SQ\tSN:chr1\tLN:20\n\
              read\t0\tchr1\t3\t60\t1M\t*\t0\t0\tA\tI\n"
                .as_slice(),
        );
        let header = source.read_header().unwrap();
        let records = source
            .record_bufs(&header)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        let mut writer = bam::io::Writer::new(File::create(&input).unwrap());
        writer.write_header(&header).unwrap();
        for record in records {
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish().unwrap();
        let index = bam::fs::index(&input).unwrap();
        bam::bai::fs::write(input.with_extension("bai"), &index).unwrap();

        let mut reader = open_indexed_alignment(&input, None).unwrap();
        let header = reader.read_header().unwrap();
        let region: Region = "chr1:1-5".parse().unwrap();
        let records = reader
            .query(&header, &region)
            .unwrap()
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(records.len(), 1);

        let mut reader = open_indexed_alignment(&input, None).unwrap();
        let header = reader.read_header().unwrap();
        let mut records = Vec::new();
        visit_raw_alignment_records(&mut reader, &header, |record| {
            records.push(record);
            Ok(())
        })
        .unwrap();
        assert_eq!(records.len(), 1);
        let record = &records[0];
        assert_eq!(record.name(), b"read");
        assert_eq!(record.reference_sequence_id(), 0);
        assert_eq!(record.alignment_start(), 2);
        assert_eq!(record.mapping_quality(), 60);
        assert_eq!(record.cigar_ops().collect::<Vec<_>>(), [(0, 1)]);
        assert_eq!(record.sequence_len(), 1);
        assert_eq!(record.quality_scores(), [40]);
    }

    #[test]
    fn visits_cram_records_from_slices() {
        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("reference.fa");
        std::fs::write(&reference, b">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
        std::fs::write(
            append_extension(&reference, "fai"),
            b"chr1\t20\t6\t20\t21\n",
        )
        .unwrap();

        let input = directory.path().join("records.cram");
        let mut source = sam::io::Reader::new(
            b"@HD\tVN:1.6\tSO:coordinate\n\
              @SQ\tSN:chr1\tLN:20\n\
              read\t0\tchr1\t3\t60\t4M\t*\t0\t0\tGTAC\tIIII\n"
                .as_slice(),
        );
        let header = source.read_header().unwrap();
        let records = source
            .record_bufs(&header)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        let repository = reference_repository(&reference).unwrap();
        let mut writer = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(&input)
            .unwrap();
        writer.write_header(&header).unwrap();
        for record in records {
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish(&header).unwrap();

        let index = cram::fs::index(&input).unwrap();
        cram::crai::fs::write(append_extension(&input, "crai"), &index).unwrap();

        let mut reader = open_indexed_alignment(&input, Some(&reference)).unwrap();
        let header = reader.read_header().unwrap();
        let mut records = Vec::new();
        visit_raw_alignment_records(&mut reader, &header, |record| {
            records.push(record);
            Ok(())
        })
        .unwrap();
        assert_eq!(records.len(), 1);
        let record = &records[0];
        assert_eq!(record.name(), b"read");
        assert_eq!(record.reference_sequence_id(), 0);
        assert_eq!(record.alignment_start(), 2);
        assert_eq!(record.mapping_quality(), 60);
        assert_eq!(record.cigar_ops().collect::<Vec<_>>(), [(0, 4)]);
        assert_eq!(record.sequence_len(), 4);
        assert_eq!(record.quality_scores(), [40; 4]);
    }

    #[test]
    fn rejects_uncompressed_alignment_input() {
        let input = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            input.path(),
            b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20\nread\t0\tchr1\t3\t60\t1M\t*\t0\t0\tA\tI\n",
        )
        .unwrap();
        assert!(matches!(
            open_indexed_alignment(input.path(), None),
            Err(RsomicsError::InvalidInput(message))
                if message.contains("require BGZF SAM, BAM, or CRAM")
        ));
    }
}
