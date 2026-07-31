use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use noodles::{bam, bgzf, cram, csi, fasta};
use noodles_util::alignment;
use rsomics_common::{Result, RsomicsError};

/// A format-independent indexed SAM, BAM, or CRAM reader.
pub type IndexedAlignmentReader = alignment::io::IndexedReader<File>;

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
    builder.build_from_path(input).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "opening indexed alignment input {}: {error}",
            input.display()
        ))
    })
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
