//! Lightweight shared IO helpers for Dino Tools.

use std::path::Path;

/// Common biological sequence file families.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SequenceFormat {
    Fasta,
    Fastq,
    Unknown,
}

impl SequenceFormat {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Fasta => "fasta",
            Self::Fastq => "fastq",
            Self::Unknown => "unknown",
        }
    }
}

/// Guess a sequence format from a filename extension.
///
/// This intentionally stays conservative. Content sniffing should be added only
/// when a promoted tool needs it.
#[must_use]
pub fn guess_sequence_format(path: impl AsRef<Path>) -> SequenceFormat {
    let Some(file_name) = path.as_ref().file_name().and_then(|name| name.to_str()) else {
        return SequenceFormat::Unknown;
    };

    let lower = file_name.to_ascii_lowercase();
    let trimmed = lower
        .strip_suffix(".gz")
        .or_else(|| lower.strip_suffix(".zst"))
        .or_else(|| lower.strip_suffix(".bz2"))
        .unwrap_or(&lower);

    if matches_fasta(trimmed) {
        SequenceFormat::Fasta
    } else if matches_fastq(trimmed) {
        SequenceFormat::Fastq
    } else {
        SequenceFormat::Unknown
    }
}

fn matches_fasta(path: &str) -> bool {
    path.ends_with(".fa") || path.ends_with(".fasta") || path.ends_with(".fna")
}

fn matches_fastq(path: &str) -> bool {
    path.ends_with(".fq") || path.ends_with(".fastq")
}

#[cfg(test)]
mod tests {
    use super::{guess_sequence_format, SequenceFormat};

    #[test]
    fn guesses_plain_sequence_extensions() {
        assert_eq!(guess_sequence_format("reads.fq"), SequenceFormat::Fastq);
        assert_eq!(guess_sequence_format("ref.fasta"), SequenceFormat::Fasta);
    }

    #[test]
    fn guesses_compressed_sequence_extensions() {
        assert_eq!(
            guess_sequence_format("reads.fastq.gz"),
            SequenceFormat::Fastq
        );
        assert_eq!(guess_sequence_format("ref.fa.zst"), SequenceFormat::Fasta);
    }

    #[test]
    fn unknown_for_unrecognized_extension() {
        assert_eq!(guess_sequence_format("notes.txt"), SequenceFormat::Unknown);
    }
}
