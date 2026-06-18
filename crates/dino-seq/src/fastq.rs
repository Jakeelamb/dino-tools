use std::io::Read;

use crate::error::{FastqError, Result};
use crate::fastq_frame;

mod chunk;
mod frame;
mod pair;
mod record;

pub use chunk::{FastqChunkConfig, FastqChunkSinkExt, FastqChunkStats, FastqRecordSink};
use frame::{ChunkVisitContext, frame_records, visit_chunk_in_slab};
use frame::{count_records_in_slab, visit_records_in_slab};
pub use pair::{
    FastqPair, InterleavedPairs, PairValidation, PairedFastqBatch, PairedFastqPairs,
    PairedFastqReader, PairedRecords, PairingMode, paired_records,
};
use pair::{validate_even_pair_count, validate_interleaved_pair_ids, validate_paired_batches};
pub use record::{FastqRecord, FastqVisitRecord, RecordRef, strip_pair_suffix};

const DEFAULT_SLAB_SIZE: usize = 1024 * 1024;

/// Configuration for FASTQ batch readers.
///
/// The default is a validated, unpaired reader with a 1 MiB slab and full pair
/// identifier validation when pairing APIs are used.
#[derive(Debug, Clone)]
pub struct FastqConfig {
    /// Target slab size in bytes.
    ///
    /// Values below 1024 are raised to 1024. Larger slabs reduce carry
    /// frequency for long records but increase the reusable buffer size.
    pub slab_size: usize,
    /// Validate FASTQ structure while framing records.
    ///
    /// Validation checks the leading `@`, leading `+`, and sequence/quality
    /// length equality. Disable only for trusted inputs.
    pub validate: bool,
    /// Whether a single stream should be interpreted as unpaired or interleaved.
    pub pairing: PairingMode,
    /// Identifier validation policy for paired APIs.
    pub pair_validation: PairValidation,
}

impl Default for FastqConfig {
    fn default() -> Self {
        Self {
            slab_size: DEFAULT_SLAB_SIZE,
            validate: true,
            pairing: PairingMode::None,
            pair_validation: PairValidation::Full,
        }
    }
}

impl FastqConfig {
    /// Treat a single FASTQ stream as adjacent interleaved read pairs.
    pub fn interleaved(mut self) -> Self {
        self.pairing = PairingMode::Interleaved;
        self
    }

    /// Set the paired-read identifier validation policy.
    pub fn pair_validation(mut self, pair_validation: PairValidation) -> Self {
        self.pair_validation = pair_validation;
        self
    }
}

/// Aggregate statistics for sequence-only FASTQ workloads.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FastqStats {
    /// Number of records observed.
    pub records: u64,
    /// Number of sequence bases observed.
    pub bases: u64,
    /// Number of quality bytes observed.
    pub qualities: u64,
    /// Number of header bytes observed, including leading `@`.
    pub name_bytes: u64,
    /// Lightweight deterministic checksum over record shape and first bases.
    pub checksum: u64,
}

/// A batch of borrowed FASTQ records from one reader slab.
///
/// The batch is invalidated by the next mutable call on the reader that
/// produced it. Process records before requesting another batch.
#[derive(Debug)]
pub struct FastqBatch<'a> {
    bytes: &'a [u8],
    records: &'a [RecordRef],
    base_offset: u64,
    first_record_index: u64,
    pair_validation: PairValidation,
}

impl<'a> FastqBatch<'a> {
    /// Number of records in the batch.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the batch has no records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Raw slab bytes backing this batch.
    pub fn bytes(&self) -> &'a [u8] {
        self.bytes
    }

    /// Record ranges within [`bytes`](Self::bytes).
    pub fn record_refs(&self) -> &'a [RecordRef] {
        self.records
    }

    /// Absolute byte offset of `bytes()[0]` in the original stream.
    pub fn base_offset(&self) -> u64 {
        self.base_offset
    }

    /// Zero-based index of the first record in this batch.
    pub fn first_record_index(&self) -> u64 {
        self.first_record_index
    }

    /// Pair validation mode inherited from the producing reader.
    pub fn pair_validation(&self) -> PairValidation {
        self.pair_validation
    }

    /// Iterate borrowed record views.
    pub fn records(&self) -> impl Iterator<Item = FastqRecord<'a>> + 'a {
        let bytes = self.bytes;
        let records: &'a [RecordRef] = self.records;
        records
            .iter()
            .map(move |record| FastqRecord { bytes, record })
    }

    /// Validate and iterate adjacent interleaved read pairs.
    ///
    /// Returns an error if the batch has an odd record count or mate
    /// identifiers fail the configured [`PairValidation`] mode.
    pub fn interleaved_pairs(&'a self) -> Result<InterleavedPairs<'a>> {
        validate_even_pair_count(self)?;
        validate_interleaved_pair_ids(self, self.pair_validation)?;
        Ok(InterleavedPairs {
            batch: self,
            next: 0,
        })
    }

    /// Validate and zip this batch with a mate batch from a separate reader.
    ///
    /// Identifier checks use this batch's configured [`PairValidation`] mode.
    pub fn paired_with(&'a self, mate: &'a FastqBatch<'a>) -> Result<PairedRecords<'a>> {
        validate_paired_batches(self, mate, self.pair_validation)?;
        Ok(PairedRecords {
            first: self,
            second: mate,
            next: 0,
        })
    }

    fn record_at(&self, index: usize) -> FastqRecord<'a> {
        let records: &'a [RecordRef] = self.records;
        FastqRecord {
            bytes: self.bytes,
            record: &records[index],
        }
    }
}

/// Slab-based FASTQ reader over any [`Read`] input.
///
/// The reader frames four-line FASTQ records from a reusable byte slab. Records
/// crossing slab boundaries are carried into the next read. Returned batches
/// borrow from the reader and must be consumed before calling [`next_batch`]
/// again.
///
/// [`next_batch`]: Self::next_batch
#[derive(Debug)]
pub struct FastqReader<R> {
    reader: R,
    config: FastqConfig,
    buf: Vec<u8>,
    len: usize,
    next_start: usize,
    chunk_resume_start: usize,
    base_offset: u64,
    record_index: u64,
    eof: bool,
    records: Vec<RecordRef>,
}

impl<R: Read> FastqReader<R> {
    /// Create a reader with [`FastqConfig::default`].
    pub fn new(reader: R) -> Self {
        Self::with_config(reader, FastqConfig::default())
    }

    /// Create a reader with explicit configuration.
    pub fn with_config(reader: R, config: FastqConfig) -> Self {
        let slab_size = config.slab_size.max(1024);
        let buf = vec![0; slab_size];
        Self {
            reader,
            config: FastqConfig {
                slab_size,
                ..config
            },
            buf,
            len: 0,
            next_start: 0,
            chunk_resume_start: 0,
            base_offset: 0,
            record_index: 0,
            eof: false,
            records: Vec::new(),
        }
    }

    /// Return the wrapped reader.
    pub fn into_inner(self) -> R {
        self.reader
    }

    /// Read the next batch of FASTQ records.
    ///
    /// Returns `Ok(None)` at EOF. Format errors include byte offset, record
    /// index, and line index when the failing location is known.
    pub fn next_batch(&mut self) -> Result<Option<FastqBatch<'_>>> {
        self.compact_carry();
        self.fill_slab()?;

        self.records.clear();
        let first_record_index = self.record_index;
        self.next_start = frame_records(
            &self.buf[..self.len],
            self.eof,
            self.config.validate,
            self.base_offset,
            first_record_index,
            &mut self.records,
        )?;
        self.align_interleaved_batch(first_record_index)?;

        self.record_index += self.records.len() as u64;

        if self.records.is_empty() {
            if self.len == 0 && self.eof {
                return Ok(None);
            }
            return Err(FastqError::RecordTooLarge {
                slab_size: self.config.slab_size,
            });
        }

        Ok(Some(FastqBatch {
            bytes: &self.buf[..self.len],
            records: &self.records,
            base_offset: self.base_offset,
            first_record_index,
            pair_validation: self.config.pair_validation,
        }))
    }

    /// Visit every record in the stream without building a batch side table.
    ///
    /// This is the fastest public parse-only path for single-pass consumers.
    /// It still honors [`FastqConfig::validate`] for record structure, but it
    /// yields individual records and does not perform paired-end identifier
    /// validation. Use [`next_batch`](Self::next_batch),
    /// [`interleaved_pairs`](FastqBatch::interleaved_pairs), or
    /// [`PairedFastqReader`] when pair validation is required.
    pub fn visit_records<F>(&mut self, mut visit: F) -> Result<()>
    where
        F: FnMut(FastqVisitRecord<'_>) -> Result<()>,
    {
        loop {
            self.compact_carry();
            self.fill_slab()?;

            let (next_start, records) = visit_records_in_slab(
                &self.buf[..self.len],
                self.eof,
                self.config.validate,
                self.base_offset,
                self.record_index,
                &mut visit,
            )?;
            self.next_start = next_start;
            self.record_index += records;

            if records == 0 {
                if self.len == 0 && self.eof {
                    return Ok(());
                }
                return Err(FastqError::RecordTooLarge {
                    slab_size: self.config.slab_size,
                });
            }
        }
    }

    /// Count remaining records and bases without building record views.
    ///
    /// This consumes the reader to EOF. It validates the same four-line FASTQ
    /// structure as [`visit_records`](Self::visit_records) when
    /// [`FastqConfig::validate`] is enabled, but it avoids the per-record
    /// visitor callback and borrowed-record construction.
    #[inline]
    pub fn count_records(&mut self) -> Result<FastqStats> {
        let mut stats = FastqStats::default();
        loop {
            self.compact_carry();
            self.fill_slab()?;

            let (next_start, records) = count_records_in_slab(
                &self.buf[..self.len],
                self.eof,
                self.config.validate,
                self.base_offset,
                self.record_index,
                &mut stats,
            )?;
            self.next_start = next_start;
            self.record_index += records;

            if records == 0 {
                if self.len == 0 && self.eof {
                    return Ok(stats);
                }
                return Err(FastqError::RecordTooLarge {
                    slab_size: self.config.slab_size,
                });
            }
        }
    }

    /// Emit one resumable chunk of FASTQ records into a caller-owned sink.
    ///
    /// This is the low-overhead path for downstream tools that already have a
    /// target output representation. It does not build the [`RecordRef`] side
    /// table used by [`next_batch`](Self::next_batch), and it returns after the
    /// configured chunk limit so callers can interleave parsing with their own
    /// pipeline stages.
    ///
    /// Returns `Ok(None)` only when EOF is reached before emitting another
    /// record. A subsequent call resumes at the first unconsumed record.
    pub fn next_chunk_with_sink<S>(
        &mut self,
        config: FastqChunkConfig,
        sink: &mut S,
    ) -> Result<Option<FastqChunkStats>>
    where
        S: FastqRecordSink,
    {
        let config = config.normalized();
        let first_record_index = self.record_index;
        let mut chunk = FastqChunkStats::new(first_record_index);

        loop {
            let parse_start = if self.chunk_resume_start == 0 {
                self.compact_carry();
                self.fill_slab()?;
                0
            } else {
                let parse_start = self.chunk_resume_start;
                self.chunk_resume_start = 0;
                parse_start
            };

            let (next_start, records, bases, stopped) = visit_chunk_in_slab(
                &self.buf[parse_start..self.len],
                ChunkVisitContext {
                    eof: self.eof,
                    validate: self.config.validate,
                    base_offset: self.base_offset + parse_start as u64,
                    first_record_index: self.record_index,
                    config,
                },
                &chunk,
                sink,
            )?;
            let next_start = parse_start + next_start;
            if stopped && next_start < self.len {
                self.chunk_resume_start = next_start;
                self.next_start = 0;
            } else {
                self.next_start = next_start;
            }
            self.record_index += records;
            chunk.records += records;
            chunk.bases += bases;

            if chunk.records > 0 && (stopped || (self.len == 0 && self.eof)) {
                return Ok(Some(chunk));
            }
            if records == 0 {
                if self.len == 0 && self.eof {
                    return if chunk.records == 0 {
                        Ok(None)
                    } else {
                        Ok(Some(chunk))
                    };
                }
                return Err(FastqError::RecordTooLarge {
                    slab_size: self.config.slab_size,
                });
            }
        }
    }

    fn compact_carry(&mut self) {
        if self.chunk_resume_start != 0 {
            self.next_start = self.chunk_resume_start;
            self.chunk_resume_start = 0;
        }
        if self.next_start == 0 {
            return;
        }
        if self.next_start >= self.len {
            self.base_offset += self.len as u64;
            self.len = 0;
            self.next_start = 0;
            return;
        }
        let carry = self.len - self.next_start;
        self.buf.copy_within(self.next_start..self.len, 0);
        self.base_offset += self.next_start as u64;
        self.len = carry;
        self.next_start = 0;
    }

    fn fill_slab(&mut self) -> Result<()> {
        while !self.eof && self.len < self.config.slab_size {
            let n = self
                .reader
                .read(&mut self.buf[self.len..self.config.slab_size])?;
            if n == 0 {
                self.eof = true;
                break;
            }
            self.len += n;
        }
        Ok(())
    }

    fn align_interleaved_batch(&mut self, first_record_index: u64) -> Result<()> {
        if self.config.pairing != PairingMode::Interleaved || self.records.len().is_multiple_of(2) {
            return Ok(());
        }

        let Some(last) = self.records.last() else {
            return Ok(());
        };
        if self.eof {
            return Err(fastq_frame::format_at(
                "interleaved FASTQ ended with an unpaired record",
                self.base_offset,
                last.name.start as usize,
                first_record_index + (self.records.len() - 1) as u64,
                0,
            ));
        }

        self.next_start = last.name.start as usize;
        self.records.pop();
        Ok(())
    }

    fn retain_records_from_parts(&mut self, next_start: usize, record_count: usize) {
        self.chunk_resume_start = 0;
        self.next_start = next_start;
        self.record_index -= record_count as u64;
    }
}

/// Visit records from an already resident FASTQ byte slice.
///
/// This path is intended for memory-mapped files, cached datasets, and other
/// callers that already own a complete FASTQ byte buffer. It validates the same
/// record structure as [`FastqReader::visit_records`] when
/// [`FastqConfig::validate`] is enabled, but it does not copy the input into a
/// streaming slab and does not perform paired-end identifier validation.
///
/// Returns the number of visited records.
pub fn visit_fastq_bytes<F>(bytes: &[u8], config: FastqConfig, mut visit: F) -> Result<u64>
where
    F: FnMut(FastqVisitRecord<'_>) -> Result<()>,
{
    let (_, records) = visit_records_in_slab(bytes, true, config.validate, 0, 0, &mut visit)?;

    Ok(records)
}

/// Count records and bases from a FASTQ reader without building record views.
#[inline]
pub fn count_fastq_read<R: Read>(reader: R) -> Result<FastqStats> {
    count_fastq_read_with_config(reader, FastqConfig::default())
}

/// Count records and bases from a FASTQ reader with explicit parser settings.
#[inline]
pub fn count_fastq_read_with_config<R: Read>(reader: R, config: FastqConfig) -> Result<FastqStats> {
    let mut reader = FastqReader::with_config(reader, config);
    reader.count_records()
}

/// Count records and bases from an already resident FASTQ byte slice.
#[inline]
pub fn count_fastq_bytes(bytes: &[u8], config: FastqConfig) -> Result<FastqStats> {
    let mut stats = FastqStats::default();
    count_records_in_slab(bytes, true, config.validate, 0, 0, &mut stats)?;
    Ok(stats)
}

#[cfg(test)]
mod tests;
