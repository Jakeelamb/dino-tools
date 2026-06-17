use std::io::Read;

use crate::error::{FastqError, Result};
use crate::fastq_frame;

use super::record::{FastqRecord, RecordRef};
use super::{FastqBatch, FastqConfig, FastqReader};

/// Pairing interpretation for a single FASTQ stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairingMode {
    /// Records are yielded independently.
    None,
    /// Adjacent records are expected to be ordered mates.
    Interleaved,
}

/// Identifier validation policy for paired-end data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairValidation {
    /// Compare normalized first-token identifiers after removing `/1` or `/2`.
    Full,
    /// Fast path for ordered `/1` and `/2` mate suffixes, with fallback to
    /// [`Full`](Self::Full) when that shape is not present.
    FastSlash,
    /// Skip identifier checks for trusted, already synchronized inputs.
    None,
}

/// Borrowed view of an ordered read pair.
#[derive(Debug, Clone, Copy)]
pub struct FastqPair<'a> {
    pub(crate) first: FastqRecord<'a>,
    pub(crate) second: FastqRecord<'a>,
}

impl<'a> FastqPair<'a> {
    /// First mate.
    pub fn first(&self) -> FastqRecord<'a> {
        self.first
    }

    /// Second mate.
    pub fn second(&self) -> FastqRecord<'a> {
        self.second
    }

    /// Normalized pair identifier from the first mate.
    pub fn pair_id(&self) -> &'a [u8] {
        self.first.pair_normalized_id()
    }
}

/// Iterator over adjacent pairs in an interleaved [`FastqBatch`].
#[derive(Debug)]
pub struct InterleavedPairs<'a> {
    pub(crate) batch: &'a FastqBatch<'a>,
    pub(crate) next: usize,
}

impl<'a> Iterator for InterleavedPairs<'a> {
    type Item = FastqPair<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.batch.records.len() {
            return None;
        }
        let pair = FastqPair {
            first: self.batch.record_at(self.next),
            second: self.batch.record_at(self.next + 1),
        };
        self.next += 2;
        Some(pair)
    }
}

/// Iterator over paired records from two separate [`FastqBatch`] values.
#[derive(Debug)]
pub struct PairedRecords<'a> {
    pub(crate) first: &'a FastqBatch<'a>,
    pub(crate) second: &'a FastqBatch<'a>,
    pub(crate) next: usize,
}

impl<'a> Iterator for PairedRecords<'a> {
    type Item = FastqPair<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.first.records.len() {
            return None;
        }
        let pair = FastqPair {
            first: self.first.record_at(self.next),
            second: self.second.record_at(self.next),
        };
        self.next += 1;
        Some(pair)
    }
}

/// Validate and zip two separate FASTQ batches by ordered mate identifiers.
pub fn paired_records<'a>(
    first: &'a FastqBatch<'a>,
    second: &'a FastqBatch<'a>,
) -> Result<PairedRecords<'a>> {
    first.paired_with(second)
}

/// A paired-end batch produced by [`PairedFastqReader`].
///
/// The batch contains the validated common prefix from the two underlying
/// readers. If one reader yields more records than the other in a slab, the
/// extra records are retained for the next call.
#[derive(Debug)]
pub struct PairedFastqBatch<'a> {
    pub(crate) first_bytes: &'a [u8],
    pub(crate) first_records: &'a [RecordRef],
    pub(crate) second_bytes: &'a [u8],
    pub(crate) second_records: &'a [RecordRef],
}

impl<'a> PairedFastqBatch<'a> {
    /// Number of read pairs in the batch.
    pub fn len(&self) -> usize {
        self.first_records.len()
    }

    /// Whether this batch contains no read pairs.
    pub fn is_empty(&self) -> bool {
        self.first_records.is_empty()
    }

    /// Iterate borrowed read-pair views.
    pub fn pairs(&'a self) -> PairedFastqPairs<'a> {
        PairedFastqPairs {
            batch: self,
            next: 0,
        }
    }

    fn first_record_at(&self, index: usize) -> FastqRecord<'a> {
        FastqRecord {
            bytes: self.first_bytes,
            record: &self.first_records[index],
        }
    }

    fn second_record_at(&self, index: usize) -> FastqRecord<'a> {
        FastqRecord {
            bytes: self.second_bytes,
            record: &self.second_records[index],
        }
    }
}

/// Iterator over read pairs in a [`PairedFastqBatch`].
#[derive(Debug)]
pub struct PairedFastqPairs<'a> {
    batch: &'a PairedFastqBatch<'a>,
    next: usize,
}

impl<'a> Iterator for PairedFastqPairs<'a> {
    type Item = FastqPair<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.batch.first_records.len() {
            return None;
        }
        let pair = FastqPair {
            first: self.batch.first_record_at(self.next),
            second: self.batch.second_record_at(self.next),
        };
        self.next += 1;
        Some(pair)
    }
}

/// Stateful reader for ordered separate-file paired-end FASTQ streams.
///
/// This reader validates matching records as it advances both streams. It does
/// not reorder or synchronize mates that appear in different orders.
#[derive(Debug)]
pub struct PairedFastqReader<R1, R2> {
    first: FastqReader<R1>,
    second: FastqReader<R2>,
    pair_validation: PairValidation,
    first_pending_retain: Option<(usize, usize)>,
    second_pending_retain: Option<(usize, usize)>,
}

impl<R1: Read, R2: Read> PairedFastqReader<R1, R2> {
    /// Create a paired reader with default FASTQ configuration.
    pub fn new(first: R1, second: R2) -> Self {
        Self::with_config(first, second, FastqConfig::default())
    }

    /// Create a paired reader using the same configuration for both streams.
    pub fn with_config(first: R1, second: R2, config: FastqConfig) -> Self {
        Self::with_configs(
            first,
            FastqConfig {
                pairing: PairingMode::None,
                ..config.clone()
            },
            second,
            FastqConfig {
                pairing: PairingMode::None,
                ..config
            },
        )
    }

    /// Create a paired reader with separate per-stream configurations.
    pub fn with_configs(
        first: R1,
        first_config: FastqConfig,
        second: R2,
        second_config: FastqConfig,
    ) -> Self {
        let pair_validation = first_config.pair_validation;
        Self {
            first: FastqReader::with_config(
                first,
                FastqConfig {
                    pairing: PairingMode::None,
                    ..first_config
                },
            ),
            second: FastqReader::with_config(
                second,
                FastqConfig {
                    pairing: PairingMode::None,
                    ..second_config
                },
            ),
            pair_validation,
            first_pending_retain: None,
            second_pending_retain: None,
        }
    }

    /// Create a paired reader from two existing [`FastqReader`] values.
    pub fn from_fastq_readers(first: FastqReader<R1>, second: FastqReader<R2>) -> Self {
        Self {
            pair_validation: first.config.pair_validation,
            first,
            second,
            first_pending_retain: None,
            second_pending_retain: None,
        }
    }

    /// Override the paired identifier validation policy.
    pub fn with_pair_validation(mut self, pair_validation: PairValidation) -> Self {
        self.pair_validation = pair_validation;
        self
    }

    /// Read the next validated paired batch.
    ///
    /// Returns `Ok(None)` when both streams end together.
    pub fn next_pair_batch(&mut self) -> Result<Option<PairedFastqBatch<'_>>> {
        if let Some((next_start, record_count)) = self.first_pending_retain.take() {
            self.first
                .retain_records_from_parts(next_start, record_count);
        }
        if let Some((next_start, record_count)) = self.second_pending_retain.take() {
            self.second
                .retain_records_from_parts(next_start, record_count);
        }

        let first = self.first.next_batch()?;
        let second = self.second.next_batch()?;

        match (first, second) {
            (None, None) => Ok(None),
            (Some(first), None) => Err(extra_record_error(&first, 0)),
            (None, Some(second)) => Err(extra_record_error(&second, 0)),
            (Some(first), Some(second)) => {
                let pair_count = first.len().min(second.len());
                validate_pair_ids_prefix(&first, &second, pair_count, self.pair_validation)?;

                self.first_pending_retain = retain_from(&first, pair_count);
                self.second_pending_retain = retain_from(&second, pair_count);

                Ok(Some(PairedFastqBatch {
                    first_bytes: first.bytes,
                    first_records: &first.records[..pair_count],
                    second_bytes: second.bytes,
                    second_records: &second.records[..pair_count],
                }))
            }
        }
    }
}

fn retain_from(batch: &FastqBatch<'_>, index: usize) -> Option<(usize, usize)> {
    if index >= batch.records.len() {
        return None;
    }
    Some((
        batch.records[index].name.start as usize,
        batch.records.len() - index,
    ))
}

pub(crate) fn validate_even_pair_count(batch: &FastqBatch<'_>) -> Result<()> {
    if batch.records.len().is_multiple_of(2) {
        return Ok(());
    }
    let Some(last) = batch.records.last() else {
        return Ok(());
    };
    Err(fastq_frame::format_at(
        "interleaved FASTQ batch has an odd record count",
        batch.base_offset,
        last.name.start as usize,
        batch.first_record_index + (batch.records.len() - 1) as u64,
        0,
    ))
}

pub(crate) fn validate_paired_batches(
    first: &FastqBatch<'_>,
    second: &FastqBatch<'_>,
    pair_validation: PairValidation,
) -> Result<()> {
    if first.records.len() == second.records.len() {
        return validate_pair_ids(first, second, pair_validation);
    }

    let (batch, index) = if first.records.len() > second.records.len() {
        (first, second.records.len())
    } else {
        (second, first.records.len())
    };
    let record = &batch.records[index];
    Err(fastq_frame::format_at(
        "paired FASTQ batches have different record counts",
        batch.base_offset,
        record.name.start as usize,
        batch.first_record_index + index as u64,
        0,
    ))
}

fn validate_pair_ids(
    first: &FastqBatch<'_>,
    second: &FastqBatch<'_>,
    pair_validation: PairValidation,
) -> Result<()> {
    if pair_validation == PairValidation::None {
        return Ok(());
    }
    for index in 0..first.records.len() {
        let r1 = first.record_at(index);
        let r2 = second.record_at(index);
        if !pair_ids_match(r1, r2, pair_validation) {
            return Err(pair_id_mismatch(second, index));
        }
    }
    Ok(())
}

pub(crate) fn validate_interleaved_pair_ids(
    batch: &FastqBatch<'_>,
    pair_validation: PairValidation,
) -> Result<()> {
    if pair_validation == PairValidation::None {
        return Ok(());
    }
    for index in (0..batch.records.len()).step_by(2) {
        let r1 = batch.record_at(index);
        let r2 = batch.record_at(index + 1);
        if !pair_ids_match(r1, r2, pair_validation) {
            return Err(pair_id_mismatch(batch, index + 1));
        }
    }
    Ok(())
}

fn validate_pair_ids_prefix(
    first: &FastqBatch<'_>,
    second: &FastqBatch<'_>,
    len: usize,
    pair_validation: PairValidation,
) -> Result<()> {
    if pair_validation == PairValidation::None {
        return Ok(());
    }
    for index in 0..len {
        let r1 = first.record_at(index);
        let r2 = second.record_at(index);
        if !pair_ids_match(r1, r2, pair_validation) {
            return Err(pair_id_mismatch(second, index));
        }
    }
    Ok(())
}

fn pair_ids_match(first: FastqRecord<'_>, second: FastqRecord<'_>, mode: PairValidation) -> bool {
    match mode {
        PairValidation::None => true,
        PairValidation::Full => first.pair_normalized_id() == second.pair_normalized_id(),
        PairValidation::FastSlash => {
            fastq_frame::fast_slash_pair_ids_match(first.name(), second.name())
                .unwrap_or_else(|| first.pair_normalized_id() == second.pair_normalized_id())
        }
    }
}

fn pair_id_mismatch(batch: &FastqBatch<'_>, index: usize) -> FastqError {
    let record = &batch.records[index];
    fastq_frame::format_at(
        "paired FASTQ record identifiers do not match",
        batch.base_offset,
        record.name.start as usize,
        batch.first_record_index + index as u64,
        0,
    )
}

pub(crate) fn extra_record_error(batch: &FastqBatch<'_>, index: usize) -> FastqError {
    let record = &batch.records[index];
    fastq_frame::format_at(
        "paired FASTQ inputs have different record counts",
        batch.base_offset,
        record.name.start as usize,
        batch.first_record_index + index as u64,
        0,
    )
}
