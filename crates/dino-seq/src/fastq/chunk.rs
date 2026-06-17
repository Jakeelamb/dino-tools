use crate::Result;

use super::FastqVisitRecord;

/// Caller-provided sink for single-pass FASTQ record processing.
///
/// This trait is the allocation-control counterpart to
/// [`FastqReader::next_batch`](super::FastqReader::next_batch). Use it when a
/// downstream tool wants to fill its own output buffers directly and does not
/// need Dino Seq to build a batch side table.
///
/// If `record` returns an error, [`FastqReader::next_chunk_with_sink`](super::FastqReader::next_chunk_with_sink)
/// returns that error immediately. The sink may already have received records
/// from the current chunk, and the reader should be dropped rather than reused
/// after a sink error.
pub trait FastqRecordSink {
    /// Consume one borrowed FASTQ record.
    fn record(&mut self, record: FastqVisitRecord<'_>) -> Result<()>;
}

impl<F> FastqRecordSink for F
where
    F: for<'a> FnMut(FastqVisitRecord<'a>) -> Result<()>,
{
    fn record(&mut self, record: FastqVisitRecord<'_>) -> Result<()> {
        self(record)
    }
}

/// Optional extension point for sinks that can preallocate per chunk.
///
/// Dino Seq calls [`FastqRecordSink::record`] for correctness. Downstream tools
/// that own a growable output buffer can also implement this trait locally and
/// call [`FastqChunkConfig::estimated_records`] before invoking
/// [`FastqReader::next_chunk_with_sink`](super::FastqReader::next_chunk_with_sink).
pub trait FastqChunkSinkExt: FastqRecordSink {
    /// Reserve space for approximately `records` upcoming records.
    fn reserve_records(&mut self, _records: usize) {}
}

/// Chunking policy for [`FastqReader::next_chunk_with_sink`](super::FastqReader::next_chunk_with_sink).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FastqChunkConfig {
    /// Target number of sequence bases to emit before returning a chunk.
    ///
    /// A value of zero disables the target limit.
    pub target_bases: u64,
    /// Minimum number of records required before the target-base limit may stop
    /// a chunk.
    ///
    /// Values below one are raised to one.
    pub min_records: usize,
    /// Hard maximum number of sequence bases to emit before returning a chunk.
    ///
    /// A value of zero disables the hard limit. Unlike [`target_bases`](Self::target_bases),
    /// this limit is not gated by [`min_records`](Self::min_records).
    pub max_bases: u64,
}

impl FastqChunkConfig {
    /// Create a chunking policy with a target base count.
    pub fn new(target_bases: u64) -> Self {
        Self {
            target_bases,
            min_records: 1,
            max_bases: 0,
        }
    }

    /// Return a copy with a minimum record count.
    pub fn min_records(mut self, min_records: usize) -> Self {
        self.min_records = min_records.max(1);
        self
    }

    /// Return a copy with a hard maximum base count.
    pub fn max_bases(mut self, max_bases: u64) -> Self {
        self.max_bases = max_bases;
        self
    }

    /// Target bases per chunk. Zero means no target-base limit.
    pub fn target_bases(self) -> u64 {
        self.target_bases
    }

    /// Minimum records required before the target-base limit may stop a chunk.
    pub fn min_records_value(self) -> usize {
        self.min_records
    }

    /// Hard maximum bases per chunk. Zero means no hard limit.
    pub fn max_bases_value(self) -> u64 {
        self.max_bases
    }

    /// Estimate records per chunk for a fixed or representative read length.
    pub fn estimated_records(self, read_len: usize) -> usize {
        let read_len = read_len.max(1) as u64;
        let target = if self.max_bases != 0 {
            self.max_bases
        } else {
            self.target_bases
        };
        if target == 0 {
            return self.min_records.max(1);
        }
        target
            .div_ceil(read_len)
            .max(self.min_records.max(1) as u64)
            .min(usize::MAX as u64) as usize
    }

    pub(crate) fn normalized(self) -> Self {
        Self {
            min_records: self.min_records.max(1),
            ..self
        }
    }

    pub(crate) fn should_stop(self, records: u64, bases: u64) -> bool {
        if self.max_bases != 0 && bases >= self.max_bases {
            return true;
        }
        self.target_bases != 0 && bases >= self.target_bases && records >= self.min_records as u64
    }
}

/// Summary for one chunk emitted by [`FastqReader::next_chunk_with_sink`](super::FastqReader::next_chunk_with_sink).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FastqChunkStats {
    /// Zero-based index of the first record in this chunk.
    pub first_record_index: u64,
    /// Number of records emitted into the sink.
    pub records: u64,
    /// Number of sequence bases emitted into the sink.
    pub bases: u64,
}

impl FastqChunkStats {
    pub(crate) fn new(first_record_index: u64) -> Self {
        Self {
            first_record_index,
            records: 0,
            bases: 0,
        }
    }

    /// Zero-based index of the first record in this chunk.
    pub fn first_record_index(self) -> u64 {
        self.first_record_index
    }

    /// Number of records emitted into the sink.
    pub fn records(self) -> u64 {
        self.records
    }

    /// Number of sequence bases emitted into the sink.
    pub fn bases(self) -> u64 {
        self.bases
    }

    /// Whether the chunk contains no records.
    pub fn is_empty(self) -> bool {
        self.records == 0
    }
}
