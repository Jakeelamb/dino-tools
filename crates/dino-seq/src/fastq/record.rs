use std::ops::Range;

use crate::fastq_frame;

/// Byte ranges for a four-line FASTQ record within a batch slab.
///
/// Ranges are local to [`FastqBatch::bytes`](super::FastqBatch::bytes). They
/// are exposed for callers that want to build their own zero-copy views without
/// using [`FastqRecord`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordRef {
    /// Header line, including the leading `@` and excluding the newline.
    pub name: Range<u32>,
    /// Sequence line, excluding the newline.
    pub seq: Range<u32>,
    /// Plus line, excluding the newline.
    pub plus: Range<u32>,
    /// Quality line, excluding the newline.
    pub qual: Range<u32>,
}

/// Borrowed view of one FASTQ record inside a [`FastqBatch`](super::FastqBatch).
///
/// The record borrows from the batch's slab buffer and cannot outlive the
/// batch. Accessors return byte slices rather than UTF-8 strings because FASTQ
/// names and qualities are byte-oriented.
#[derive(Debug, Clone, Copy)]
pub struct FastqRecord<'a> {
    pub(crate) bytes: &'a [u8],
    pub(crate) record: &'a RecordRef,
}

impl<'a> FastqRecord<'a> {
    /// Return the header line including the leading `@`.
    #[inline]
    pub fn name(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.name.clone())]
    }

    /// Return the header line without a leading `@`.
    #[inline]
    pub fn name_without_at(self) -> &'a [u8] {
        let name = self.name();
        name.strip_prefix(b"@").unwrap_or(name)
    }

    /// Return the first whitespace-delimited identifier token.
    #[inline]
    pub fn id_token(self) -> &'a [u8] {
        let name = self.name_without_at();
        let end = name
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(name.len());
        &name[..end]
    }

    /// Return the identifier token with a trailing `/1` or `/2` removed.
    #[inline]
    pub fn pair_normalized_id(self) -> &'a [u8] {
        strip_pair_suffix(self.id_token())
    }

    /// Return the sequence line.
    #[inline]
    pub fn seq(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.seq.clone())]
    }

    /// Return the plus line.
    #[inline]
    pub fn plus(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.plus.clone())]
    }

    /// Return the quality line.
    #[inline]
    pub fn qual(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.qual.clone())]
    }
}

/// Borrowed FASTQ record passed to [`FastqReader::visit_records`](super::FastqReader::visit_records).
///
/// This view is optimized for single-pass consumers that do not need the batch
/// side table. The slices point directly into the reader slab and are valid
/// only for the duration of the visitor callback.
#[derive(Debug, Clone, Copy)]
pub struct FastqVisitRecord<'a> {
    pub(crate) name: &'a [u8],
    pub(crate) seq: &'a [u8],
    pub(crate) plus: &'a [u8],
    pub(crate) qual: &'a [u8],
}

impl<'a> FastqVisitRecord<'a> {
    /// Return the header line including the leading `@`.
    #[inline]
    pub fn name(self) -> &'a [u8] {
        self.name
    }

    /// Return the sequence line.
    #[inline]
    pub fn seq(self) -> &'a [u8] {
        self.seq
    }

    /// Return the plus line.
    #[inline]
    pub fn plus(self) -> &'a [u8] {
        self.plus
    }

    /// Return the quality line.
    #[inline]
    pub fn qual(self) -> &'a [u8] {
        self.qual
    }
}

/// Strip a terminal `/1` or `/2` pair suffix from an identifier token.
#[inline]
pub fn strip_pair_suffix(id: &[u8]) -> &[u8] {
    fastq_frame::strip_pair_suffix(id)
}

#[inline]
pub(crate) fn to_u32_range(range: Range<usize>) -> Range<u32> {
    debug_assert!(u32::try_from(range.end).is_ok());
    range.start as u32..range.end as u32
}

#[inline]
fn to_usize(range: Range<u32>) -> Range<usize> {
    range.start as usize..range.end as usize
}
