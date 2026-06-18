use crate::error::{FastqError, Result};
use crate::fastq_frame;

use super::FastqStats;
use super::chunk::{FastqChunkConfig, FastqChunkStats, FastqRecordSink};
use super::record::{FastqVisitRecord, RecordRef, to_u32_range};

pub(super) fn frame_records(
    bytes: &[u8],
    eof: bool,
    validate: bool,
    base_offset: u64,
    first_record_index: u64,
    records: &mut Vec<RecordRef>,
) -> Result<usize> {
    if bytes.len() > u32::MAX as usize {
        return Err(FastqError::Format(
            "FASTQ slab byte offsets exceed u32 range".into(),
        ));
    }
    records.reserve(bytes.len() / 128);

    let mut cursor = 0;
    let mut newlines = memchr::memchr_iter(b'\n', bytes);

    while cursor < bytes.len() {
        let record_start = cursor;
        let Some(name) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return Ok(record_start);
        };
        let Some(seq) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return incomplete_or_truncated_frame(
                eof,
                base_offset,
                record_start,
                first_record_index + records.len() as u64,
                1,
            );
        };
        let Some(plus) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return incomplete_or_truncated_frame(
                eof,
                base_offset,
                record_start,
                first_record_index + records.len() as u64,
                2,
            );
        };
        let Some(qual) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return incomplete_or_truncated_frame(
                eof,
                base_offset,
                record_start,
                first_record_index + records.len() as u64,
                3,
            );
        };

        if validate {
            validate_record_ranges(
                bytes,
                name,
                seq,
                plus,
                qual,
                base_offset,
                first_record_index + records.len() as u64,
            )?;
        }

        records.push(RecordRef {
            name: to_u32_range(name.0..name.1),
            seq: to_u32_range(seq.0..seq.1),
            plus: to_u32_range(plus.0..plus.1),
            qual: to_u32_range(qual.0..qual.1),
        });
    }

    Ok(bytes.len())
}

fn incomplete_or_truncated_frame(
    eof: bool,
    base_offset: u64,
    record_start: usize,
    record_index: u64,
    line_index: u8,
) -> Result<usize> {
    if eof {
        return Err(fastq_frame::format_at(
            "truncated FASTQ record",
            base_offset,
            record_start,
            record_index,
            line_index,
        ));
    }
    Ok(record_start)
}

pub(super) fn visit_records_in_slab<F>(
    bytes: &[u8],
    eof: bool,
    validate: bool,
    base_offset: u64,
    first_record_index: u64,
    visit: &mut F,
) -> Result<(usize, u64)>
where
    F: FnMut(FastqVisitRecord<'_>) -> Result<()>,
{
    let mut cursor = 0;
    let mut records = 0_u64;
    let mut newlines = memchr::memchr_iter(b'\n', bytes);

    while cursor < bytes.len() {
        let record_start = cursor;
        let Some(name) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return Ok((record_start, records));
        };
        let Some(seq) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return incomplete_or_truncated_visit(
                eof,
                base_offset,
                record_start,
                first_record_index + records,
                records,
                1,
            );
        };
        let Some(plus) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return incomplete_or_truncated_visit(
                eof,
                base_offset,
                record_start,
                first_record_index + records,
                records,
                2,
            );
        };
        let Some(qual) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return incomplete_or_truncated_visit(
                eof,
                base_offset,
                record_start,
                first_record_index + records,
                records,
                3,
            );
        };

        if validate {
            let record_index = first_record_index + records;
            validate_record_ranges(bytes, name, seq, plus, qual, base_offset, record_index)?;
        }

        visit(FastqVisitRecord {
            name: &bytes[name.0..name.1],
            seq: &bytes[seq.0..seq.1],
            plus: &bytes[plus.0..plus.1],
            qual: &bytes[qual.0..qual.1],
        })?;
        records += 1;
    }

    Ok((bytes.len(), records))
}

#[inline]
pub(super) fn count_records_in_slab(
    bytes: &[u8],
    eof: bool,
    validate: bool,
    base_offset: u64,
    first_record_index: u64,
    stats: &mut FastqStats,
) -> Result<(usize, u64)> {
    if validate {
        return count_records_in_slab_impl::<true>(
            bytes,
            eof,
            base_offset,
            first_record_index,
            stats,
        );
    }
    count_records_in_slab_impl::<false>(bytes, eof, base_offset, first_record_index, stats)
}

#[inline]
fn count_records_in_slab_impl<const VALIDATE: bool>(
    bytes: &[u8],
    eof: bool,
    base_offset: u64,
    first_record_index: u64,
    stats: &mut FastqStats,
) -> Result<(usize, u64)> {
    let mut cursor = 0;
    let mut records = 0_u64;
    let mut newlines = memchr::memchr_iter(b'\n', bytes);

    while cursor < bytes.len() {
        let record_start = cursor;

        let name_start = cursor;
        let Some(name_raw_end) = newlines.next() else {
            if eof && name_start < bytes.len() {
                return incomplete_or_truncated_visit(
                    eof,
                    base_offset,
                    record_start,
                    first_record_index + records,
                    records,
                    1,
                );
            }
            return Ok((record_start, records));
        };
        let name_end = fastq_frame::trim_cr_end(bytes, name_start, name_raw_end);

        let seq_start = name_raw_end + 1;
        let Some(seq_raw_end) = newlines.next() else {
            if eof && seq_start < bytes.len() {
                return incomplete_or_truncated_visit(
                    eof,
                    base_offset,
                    record_start,
                    first_record_index + records,
                    records,
                    2,
                );
            }
            return incomplete_or_truncated_visit(
                eof,
                base_offset,
                record_start,
                first_record_index + records,
                records,
                1,
            );
        };
        let seq_end = fastq_frame::trim_cr_end(bytes, seq_start, seq_raw_end);

        let plus_start = seq_raw_end + 1;
        let Some(plus_raw_end) = newlines.next() else {
            if eof && plus_start < bytes.len() {
                return incomplete_or_truncated_visit(
                    eof,
                    base_offset,
                    record_start,
                    first_record_index + records,
                    records,
                    3,
                );
            }
            return incomplete_or_truncated_visit(
                eof,
                base_offset,
                record_start,
                first_record_index + records,
                records,
                2,
            );
        };
        let qual_start = plus_raw_end + 1;
        let qual_raw_end = if let Some(end) = newlines.next() {
            cursor = end + 1;
            end
        } else if eof && qual_start < bytes.len() {
            cursor = bytes.len();
            bytes.len()
        } else {
            return incomplete_or_truncated_visit(
                eof,
                base_offset,
                record_start,
                first_record_index + records,
                records,
                3,
            );
        };
        let qual_end = fastq_frame::trim_cr_end(bytes, qual_start, qual_raw_end);

        if VALIDATE {
            let record_index = first_record_index + records;
            validate_count_record(
                bytes,
                CountRecordRanges {
                    name_start,
                    seq_start,
                    seq_end,
                    plus_start,
                    qual_start,
                    qual_end,
                },
                base_offset,
                record_index,
            )?;
        }

        let seq_len = (seq_end - seq_start) as u64;
        let seq_first = if seq_start < seq_end {
            bytes[seq_start]
        } else {
            0
        };
        stats.records += 1;
        stats.bases += seq_len;
        stats.qualities += (qual_end - qual_start) as u64;
        stats.name_bytes += (name_end - name_start) as u64;
        stats.checksum = stats
            .checksum
            .wrapping_add(seq_first as u64)
            .wrapping_mul(1_099_511_628_211)
            .wrapping_add(seq_len);
        records += 1;
    }

    Ok((bytes.len(), records))
}

#[derive(Debug, Clone, Copy)]
struct CountRecordRanges {
    name_start: usize,
    seq_start: usize,
    seq_end: usize,
    plus_start: usize,
    qual_start: usize,
    qual_end: usize,
}

#[inline(always)]
fn validate_count_record(
    bytes: &[u8],
    ranges: CountRecordRanges,
    base_offset: u64,
    record_index: u64,
) -> Result<()> {
    if bytes.get(ranges.name_start) != Some(&b'@') {
        return Err(fastq_frame::format_at(
            "header must start with `@`",
            base_offset,
            ranges.name_start,
            record_index,
            0,
        ));
    }
    match bytes.get(ranges.name_start + 1) {
        Some(byte) if !byte.is_ascii_whitespace() => {}
        _ => {
            return Err(fastq_frame::format_at(
                "empty FASTQ id",
                base_offset,
                ranges.name_start,
                record_index,
                0,
            ));
        }
    }
    if bytes.get(ranges.plus_start) != Some(&b'+') {
        return Err(fastq_frame::format_at(
            "plus line must start with `+`",
            base_offset,
            ranges.plus_start,
            record_index,
            2,
        ));
    }
    let seq_len = ranges.seq_end - ranges.seq_start;
    let qual_len = ranges.qual_end - ranges.qual_start;
    if seq_len != qual_len {
        return Err(fastq_frame::format_at(
            format!("quality length {qual_len} != sequence length {seq_len}"),
            base_offset,
            ranges.qual_start,
            record_index,
            3,
        ));
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
pub(super) struct ChunkVisitContext {
    pub(super) eof: bool,
    pub(super) validate: bool,
    pub(super) base_offset: u64,
    pub(super) first_record_index: u64,
    pub(super) config: FastqChunkConfig,
}

pub(super) fn visit_chunk_in_slab<S>(
    bytes: &[u8],
    context: ChunkVisitContext,
    current: &FastqChunkStats,
    sink: &mut S,
) -> Result<(usize, u64, u64, bool)>
where
    S: FastqRecordSink,
{
    let mut cursor = 0;
    let mut records = 0_u64;
    let mut bases = 0_u64;
    let mut newlines = memchr::memchr_iter(b'\n', bytes);

    while cursor < bytes.len() {
        let record_start = cursor;
        let Some(name) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, context.eof)
        else {
            return Ok((record_start, records, bases, false));
        };
        let Some(seq) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, context.eof)
        else {
            let (next_start, prior_records) = incomplete_or_truncated_visit(
                context.eof,
                context.base_offset,
                record_start,
                context.first_record_index + records,
                records,
                1,
            )?;
            return Ok((next_start, prior_records, bases, false));
        };
        let Some(plus) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, context.eof)
        else {
            let (next_start, prior_records) = incomplete_or_truncated_visit(
                context.eof,
                context.base_offset,
                record_start,
                context.first_record_index + records,
                records,
                2,
            )?;
            return Ok((next_start, prior_records, bases, false));
        };
        let Some(qual) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, context.eof)
        else {
            let (next_start, prior_records) = incomplete_or_truncated_visit(
                context.eof,
                context.base_offset,
                record_start,
                context.first_record_index + records,
                records,
                3,
            )?;
            return Ok((next_start, prior_records, bases, false));
        };

        if context.validate {
            let record_index = context.first_record_index + records;
            validate_record_ranges(
                bytes,
                name,
                seq,
                plus,
                qual,
                context.base_offset,
                record_index,
            )?;
        }

        sink.record(FastqVisitRecord {
            name: &bytes[name.0..name.1],
            seq: &bytes[seq.0..seq.1],
            plus: &bytes[plus.0..plus.1],
            qual: &bytes[qual.0..qual.1],
        })?;
        records += 1;
        bases += (seq.1 - seq.0) as u64;

        if context
            .config
            .should_stop(current.records + records, current.bases + bases)
        {
            return Ok((cursor, records, bases, true));
        }
    }

    Ok((bytes.len(), records, bases, false))
}

#[inline(always)]
fn next_visit_line_bounds(
    bytes: &[u8],
    cursor: &mut usize,
    newlines: &mut impl Iterator<Item = usize>,
    eof: bool,
) -> Option<(usize, usize)> {
    let start = *cursor;
    if let Some(end) = newlines.next() {
        *cursor = end + 1;
        let end = fastq_frame::trim_cr_end(bytes, start, end);
        return Some((start, end));
    }
    if eof && start < bytes.len() {
        *cursor = bytes.len();
        let end = fastq_frame::trim_cr_end(bytes, start, bytes.len());
        return Some((start, end));
    }
    None
}

#[inline(always)]
fn validate_record_ranges(
    bytes: &[u8],
    name: (usize, usize),
    seq: (usize, usize),
    plus: (usize, usize),
    qual: (usize, usize),
    base_offset: u64,
    record_index: u64,
) -> Result<()> {
    if bytes.get(name.0) != Some(&b'@') {
        return Err(fastq_frame::format_at(
            "header must start with `@`",
            base_offset,
            name.0,
            record_index,
            0,
        ));
    }
    match bytes.get(name.0 + 1) {
        Some(byte) if !byte.is_ascii_whitespace() => {}
        _ => {
            return Err(fastq_frame::format_at(
                "empty FASTQ id",
                base_offset,
                name.0,
                record_index,
                0,
            ));
        }
    }
    if bytes.get(plus.0) != Some(&b'+') {
        return Err(fastq_frame::format_at(
            "plus line must start with `+`",
            base_offset,
            plus.0,
            record_index,
            2,
        ));
    }
    let seq_len = seq.1 - seq.0;
    let qual_len = qual.1 - qual.0;
    if seq_len != qual_len {
        return Err(fastq_frame::format_at(
            format!("quality length {qual_len} != sequence length {seq_len}"),
            base_offset,
            qual.0,
            record_index,
            3,
        ));
    }
    Ok(())
}

fn incomplete_or_truncated_visit(
    eof: bool,
    base_offset: u64,
    record_start: usize,
    record_index: u64,
    records: u64,
    line_index: u8,
) -> Result<(usize, u64)> {
    if eof {
        return Err(fastq_frame::format_at(
            "truncated FASTQ record",
            base_offset,
            record_start,
            record_index,
            line_index,
        ));
    }
    Ok((record_start, records))
}
