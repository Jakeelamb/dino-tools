use crate::error::{FastqError, Result};
use crate::fastq_frame;

use super::chunk::{FastqChunkConfig, FastqChunkStats, FastqRecordSink};
use super::record::{FastqVisitRecord, RecordRef, to_u32_range};

pub(super) fn frame_records(
    bytes: &[u8],
    newline_offsets: &[usize],
    eof: bool,
    validate: bool,
    base_offset: u64,
    first_record_index: u64,
    records: &mut Vec<RecordRef>,
) -> Result<usize> {
    let layout = fastq_frame::slab_line_layout(bytes, newline_offsets, eof);
    if eof && !layout.line_count.is_multiple_of(4) {
        let start = fastq_frame::line_start(newline_offsets, (layout.line_count / 4) * 4);
        return Err(fastq_frame::format_at(
            "truncated FASTQ record",
            base_offset,
            start,
            first_record_index + records.len() as u64,
            (layout.line_count % 4) as u8,
        ));
    }
    if bytes.len() > u32::MAX as usize {
        return Err(FastqError::Format(
            "FASTQ slab byte offsets exceed u32 range".into(),
        ));
    }
    records.reserve(layout.complete_lines / 4);

    let mut line = 0;
    let mut name_start = 0;
    while line < layout.complete_lines {
        debug_assert!(line + 2 < newline_offsets.len());
        let name_end = fastq_frame::trim_cr_end(bytes, name_start, newline_offsets[line]);
        let seq_start = newline_offsets[line] + 1;
        let seq_end = fastq_frame::trim_cr_end(bytes, seq_start, newline_offsets[line + 1]);
        let plus_start = newline_offsets[line + 1] + 1;
        let plus_end = fastq_frame::trim_cr_end(bytes, plus_start, newline_offsets[line + 2]);
        let qual_start = newline_offsets[line + 2] + 1;
        let qual_line = line + 3;
        let qual_raw_end = if qual_line < newline_offsets.len() {
            newline_offsets[qual_line]
        } else {
            bytes.len()
        };
        let qual_end = fastq_frame::trim_cr_end(bytes, qual_start, qual_raw_end);
        let name = (name_start, name_end);
        let seq = (seq_start, seq_end);
        let plus = (plus_start, plus_end);
        let qual = (qual_start, qual_end);

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
        line += 4;
        if line < layout.complete_lines {
            name_start = newline_offsets[line - 1] + 1;
        }
    }

    if layout.complete_lines == layout.line_count && !layout.has_partial_trailing_line {
        Ok(bytes.len())
    } else {
        Ok(fastq_frame::line_start(
            newline_offsets,
            layout.complete_lines,
        ))
    }
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
