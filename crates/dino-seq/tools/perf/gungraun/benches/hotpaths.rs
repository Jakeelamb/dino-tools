use dino_seq::pack::{
    TrustedPackedRecord, pack_bases_and_summarize_qualities_into, pack_trusted_fastq_direct,
    pack_trusted_fastq_sink,
};
use dino_seq::{FastqConfig, FastqReader, FastqVisitRecord};
use gungraun::prelude::*;
use std::hint::black_box;

#[derive(Clone)]
struct Fixture {
    bytes: Vec<u8>,
}

#[derive(Clone)]
struct ReadFixture {
    seq: Vec<u8>,
    qual: Vec<u8>,
}

#[derive(Clone, Copy)]
enum FastqShape {
    CanonicalFixed,
    CanonicalVariable,
    AmbiguousFixed,
}

#[derive(Default)]
struct Stats {
    checksum: u64,
}

impl Stats {
    fn observe(&mut self, name: &[u8], seq: &[u8], qual: &[u8]) {
        self.checksum = self
            .checksum
            .wrapping_add(name.len() as u64)
            .wrapping_mul(1_099_511_628_211)
            .wrapping_add(seq.len() as u64)
            .wrapping_add(qual.len() as u64);
    }
}

fn synthetic_fastq_fixed(args: (usize, usize)) -> Fixture {
    synthetic_fastq(args, FastqShape::CanonicalFixed)
}

fn synthetic_fastq_variable(args: (usize, usize)) -> Fixture {
    synthetic_fastq(args, FastqShape::CanonicalVariable)
}

fn synthetic_fastq_ambiguous(args: (usize, usize)) -> Fixture {
    synthetic_fastq(args, FastqShape::AmbiguousFixed)
}

fn synthetic_fastq(args: (usize, usize), shape: FastqShape) -> Fixture {
    let (records, read_len) = args;
    let mut bytes = Vec::with_capacity(records * (read_len * 2 + 32));
    for i in 0..records {
        bytes.extend_from_slice(b"@r");
        bytes.extend_from_slice(i.to_string().as_bytes());
        bytes.push(b'\n');
        for j in 0..read_len {
            bytes.push(base_at(i, j, shape));
        }
        bytes.extend_from_slice(b"\n+\n");
        for j in 0..read_len {
            bytes.push(quality_at(i, j, shape));
        }
        bytes.push(b'\n');
    }
    Fixture { bytes }
}

fn read_canonical_fixed(len: usize) -> ReadFixture {
    synthetic_read(len, FastqShape::CanonicalFixed)
}

fn read_canonical_variable(len: usize) -> ReadFixture {
    synthetic_read(len, FastqShape::CanonicalVariable)
}

fn read_ambiguous_fixed(len: usize) -> ReadFixture {
    synthetic_read(len, FastqShape::AmbiguousFixed)
}

fn synthetic_read(len: usize, shape: FastqShape) -> ReadFixture {
    let mut seq = Vec::with_capacity(len);
    let mut qual = Vec::with_capacity(len);
    for i in 0..len {
        seq.push(base_at(0, i, shape));
        qual.push(quality_at(0, i, shape));
    }
    ReadFixture { seq, qual }
}

fn base_at(record: usize, offset: usize, shape: FastqShape) -> u8 {
    if matches!(shape, FastqShape::AmbiguousFixed) && (record + offset).is_multiple_of(17) {
        b'N'
    } else {
        b"ACGT"[(record + offset) & 3]
    }
}

fn quality_at(record: usize, offset: usize, shape: FastqShape) -> u8 {
    if matches!(shape, FastqShape::CanonicalVariable) {
        33 + ((record + offset) % 41) as u8
    } else {
        b'I'
    }
}

#[library_benchmark]
#[bench::canonical_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_fixed)]
#[bench::canonical_variable_10k(args = ((10_000, 150)), setup = synthetic_fastq_variable)]
#[bench::ambiguous_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_ambiguous)]
fn visit_records_strict(fixture: Fixture) -> u64 {
    visit_records::<true>(fixture)
}

#[library_benchmark]
#[bench::canonical_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_fixed)]
#[bench::canonical_variable_10k(args = ((10_000, 150)), setup = synthetic_fastq_variable)]
#[bench::ambiguous_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_ambiguous)]
fn visit_records_no_validate(fixture: Fixture) -> u64 {
    visit_records::<false>(fixture)
}

#[inline(always)]
fn visit_records<const VALIDATE: bool>(fixture: Fixture) -> u64 {
    let mut reader = FastqReader::with_config(
        fixture.bytes.as_slice(),
        FastqConfig {
            slab_size: 8 * 1024 * 1024,
            validate: VALIDATE,
            ..FastqConfig::default()
        },
    );
    let mut stats = Stats::default();
    if reader
        .visit_records(|record: FastqVisitRecord<'_>| {
            stats.observe(record.name(), record.seq(), record.qual());
            Ok(())
        })
        .is_err()
    {
        return u64::MAX;
    }
    black_box(stats.checksum)
}

#[library_benchmark]
#[bench::canonical_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_fixed)]
#[bench::canonical_variable_10k(args = ((10_000, 150)), setup = synthetic_fastq_variable)]
#[bench::ambiguous_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_ambiguous)]
fn parse_batches_strict(fixture: Fixture) -> u64 {
    parse_batches::<true>(fixture)
}

#[library_benchmark]
#[bench::canonical_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_fixed)]
#[bench::canonical_variable_10k(args = ((10_000, 150)), setup = synthetic_fastq_variable)]
#[bench::ambiguous_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_ambiguous)]
fn parse_batches_no_validate(fixture: Fixture) -> u64 {
    parse_batches::<false>(fixture)
}

#[inline(always)]
fn parse_batches<const VALIDATE: bool>(fixture: Fixture) -> u64 {
    let mut reader = FastqReader::with_config(
        fixture.bytes.as_slice(),
        FastqConfig {
            slab_size: 8 * 1024 * 1024,
            validate: VALIDATE,
            ..FastqConfig::default()
        },
    );
    let mut stats = Stats::default();
    loop {
        let batch = match reader.next_batch() {
            Ok(Some(batch)) => batch,
            Ok(None) => break,
            Err(_) => return u64::MAX,
        };
        for record in batch.records() {
            stats.observe(record.name(), record.seq(), record.qual());
        }
    }
    black_box(stats.checksum)
}

#[library_benchmark]
#[bench::canonical_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_fixed)]
#[bench::canonical_variable_10k(args = ((10_000, 150)), setup = synthetic_fastq_variable)]
#[bench::ambiguous_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_ambiguous)]
fn parse_and_pack(fixture: Fixture) -> u64 {
    let mut reader = FastqReader::new(fixture.bytes.as_slice());
    let mut packed = Vec::new();
    let mut mask = Vec::new();
    let mut checksum = 0_u64;
    loop {
        let batch = match reader.next_batch() {
            Ok(Some(batch)) => batch,
            Ok(None) => break,
            Err(_) => return u64::MAX,
        };
        for record in batch.records() {
            let summary = match pack_bases_and_summarize_qualities_into(
                record.seq(),
                record.qual(),
                &mut packed,
                &mut mask,
            ) {
                Ok(summary) => summary,
                Err(_) => return u64::MAX,
            };
            checksum = checksum
                .wrapping_add(summary.bases.canonical_bases() as u64)
                .wrapping_add(summary.qualities.sum_phred);
        }
    }
    black_box(checksum)
}

#[library_benchmark]
#[bench::canonical_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_fixed)]
#[bench::canonical_variable_10k(args = ((10_000, 150)), setup = synthetic_fastq_variable)]
#[bench::ambiguous_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_ambiguous)]
fn trusted_default_pack(fixture: Fixture) -> u64 {
    let mut checksum = 0_u64;
    if pack_trusted_fastq_sink(
        fixture.bytes.as_slice(),
        |record: TrustedPackedRecord<'_>| {
            checksum = checksum
                .wrapping_add(record.summary.bases.canonical_bases() as u64)
                .wrapping_add(record.summary.bases.n as u64)
                .wrapping_add(record.summary.qualities.sum_phred);
            Ok(())
        },
    )
    .is_err()
    {
        return u64::MAX;
    }
    black_box(checksum)
}

#[library_benchmark]
#[bench::canonical_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_fixed)]
#[bench::canonical_variable_10k(args = ((10_000, 150)), setup = synthetic_fastq_variable)]
#[bench::ambiguous_fixed_10k(args = ((10_000, 150)), setup = synthetic_fastq_ambiguous)]
fn trusted_direct_pack(fixture: Fixture) -> u64 {
    let mut checksum = 0_u64;
    if pack_trusted_fastq_direct(fixture.bytes.as_slice(), |record| {
        checksum = checksum
            .wrapping_add(record.summary.bases.canonical_bases() as u64)
            .wrapping_add(record.summary.bases.n as u64)
            .wrapping_add(record.summary.qualities.sum_phred);
        Ok(())
    })
    .is_err()
    {
        return u64::MAX;
    }
    black_box(checksum)
}

#[library_benchmark]
#[bench::canonical_fixed_150k(args = (150_000), setup = read_canonical_fixed)]
#[bench::canonical_variable_150k(args = (150_000), setup = read_canonical_variable)]
#[bench::ambiguous_fixed_150k(args = (150_000), setup = read_ambiguous_fixed)]
fn pack_large_read(fixture: ReadFixture) -> u64 {
    let mut packed = Vec::new();
    let mut mask = Vec::new();
    let summary = match pack_bases_and_summarize_qualities_into(
        &fixture.seq,
        &fixture.qual,
        &mut packed,
        &mut mask,
    ) {
        Ok(summary) => summary,
        Err(_) => return u64::MAX,
    };
    black_box(
        summary.bases.canonical_bases() as u64
            + summary.bases.n as u64
            + summary.qualities.sum_phred
            + packed.len() as u64
            + mask.len() as u64,
    )
}

library_benchmark_group!(
    name = hotpaths;
    benchmarks =
        visit_records_strict,
        visit_records_no_validate,
        parse_batches_strict,
        parse_batches_no_validate,
        parse_and_pack,
        trusted_default_pack,
        trusted_direct_pack,
        pack_large_read
);

main!(library_benchmark_groups = hotpaths);
