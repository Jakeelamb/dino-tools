use dino_seq::pack::{pack_bases_and_summarize_qualities_into, pack_trusted_fastq};
use dino_seq::{
    FastaIndex, FastaReader, FastqConfig, FastqReader, FastqVisitRecord, IndexedFastaReader,
    build_fasta_index, count_fasta_bytes,
};
use std::env;
use std::hint::black_box;
use std::io::Cursor;
use std::sync::Arc;
use std::time::Instant;

const RECORDS: usize = 50_000;
const READ_LEN: usize = 150;
const ITERS: usize = 10;
const DEFAULT_SLAB_SIZE: usize = 8 * 1024 * 1024;

type BenchResult<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Default)]
struct StreamStats {
    records: u64,
    bases: u64,
    qualities: u64,
    checksum: u64,
}

impl StreamStats {
    fn observe_record(&mut self, name: &[u8], seq: &[u8], qual: &[u8]) {
        self.records += 1;
        self.bases += seq.len() as u64;
        self.qualities += qual.len() as u64;
        self.checksum = self
            .checksum
            .wrapping_add(name.len() as u64)
            .wrapping_mul(1_099_511_628_211)
            .wrapping_add(seq.len() as u64)
            .wrapping_add(qual.len() as u64);
    }
}

fn consume_fastq<R: std::io::Read>(reader: &mut FastqReader<R>) -> dino_seq::Result<StreamStats> {
    let mut stats = StreamStats::default();
    while let Some(batch) = reader.next_batch()? {
        for record in batch.records() {
            stats.observe_record(record.name(), record.seq(), record.qual());
        }
    }
    Ok(stats)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Case {
    All,
    ParseStrict,
    VisitStrict,
    VisitNoValidate,
    ParseNoValidate,
    Pack,
    TrustedPack,
    FastaCountResident,
    FastaVisitStream,
    FastaIndexBuild,
    FastaFetchIndexed,
    #[cfg(feature = "bgzf")]
    BgzfDecodeParse,
}

impl Case {
    fn name(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::ParseStrict => "parse_raw_strict_8m_slab",
            Self::VisitStrict => "visit_raw_strict_8m_slab",
            Self::VisitNoValidate => "visit_raw_no_validate_8m_slab",
            Self::ParseNoValidate => "parse_raw_no_validate_8m_slab",
            Self::Pack => "parse_and_pack_seq_qual",
            Self::TrustedPack => "trusted_default_pack_seq_qual",
            Self::FastaCountResident => "fasta_count_resident_wrapped",
            Self::FastaVisitStream => "fasta_visit_stream_wrapped",
            Self::FastaIndexBuild => "fasta_index_build_wrapped",
            Self::FastaFetchIndexed => "fasta_fetch_indexed_wrapped",
            #[cfg(feature = "bgzf")]
            Self::BgzfDecodeParse => "bgzf_parallel_stream_parse",
        }
    }

    fn parse(value: &str) -> Option<Self> {
        match value {
            "all" => Some(Self::All),
            "parse" | "parse-strict" | "parse_raw_strict_8m_slab" => Some(Self::ParseStrict),
            "visit" | "visit-strict" | "visit_raw_strict_8m_slab" => Some(Self::VisitStrict),
            "visit-no-validate" | "visit_raw_no_validate_8m_slab" => Some(Self::VisitNoValidate),
            "parse-no-validate" | "parse_raw_no_validate_8m_slab" => Some(Self::ParseNoValidate),
            "pack" | "parse_and_pack_seq_qual" => Some(Self::Pack),
            "trusted-pack" | "trusted_default_pack_seq_qual" | "trusted_direct_pack_seq_qual" => {
                Some(Self::TrustedPack)
            }
            "fasta-count" | "fasta_count_resident_wrapped" => Some(Self::FastaCountResident),
            "fasta-visit" | "fasta_visit_stream_wrapped" => Some(Self::FastaVisitStream),
            "fasta-index" | "fasta_index_build_wrapped" => Some(Self::FastaIndexBuild),
            "fasta-fetch" | "fasta_fetch_indexed_wrapped" => Some(Self::FastaFetchIndexed),
            #[cfg(feature = "bgzf")]
            "bgzf" | "bgzf_parallel_stream_parse" | "bgzf_parallel_decode_then_parse" => {
                Some(Self::BgzfDecodeParse)
            }
            _ => None,
        }
    }

    fn should_run(self, requested: Self) -> bool {
        requested == Self::All || self == requested
    }
}

struct BenchConfig {
    case: Case,
    records: usize,
    read_len: usize,
    iters: usize,
}

impl BenchConfig {
    fn from_env_and_args() -> BenchResult<Self> {
        let mut config = Self {
            case: env::var("DINO_SEQ_BENCH_CASE")
                .ok()
                .as_deref()
                .and_then(Case::parse)
                .unwrap_or(Case::All),
            records: env_parse("DINO_SEQ_BENCH_RECORDS")?.unwrap_or(RECORDS),
            read_len: env_parse("DINO_SEQ_BENCH_READ_LEN")?.unwrap_or(READ_LEN),
            iters: env_parse("DINO_SEQ_BENCH_ITERS")?.unwrap_or(ITERS),
        };

        let mut args = env::args().skip(1);
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--case" => {
                    let value = next_arg_value(&mut args, "--case")?;
                    config.case = Case::parse(&value)
                        .ok_or_else(|| format!("unknown benchmark case `{value}`"))?;
                }
                "--records" => {
                    config.records = next_arg_value(&mut args, "--records")?.parse()?;
                }
                "--read-len" => {
                    config.read_len = next_arg_value(&mut args, "--read-len")?.parse()?;
                }
                "--iters" => {
                    config.iters = next_arg_value(&mut args, "--iters")?.parse()?;
                }
                "--list" => {
                    print_cases();
                    std::process::exit(0);
                }
                "--bench" => {}
                other => return Err(format!("unknown argument `{other}`").into()),
            }
        }

        Ok(config)
    }
}

fn env_parse(name: &str) -> BenchResult<Option<usize>> {
    match env::var(name) {
        Ok(value) => Ok(Some(value.parse()?)),
        Err(env::VarError::NotPresent) => Ok(None),
        Err(err) => Err(Box::new(err)),
    }
}

fn next_arg_value(args: &mut impl Iterator<Item = String>, flag: &str) -> BenchResult<String> {
    args.next()
        .ok_or_else(|| format!("{flag} requires a value").into())
}

fn print_cases() {
    println!("all");
    println!("parse");
    println!("visit");
    println!("visit-no-validate");
    println!("parse-no-validate");
    println!("pack");
    println!("trusted-pack");
    println!("fasta-count");
    println!("fasta-visit");
    println!("fasta-index");
    println!("fasta-fetch");
    #[cfg(feature = "bgzf")]
    println!("bgzf");
}

fn synthetic_fastq(records: usize, read_len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(records * (read_len * 2 + 32));
    for i in 0..records {
        out.extend_from_slice(b"@r");
        out.extend_from_slice(i.to_string().as_bytes());
        out.extend_from_slice(b"\n");
        for j in 0..read_len {
            out.push(b"ACGT"[(i + j) & 3]);
        }
        out.extend_from_slice(b"\n+\n");
        out.extend(std::iter::repeat_n(b'I', read_len));
        out.extend_from_slice(b"\n");
    }
    out
}

fn synthetic_fasta(records: usize, read_len: usize, wrap: usize) -> Vec<u8> {
    let wrap = wrap.max(1);
    let mut out = Vec::with_capacity(records * (read_len + read_len / wrap + 32));
    for i in 0..records {
        out.extend_from_slice(b">seq");
        out.extend_from_slice(i.to_string().as_bytes());
        out.extend_from_slice(b"\n");
        write_wrapped_sequence(&mut out, i, read_len, wrap);
    }
    out
}

fn synthetic_reference_fasta(name: &[u8], bases: usize, wrap: usize) -> Vec<u8> {
    let wrap = wrap.max(1);
    let mut out = Vec::with_capacity(bases + bases / wrap + name.len() + 8);
    out.extend_from_slice(b">");
    out.extend_from_slice(name);
    out.extend_from_slice(b"\n");
    write_wrapped_sequence(&mut out, 0, bases, wrap);
    out
}

fn write_wrapped_sequence(out: &mut Vec<u8>, seed: usize, len: usize, wrap: usize) {
    let mut written = 0;
    while written < len {
        let take = (len - written).min(wrap);
        for j in 0..take {
            out.push(b"ACGT"[(seed + written + j) & 3]);
        }
        out.push(b'\n');
        written += take;
    }
}

fn run_case<F>(name: &str, bytes: u64, iters: usize, mut f: F) -> BenchResult<()>
where
    F: FnMut() -> BenchResult<u64>,
{
    let start = Instant::now();
    let mut checksum = 0_u64;
    for _ in 0..iters {
        checksum = checksum.wrapping_add(black_box(f()?));
    }
    let elapsed = start.elapsed();
    let processed = bytes * iters as u64;
    let mib = processed as f64 / (1024.0 * 1024.0);
    let seconds = elapsed.as_secs_f64();
    println!(
        "{name}\t{mib:.2} MiB\t{seconds:.4} s\t{:.2} MiB/s",
        mib / seconds
    );
    black_box(checksum);
    Ok(())
}

fn parse_raw_strict_8m_slab(input: &[u8]) -> BenchResult<u64> {
    let mut reader = FastqReader::with_config(
        input,
        FastqConfig {
            slab_size: DEFAULT_SLAB_SIZE,
            validate: true,
            ..FastqConfig::default()
        },
    );
    Ok(consume_fastq(&mut reader)?.checksum)
}

fn visit_raw_strict_8m_slab(input: &[u8]) -> BenchResult<u64> {
    let mut reader = FastqReader::with_config(
        input,
        FastqConfig {
            slab_size: DEFAULT_SLAB_SIZE,
            validate: true,
            ..FastqConfig::default()
        },
    );
    let mut stats = StreamStats::default();
    reader.visit_records(|record: FastqVisitRecord<'_>| {
        stats.observe_record(record.name(), record.seq(), record.qual());
        Ok(())
    })?;
    Ok(stats.checksum)
}

fn visit_raw_no_validate_8m_slab(input: &[u8]) -> BenchResult<u64> {
    let mut reader = FastqReader::with_config(
        input,
        FastqConfig {
            slab_size: DEFAULT_SLAB_SIZE,
            validate: false,
            ..FastqConfig::default()
        },
    );
    let mut stats = StreamStats::default();
    reader.visit_records(|record: FastqVisitRecord<'_>| {
        stats.observe_record(record.name(), record.seq(), record.qual());
        Ok(())
    })?;
    Ok(stats.checksum)
}

fn parse_raw_no_validate_8m_slab(input: &[u8]) -> BenchResult<u64> {
    let mut reader = FastqReader::with_config(
        input,
        FastqConfig {
            slab_size: DEFAULT_SLAB_SIZE,
            validate: false,
            ..FastqConfig::default()
        },
    );
    Ok(consume_fastq(&mut reader)?.checksum)
}

fn parse_and_pack_seq_qual(input: &[u8]) -> BenchResult<u64> {
    let mut reader = FastqReader::new(input);
    let mut packed = Vec::new();
    let mut mask = Vec::new();
    let mut checksum = 0_u64;
    while let Some(batch) = reader.next_batch()? {
        for record in batch.records() {
            let summary = pack_bases_and_summarize_qualities_into(
                record.seq(),
                record.qual(),
                &mut packed,
                &mut mask,
            )?;
            checksum = checksum
                .wrapping_add(summary.bases.canonical_bases() as u64)
                .wrapping_add(summary.qualities.sum_phred);
        }
    }
    Ok(checksum)
}

fn trusted_default_pack_seq_qual(input: &[u8]) -> BenchResult<u64> {
    let mut checksum = 0_u64;
    pack_trusted_fastq(input, |record| {
        checksum = checksum
            .wrapping_add(record.summary.bases.canonical_bases() as u64)
            .wrapping_add(record.summary.qualities.sum_phred);
        Ok(())
    })?;
    Ok(checksum)
}

fn fasta_count_resident_wrapped(input: &[u8]) -> BenchResult<u64> {
    Ok(count_fasta_bytes(input)?.checksum)
}

fn fasta_visit_stream_wrapped(input: &[u8]) -> BenchResult<u64> {
    let mut reader = FastaReader::new(input);
    let mut checksum = 0_u64;
    reader.visit_records(|record| {
        checksum = checksum
            .wrapping_add(record.name().len() as u64)
            .wrapping_mul(1_099_511_628_211)
            .wrapping_add(record.seq().len() as u64);
        Ok(())
    })?;
    Ok(checksum)
}

fn fasta_index_build_wrapped(input: &[u8]) -> BenchResult<u64> {
    let index = build_fasta_index(input)?;
    Ok(index
        .entries()
        .iter()
        .fold(index.len() as u64, |acc, entry| acc.wrapping_add(entry.len)))
}

fn fasta_fetch_indexed_wrapped(input: &[u8], index: &FastaIndex, bases: usize) -> BenchResult<u64> {
    let mut reader = IndexedFastaReader::new(Cursor::new(input), index.clone());
    let fetched = reader.fetch(b"chr1", 0..bases as u64)?;
    Ok(fetched.first().copied().unwrap_or_default() as u64
        + fetched.last().copied().unwrap_or_default() as u64
        + fetched.len() as u64)
}

#[cfg(feature = "bgzf")]
fn bgzf_parallel_stream_parse(encoded: Arc<[u8]>) -> BenchResult<u64> {
    let bgzf = dino_seq::BgzfParallelReader::new(Cursor::new(encoded), 4)?;
    let mut reader = FastqReader::new(bgzf);
    Ok(consume_fastq(&mut reader)?.checksum)
}

fn main() -> BenchResult<()> {
    let config = BenchConfig::from_env_and_args()?;
    let input = synthetic_fastq(config.records, config.read_len);
    let bytes = input.len() as u64;
    let fasta_input = synthetic_fasta(config.records, config.read_len, 80);
    let fasta_bytes = fasta_input.len() as u64;
    let fetch_bases = config.records.saturating_mul(config.read_len);
    let fasta_reference = synthetic_reference_fasta(b"chr1", fetch_bases, 80);
    if Case::ParseStrict.should_run(config.case) {
        run_case(Case::ParseStrict.name(), bytes, config.iters, || {
            parse_raw_strict_8m_slab(&input)
        })?;
    }
    if Case::VisitStrict.should_run(config.case) {
        run_case(Case::VisitStrict.name(), bytes, config.iters, || {
            visit_raw_strict_8m_slab(&input)
        })?;
    }
    if Case::VisitNoValidate.should_run(config.case) {
        run_case(Case::VisitNoValidate.name(), bytes, config.iters, || {
            visit_raw_no_validate_8m_slab(&input)
        })?;
    }
    if Case::ParseNoValidate.should_run(config.case) {
        run_case(Case::ParseNoValidate.name(), bytes, config.iters, || {
            parse_raw_no_validate_8m_slab(&input)
        })?;
    }
    if Case::Pack.should_run(config.case) {
        run_case(Case::Pack.name(), bytes, config.iters, || {
            parse_and_pack_seq_qual(&input)
        })?;
    }
    if Case::TrustedPack.should_run(config.case) {
        run_case(Case::TrustedPack.name(), bytes, config.iters, || {
            trusted_default_pack_seq_qual(&input)
        })?;
    }
    if Case::FastaCountResident.should_run(config.case) {
        run_case(
            Case::FastaCountResident.name(),
            fasta_bytes,
            config.iters,
            || fasta_count_resident_wrapped(&fasta_input),
        )?;
    }
    if Case::FastaVisitStream.should_run(config.case) {
        run_case(
            Case::FastaVisitStream.name(),
            fasta_bytes,
            config.iters,
            || fasta_visit_stream_wrapped(&fasta_input),
        )?;
    }
    if Case::FastaIndexBuild.should_run(config.case) {
        run_case(
            Case::FastaIndexBuild.name(),
            fasta_bytes,
            config.iters,
            || fasta_index_build_wrapped(&fasta_input),
        )?;
    }
    if Case::FastaFetchIndexed.should_run(config.case) {
        let fasta_reference_index = build_fasta_index(fasta_reference.as_slice())?;
        run_case(
            Case::FastaFetchIndexed.name(),
            fetch_bases as u64,
            config.iters,
            || fasta_fetch_indexed_wrapped(&fasta_reference, &fasta_reference_index, fetch_bases),
        )?;
    }

    #[cfg(feature = "bgzf")]
    if Case::BgzfDecodeParse.should_run(config.case) {
        let encoded: Arc<[u8]> = dino_seq::compress_bgzf_parallel(&input, 4)?.into();
        run_case(Case::BgzfDecodeParse.name(), bytes, config.iters, || {
            bgzf_parallel_stream_parse(encoded.clone())
        })?;
    }

    Ok(())
}
