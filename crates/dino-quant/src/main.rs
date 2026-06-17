use dino_quant::{
    QuantizerConfig, ReferenceIndexConfig, ReferenceWindowIndex, SearchHit, concatenate_records,
    cosine_similarity, dna_kmer_sketch, dot, intervals_overlap, mutate_dna, read_fasta_file,
    reconstruction_metrics, synthetic_dna, visit_fastq_slices_file,
};
use std::time::{Duration, Instant};

fn main() {
    let mut args = std::env::args();
    let program = args.next().unwrap_or_else(|| "dino-quant".to_owned());
    match args.next().as_deref() {
        Some("demo") => {
            if let Err(err) = run_demo() {
                eprintln!("error: {err}");
                std::process::exit(1);
            }
        }
        Some("bench") => {
            if let Err(err) = run_bench() {
                eprintln!("error: {err}");
                std::process::exit(1);
            }
        }
        Some("bench-fasta") => {
            if let Err(err) = run_bench_fasta(args.collect()) {
                eprintln!("error: {err}");
                std::process::exit(1);
            }
        }
        Some("sweep-fasta") => {
            if let Err(err) = run_sweep_fasta(args.collect()) {
                eprintln!("error: {err}");
                std::process::exit(1);
            }
        }
        Some("emit-candidates") => {
            if let Err(err) = run_emit_candidates(args.collect()) {
                eprintln!("error: {err}");
                std::process::exit(1);
            }
        }
        Some("-h" | "--help") | None => print_help(&program),
        Some(other) => {
            eprintln!("unknown command: {other}");
            print_help(&program);
            std::process::exit(2);
        }
    }
}

fn print_help(program: &str) {
    println!("{program} demo");
    println!("{program} bench");
    println!(
        "{program} bench-fasta <reference.fa> [--queries N] [--top-k N] [--k N] [--dim N] [--window N] [--stride N] [--bits N] [--max-bases N] [--qjl]"
    );
    println!("{program} sweep-fasta <reference.fa> [bench-fasta options]");
    println!(
        "{program} emit-candidates <reference.fa> <reads.fastq> [--top-k N] [--k N] [--dim N] [--window N] [--stride N] [--bits N] [--max-bases N] [--qjl]"
    );
}

fn run_demo() -> Result<(), String> {
    let dim = 256;
    let k = 17;
    let reference_seq = synthetic_dna(20_000, 0xd1a0_5eed);
    let query_seq = mutate_dna(&reference_seq, 29);
    let reference = dna_kmer_sketch(&reference_seq, k, dim)?;
    let query = dna_kmer_sketch(&query_seq, k, dim)?;
    let raw_dot = dot(&reference, &query)?;
    let raw_cos = cosine_similarity(&reference, &query)?;

    println!("dino-quant demo: TurboQuant-inspired DNA sketch compression");
    println!("window_bases={} k={} dim={dim}", reference_seq.len(), k);
    println!("raw_query_cos={raw_cos:.6} raw_query_dot={raw_dot:.6}");
    println!();
    println!(
        "{:<5} {:<5} {:>8} {:>12} {:>10} {:>12} {:>15}",
        "bits", "qjl", "ratio", "mse", "cos", "dot_err", "query_dot_err"
    );

    for bits in [2_u8, 3, 4, 5] {
        for use_qjl_residual in [false, true] {
            let config = QuantizerConfig {
                bits,
                use_qjl_residual,
                ..QuantizerConfig::default()
            };
            let quantized = config.encode(&reference)?;
            let decoded = quantized.decode()?;
            let metrics = reconstruction_metrics(&reference, &decoded)?;
            let decoded_dot = dot(&decoded, &query)?;
            let query_dot_err = (raw_dot - decoded_dot).abs();
            let raw_bits = reference.len() * 32;
            let ratio = raw_bits as f32 / quantized.compressed_bits() as f32;

            println!(
                "{:<5} {:<5} {:>8.2} {:>12.6} {:>10.6} {:>12.6} {:>15.6}",
                bits,
                if use_qjl_residual { "yes" } else { "no" },
                ratio,
                metrics.mse,
                metrics.cosine,
                metrics.dot_error,
                query_dot_err
            );
        }
    }

    Ok(())
}

#[derive(Clone, Debug)]
struct ExactWindow {
    start: usize,
    end: usize,
    sketch: Vec<f32>,
}

#[derive(Clone, Debug)]
struct QueryCase {
    start: usize,
    end: usize,
    sketch: Vec<f32>,
}

#[derive(Clone, Copy, Debug)]
struct RecallCounts {
    top1: usize,
    top5: usize,
    top10: usize,
}

impl RecallCounts {
    fn rates(self, total: usize) -> (f32, f32, f32) {
        let total = total as f32;
        (
            self.top1 as f32 / total,
            self.top5 as f32 / total,
            self.top10 as f32 / total,
        )
    }
}

#[derive(Clone, Copy, Debug)]
struct BenchOptions {
    query_count: usize,
    top_k: usize,
    k: usize,
    dim: usize,
    window_len: usize,
    stride: usize,
    bits: u8,
    max_bases: Option<usize>,
    use_qjl_residual: bool,
}

#[derive(Clone, Debug)]
struct SweepOptions {
    bench: BenchOptions,
    dims: Vec<usize>,
    bits: Vec<u8>,
}

#[derive(Clone, Debug)]
struct EmitOptions {
    reference_path: String,
    reads_path: String,
    bench: BenchOptions,
}

#[derive(Clone, Copy, Debug)]
struct PrefilterMetrics {
    windows: usize,
    exact_bytes: usize,
    compressed_bytes: usize,
    exact_build: Duration,
    compressed_build: Duration,
    exact_search: Duration,
    compressed_search: Duration,
    exact_recall: RecallCounts,
    compressed_recall: RecallCounts,
}

impl PrefilterMetrics {
    fn compression_ratio(self) -> f32 {
        self.exact_bytes as f32 / self.compressed_bytes.max(1) as f32
    }

    fn candidate_reduction(self, top_k: usize) -> f32 {
        self.windows as f32 / top_k.max(1) as f32
    }

    fn exact_qps(self, query_count: usize) -> f32 {
        qps(query_count, self.exact_search)
    }

    fn compressed_qps(self, query_count: usize) -> f32 {
        qps(query_count, self.compressed_search)
    }
}

impl Default for BenchOptions {
    fn default() -> Self {
        Self {
            query_count: 256,
            top_k: 10,
            k: 17,
            dim: 256,
            window_len: 1024,
            stride: 256,
            bits: 4,
            max_bases: None,
            use_qjl_residual: false,
        }
    }
}

fn run_bench() -> Result<(), String> {
    let reference_len = 1_000_000;
    let reference = synthetic_dna(reference_len, 0x571a_9e2d);
    run_prefilter_bench("synthetic", &reference, BenchOptions::default())
}

fn run_bench_fasta(args: Vec<String>) -> Result<(), String> {
    let (path, options) = parse_bench_fasta_args(&args)?;
    let reference = load_reference(&path, options.max_bases)?;

    run_prefilter_bench(&format!("fasta:{path}"), &reference, options)
}

fn run_sweep_fasta(args: Vec<String>) -> Result<(), String> {
    let (path, sweep) = parse_sweep_fasta_args(&args)?;
    let reference = load_reference(&path, sweep.bench.max_bases)?;
    println!(
        "source\treference_bases\twindows\tqueries\ttop_k\tk\tdim\twindow\tstride\tbits\tqjl\tcompression_ratio\tcandidate_reduction\texact_qps\tcompressed_qps\tcompressed_vs_exact\trecall1\trecall5\trecall10"
    );

    for dim in &sweep.dims {
        for bits in &sweep.bits {
            let options = BenchOptions {
                dim: *dim,
                bits: *bits,
                ..sweep.bench
            };
            let config = build_index_config(options)?;
            let metrics =
                measure_prefilter(&reference, config, options.query_count, options.top_k)?;
            let compressed_qps = metrics.compressed_qps(options.query_count);
            let exact_qps = metrics.exact_qps(options.query_count);
            let (_, _, _) = metrics.exact_recall.rates(options.query_count);
            let (recall1, recall5, recall10) = metrics.compressed_recall.rates(options.query_count);
            println!(
                "fasta:{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.2}\t{:.2}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                path,
                reference.len(),
                metrics.windows,
                options.query_count,
                options.top_k,
                config.k,
                config.dim,
                config.window_len,
                config.stride,
                config.quantizer.bits,
                config.quantizer.use_qjl_residual,
                metrics.compression_ratio(),
                metrics.candidate_reduction(options.top_k),
                exact_qps,
                compressed_qps,
                compressed_qps / exact_qps.max(f32::EPSILON),
                recall1,
                recall5,
                recall10
            );
        }
    }

    Ok(())
}

fn run_emit_candidates(args: Vec<String>) -> Result<(), String> {
    let options = parse_emit_args(&args)?;
    let reference = load_reference(&options.reference_path, options.bench.max_bases)?;
    let config = build_index_config(options.bench)?;
    let build_start = Instant::now();
    let index = ReferenceWindowIndex::build(&reference, config)?;
    let build_elapsed = build_start.elapsed();

    println!("read_name\tread_len\trank\tref_linear_start\tref_linear_end\tscore");
    let mut reads = 0_usize;
    let mut candidates = 0_usize;
    let search_start = Instant::now();
    visit_fastq_slices_file(&options.reads_path, |record| {
        let hits = index.search_sequence(record.bases, options.bench.top_k)?;
        reads += 1;
        let read_name = String::from_utf8_lossy(record.name);
        for (rank, hit) in hits.iter().enumerate() {
            candidates += 1;
            println!(
                "{}\t{}\t{}\t{}\t{}\t{:.6}",
                read_name,
                record.bases.len(),
                rank + 1,
                hit.start,
                hit.end,
                hit.score
            );
        }
        Ok(())
    })?;
    let search_elapsed = search_start.elapsed();

    eprintln!(
        "dino-quant emit-candidates: reference_bases={} windows={} reads={} candidates={} build_ms={:.2} read_qps={:.2} compression_ratio={:.2}",
        reference.len(),
        index.window_count(),
        reads,
        candidates,
        ms(build_elapsed),
        qps(reads, search_elapsed),
        index.compression_ratio()
    );
    Ok(())
}

fn run_prefilter_bench(label: &str, reference: &[u8], options: BenchOptions) -> Result<(), String> {
    if options.query_count == 0 {
        return Err("query count must be non-zero".to_owned());
    }
    if options.top_k == 0 {
        return Err("top-k must be non-zero".to_owned());
    }

    let config = ReferenceIndexConfig {
        k: options.k,
        dim: options.dim,
        window_len: options.window_len,
        stride: options.stride,
        quantizer: QuantizerConfig {
            bits: options.bits,
            use_qjl_residual: options.use_qjl_residual,
            ..QuantizerConfig::default()
        },
    };
    config.validate()?;
    if reference.len() < config.window_len {
        return Err(format!(
            "reference has {} bases but window length is {}",
            reference.len(),
            config.window_len
        ));
    }

    let query_count = options.query_count;
    let top_k = options.top_k;
    let metrics = measure_prefilter(reference, config, query_count, top_k)?;
    let exact_qps = metrics.exact_qps(query_count);
    let compressed_qps = metrics.compressed_qps(query_count);
    let (exact_r1, exact_r5, exact_r10) = metrics.exact_recall.rates(query_count);
    let (comp_r1, comp_r5, comp_r10) = metrics.compressed_recall.rates(query_count);

    println!("dino-quant bench: compressed DNA window prefilter");
    println!(
        "source={label} reference_bases={} windows={} queries={query_count} k={} dim={} window={} stride={} bits={} qjl={}",
        reference.len(),
        metrics.windows,
        config.k,
        config.dim,
        config.window_len,
        config.stride,
        config.quantizer.bits,
        config.quantizer.use_qjl_residual
    );
    println!();
    println!(
        "index_bytes_exact={} index_bytes_compressed={} compression_ratio={:.2}",
        metrics.exact_bytes,
        metrics.compressed_bytes,
        metrics.compression_ratio()
    );
    println!(
        "candidate_reduction_at_top{top_k}={:.2}x",
        metrics.candidate_reduction(top_k)
    );
    println!();
    println!(
        "{:<11} {:>12} {:>12} {:>10} {:>10} {:>10}",
        "mode", "build_ms", "query_qps", "recall@1", "recall@5", "recall@10"
    );
    println!(
        "{:<11} {:>12.2} {:>12.2} {:>10.3} {:>10.3} {:>10.3}",
        "exact",
        ms(metrics.exact_build),
        exact_qps,
        exact_r1,
        exact_r5,
        exact_r10
    );
    println!(
        "{:<11} {:>12.2} {:>12.2} {:>10.3} {:>10.3} {:>10.3}",
        "compressed",
        ms(metrics.compressed_build),
        compressed_qps,
        comp_r1,
        comp_r5,
        comp_r10
    );
    println!(
        "compressed_qps_vs_exact={:.2}x",
        compressed_qps / exact_qps.max(f32::EPSILON)
    );

    Ok(())
}

fn parse_bench_fasta_args(args: &[String]) -> Result<(String, BenchOptions), String> {
    let path = args
        .first()
        .ok_or_else(|| "bench-fasta requires a reference FASTA path".to_owned())?
        .clone();
    let mut options = BenchOptions::default();
    let mut idx = 1;
    while idx < args.len() {
        match args[idx].as_str() {
            "--queries" => {
                options.query_count = parse_next_usize(args, &mut idx, "--queries")?;
            }
            "--top-k" => {
                options.top_k = parse_next_usize(args, &mut idx, "--top-k")?;
            }
            "--k" => {
                options.k = parse_next_usize(args, &mut idx, "--k")?;
            }
            "--dim" => {
                options.dim = parse_next_usize(args, &mut idx, "--dim")?;
            }
            "--window" => {
                options.window_len = parse_next_usize(args, &mut idx, "--window")?;
            }
            "--stride" => {
                options.stride = parse_next_usize(args, &mut idx, "--stride")?;
            }
            "--bits" => {
                let bits = parse_next_usize(args, &mut idx, "--bits")?;
                options.bits = u8::try_from(bits)
                    .map_err(|_| "--bits must fit in an unsigned 8-bit integer".to_owned())?;
            }
            "--max-bases" => {
                options.max_bases = Some(parse_next_usize(args, &mut idx, "--max-bases")?);
            }
            "--qjl" => {
                options.use_qjl_residual = true;
                idx += 1;
            }
            other => return Err(format!("unknown bench-fasta option: {other}")),
        }
    }
    Ok((path, options))
}

fn parse_sweep_fasta_args(args: &[String]) -> Result<(String, SweepOptions), String> {
    let path = args
        .first()
        .ok_or_else(|| "sweep-fasta requires a reference FASTA path".to_owned())?
        .clone();
    let mut bench = BenchOptions::default();
    let mut dims = vec![128, 256, 512];
    let mut bits = vec![2, 3, 4, 5];
    let mut idx = 1;
    while idx < args.len() {
        match args[idx].as_str() {
            "--dims" => {
                dims = parse_next_usize_list(args, &mut idx, "--dims")?;
            }
            "--bits-list" => {
                bits = parse_next_u8_list(args, &mut idx, "--bits-list")?;
            }
            "--queries" => {
                bench.query_count = parse_next_usize(args, &mut idx, "--queries")?;
            }
            "--top-k" => {
                bench.top_k = parse_next_usize(args, &mut idx, "--top-k")?;
            }
            "--k" => {
                bench.k = parse_next_usize(args, &mut idx, "--k")?;
            }
            "--window" => {
                bench.window_len = parse_next_usize(args, &mut idx, "--window")?;
            }
            "--stride" => {
                bench.stride = parse_next_usize(args, &mut idx, "--stride")?;
            }
            "--max-bases" => {
                bench.max_bases = Some(parse_next_usize(args, &mut idx, "--max-bases")?);
            }
            "--qjl" => {
                bench.use_qjl_residual = true;
                idx += 1;
            }
            other => return Err(format!("unknown sweep-fasta option: {other}")),
        }
    }
    Ok((path, SweepOptions { bench, dims, bits }))
}

fn parse_emit_args(args: &[String]) -> Result<EmitOptions, String> {
    let reference_path = args
        .first()
        .ok_or_else(|| "emit-candidates requires a reference FASTA path".to_owned())?
        .clone();
    let reads_path = args
        .get(1)
        .ok_or_else(|| "emit-candidates requires a reads FASTQ path".to_owned())?
        .clone();
    let mut bench = BenchOptions::default();
    let mut idx = 2;
    while idx < args.len() {
        match args[idx].as_str() {
            "--top-k" => {
                bench.top_k = parse_next_usize(args, &mut idx, "--top-k")?;
            }
            "--k" => {
                bench.k = parse_next_usize(args, &mut idx, "--k")?;
            }
            "--dim" => {
                bench.dim = parse_next_usize(args, &mut idx, "--dim")?;
            }
            "--window" => {
                bench.window_len = parse_next_usize(args, &mut idx, "--window")?;
            }
            "--stride" => {
                bench.stride = parse_next_usize(args, &mut idx, "--stride")?;
            }
            "--bits" => {
                let bits = parse_next_usize(args, &mut idx, "--bits")?;
                bench.bits = u8::try_from(bits)
                    .map_err(|_| "--bits must fit in an unsigned 8-bit integer".to_owned())?;
            }
            "--max-bases" => {
                bench.max_bases = Some(parse_next_usize(args, &mut idx, "--max-bases")?);
            }
            "--qjl" => {
                bench.use_qjl_residual = true;
                idx += 1;
            }
            other => return Err(format!("unknown emit-candidates option: {other}")),
        }
    }
    Ok(EmitOptions {
        reference_path,
        reads_path,
        bench,
    })
}

fn parse_next_usize(args: &[String], idx: &mut usize, flag: &str) -> Result<usize, String> {
    let value_idx = *idx + 1;
    let value = args
        .get(value_idx)
        .ok_or_else(|| format!("{flag} requires a value"))?;
    let parsed = value
        .parse::<usize>()
        .map_err(|err| format!("{flag} value must be a positive integer: {err}"))?;
    *idx += 2;
    Ok(parsed)
}

fn parse_next_usize_list(
    args: &[String],
    idx: &mut usize,
    flag: &str,
) -> Result<Vec<usize>, String> {
    let value_idx = *idx + 1;
    let value = args
        .get(value_idx)
        .ok_or_else(|| format!("{flag} requires a comma-separated value"))?;
    let mut values = Vec::new();
    for part in value.split(',') {
        values.push(
            part.parse::<usize>()
                .map_err(|err| format!("{flag} values must be positive integers: {err}"))?,
        );
    }
    if values.is_empty() {
        return Err(format!("{flag} requires at least one value"));
    }
    *idx += 2;
    Ok(values)
}

fn parse_next_u8_list(args: &[String], idx: &mut usize, flag: &str) -> Result<Vec<u8>, String> {
    parse_next_usize_list(args, idx, flag)?
        .into_iter()
        .map(|value| u8::try_from(value).map_err(|_| format!("{flag} values must fit in u8")))
        .collect()
}

fn load_reference(path: &str, max_bases: Option<usize>) -> Result<Vec<u8>, String> {
    let records = read_fasta_file(path)?;
    let mut reference = concatenate_records(&records)?;
    if let Some(max_bases) = max_bases {
        reference.truncate(max_bases.min(reference.len()));
    }
    Ok(reference)
}

fn build_index_config(options: BenchOptions) -> Result<ReferenceIndexConfig, String> {
    if options.top_k == 0 {
        return Err("top-k must be non-zero".to_owned());
    }
    if options.query_count == 0 {
        return Err("query count must be non-zero".to_owned());
    }
    let config = ReferenceIndexConfig {
        k: options.k,
        dim: options.dim,
        window_len: options.window_len,
        stride: options.stride,
        quantizer: QuantizerConfig {
            bits: options.bits,
            use_qjl_residual: options.use_qjl_residual,
            ..QuantizerConfig::default()
        },
    };
    config.validate()?;
    Ok(config)
}

fn measure_prefilter(
    reference: &[u8],
    config: ReferenceIndexConfig,
    query_count: usize,
    top_k: usize,
) -> Result<PrefilterMetrics, String> {
    let exact_build_start = Instant::now();
    let exact_windows = build_exact_windows(reference, config)?;
    let exact_build = exact_build_start.elapsed();

    let compressed_build_start = Instant::now();
    let compressed_index = ReferenceWindowIndex::build(reference, config)?;
    let compressed_build = compressed_build_start.elapsed();

    let queries = build_queries(reference, config, query_count)?;

    let exact_start = Instant::now();
    let mut exact_recall = RecallCounts {
        top1: 0,
        top5: 0,
        top10: 0,
    };
    for query in &queries {
        let hits = search_exact(&exact_windows, &query.sketch, top_k)?;
        count_recall(&mut exact_recall, &hits, query);
    }
    let exact_search = exact_start.elapsed();

    let compressed_start = Instant::now();
    let mut compressed_recall = RecallCounts {
        top1: 0,
        top5: 0,
        top10: 0,
    };
    for query in &queries {
        let hits = compressed_index.search_sketch(&query.sketch, top_k)?;
        count_recall(&mut compressed_recall, &hits, query);
    }
    let compressed_search = compressed_start.elapsed();

    Ok(PrefilterMetrics {
        windows: compressed_index.window_count(),
        exact_bytes: exact_windows.len() * config.dim * std::mem::size_of::<f32>(),
        compressed_bytes: compressed_index.compressed_bytes(),
        exact_build,
        compressed_build,
        exact_search,
        compressed_search,
        exact_recall,
        compressed_recall,
    })
}

fn build_exact_windows(
    reference: &[u8],
    config: ReferenceIndexConfig,
) -> Result<Vec<ExactWindow>, String> {
    config.validate()?;
    if reference.len() < config.window_len {
        return Err("reference is shorter than the configured window length".to_owned());
    }

    let window_count = ((reference.len() - config.window_len) / config.stride) + 1;
    let mut windows = Vec::with_capacity(window_count);
    for start in (0..=reference.len() - config.window_len).step_by(config.stride) {
        let end = start + config.window_len;
        let sketch = dna_kmer_sketch(&reference[start..end], config.k, config.dim)?;
        windows.push(ExactWindow { start, end, sketch });
    }
    Ok(windows)
}

fn build_queries(
    reference: &[u8],
    config: ReferenceIndexConfig,
    count: usize,
) -> Result<Vec<QueryCase>, String> {
    let window_count = ((reference.len() - config.window_len) / config.stride) + 1;
    let step = (window_count / count.max(1)).max(1);
    let mut queries = Vec::with_capacity(count);
    for idx in 0..count {
        let window_idx = (idx * step + idx * idx * 17) % window_count;
        let start = window_idx * config.stride;
        let end = start + config.window_len;
        let mutated = mutate_dna(&reference[start..end], 43);
        let sketch = dna_kmer_sketch(&mutated, config.k, config.dim)?;
        queries.push(QueryCase { start, end, sketch });
    }
    Ok(queries)
}

fn search_exact(
    windows: &[ExactWindow],
    query_sketch: &[f32],
    top_k: usize,
) -> Result<Vec<SearchHit>, String> {
    if top_k == 0 {
        return Ok(Vec::new());
    }

    let mut top = Vec::with_capacity(top_k.min(windows.len()));
    for window in windows {
        push_top_hit(
            &mut top,
            top_k,
            SearchHit {
                start: window.start,
                end: window.end,
                score: dot(&window.sketch, query_sketch)?,
            },
        );
    }
    top.sort_by(|left, right| right.score.total_cmp(&left.score));
    Ok(top)
}

fn push_top_hit(top: &mut Vec<SearchHit>, top_k: usize, hit: SearchHit) {
    if top.len() < top_k {
        top.push(hit);
        return;
    }

    let mut worst_idx = 0;
    let mut worst_score = top[0].score;
    for (idx, current) in top.iter().enumerate().skip(1) {
        if current.score < worst_score {
            worst_idx = idx;
            worst_score = current.score;
        }
    }

    if hit.score > worst_score {
        top[worst_idx] = hit;
    }
}

fn count_recall(counts: &mut RecallCounts, hits: &[SearchHit], query: &QueryCase) {
    if is_hit(hits.first(), query) {
        counts.top1 += 1;
    }
    if hits
        .iter()
        .take(5)
        .any(|hit| intervals_overlap(hit.start, hit.end, query.start, query.end))
    {
        counts.top5 += 1;
    }
    if hits
        .iter()
        .take(10)
        .any(|hit| intervals_overlap(hit.start, hit.end, query.start, query.end))
    {
        counts.top10 += 1;
    }
}

fn is_hit(hit: Option<&SearchHit>, query: &QueryCase) -> bool {
    hit.is_some_and(|hit| intervals_overlap(hit.start, hit.end, query.start, query.end))
}

fn qps(query_count: usize, elapsed: Duration) -> f32 {
    query_count as f32 / elapsed.as_secs_f32().max(f32::EPSILON)
}

fn ms(elapsed: Duration) -> f32 {
    elapsed.as_secs_f32() * 1000.0
}
