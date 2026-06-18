use dino_quant::{
    QjlResidualSnapshot, QuantizedVector, QuantizedVectorSnapshot, QuantizerConfig,
    ReferenceIndexConfig, ReferenceWindowIndex, SearchHit, concatenate_records, cosine_similarity,
    dna_kmer_sketch, dot, intervals_overlap, mutate_dna, protein_kmer_sketch, read_fasta_file,
    read_protein_fasta_file, reconstruction_metrics, synthetic_dna, visit_fastq_slices_file_limit,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
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
        Some("emit-candidate-reference") => {
            if let Err(err) = run_emit_candidate_reference(args.collect()) {
                eprintln!("error: {err}");
                std::process::exit(1);
            }
        }
        Some("emit-protein-candidates") => {
            if let Err(err) = run_emit_protein_candidates(args.collect()) {
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
        "{program} emit-candidates <reference.fa> <reads.fastq> [--top-k N] [--retrieval scan|minimizer|simhash|mih|ivf|hnsw] [--candidate-limit N] [--reference-cache PATH] [--k N] [--dim N] [--window N] [--stride N] [--bits N] [--max-bases N] [--max-reads N] [--qjl]"
    );
    println!(
        "{program} emit-candidate-reference <reference.fa> <candidates.tsv> [--padding N] [--merge-gap N] [--mask-reference]"
    );
    println!(
        "{program} emit-protein-candidates <proteins.faa> <queries.faa> [--top-k N] [--k N] [--dim N] [--bits N] [--max-bases N] [--qjl]"
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
    max_reads: Option<usize>,
    bench: BenchOptions,
    retrieval: RetrievalOptions,
    reference_cache: Option<String>,
}

#[derive(Clone, Debug)]
struct CandidateReferenceOptions {
    reference_path: String,
    candidates_path: String,
    padding: usize,
    merge_gap: usize,
    mask_reference: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Interval {
    start: usize,
    end: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RetrievalMode {
    Scan,
    Minimizer,
    Simhash,
    Mih,
    Ivf,
    Hnsw,
}

impl RetrievalMode {
    fn parse(value: &str) -> Result<Self, String> {
        match value {
            "scan" => Ok(Self::Scan),
            "minimizer" => Ok(Self::Minimizer),
            "simhash" => Ok(Self::Simhash),
            "mih" => Ok(Self::Mih),
            "ivf" => Ok(Self::Ivf),
            "hnsw" => Ok(Self::Hnsw),
            other => Err(format!(
                "unknown retrieval mode: {other}; expected scan|minimizer|simhash|mih|ivf|hnsw"
            )),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Scan => "scan",
            Self::Minimizer => "minimizer",
            Self::Simhash => "simhash",
            Self::Mih => "mih",
            Self::Ivf => "ivf",
            Self::Hnsw => "hnsw",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct RetrievalOptions {
    mode: RetrievalMode,
    candidate_limit: usize,
    minimizer_k: usize,
    minimizer_window: usize,
    simhash_bands: usize,
    centroids: usize,
    probes: usize,
}

impl Default for RetrievalOptions {
    fn default() -> Self {
        Self {
            mode: RetrievalMode::Scan,
            candidate_limit: 2048,
            minimizer_k: 15,
            minimizer_window: 8,
            simhash_bands: 8,
            centroids: 32,
            probes: 4,
        }
    }
}

#[derive(Clone, Debug)]
struct IndexedSearchResult {
    hits: Vec<SearchHit>,
    considered: usize,
}

#[derive(Clone, Debug)]
struct CandidateScratch {
    counts: Vec<u32>,
    touched: Vec<usize>,
    ranked: Vec<(usize, u32)>,
    candidate_ids: Vec<usize>,
    minimizer_hashes: Vec<u64>,
    minimizers: Vec<u64>,
    minimizer_deque: Vec<usize>,
}

impl CandidateScratch {
    fn new(window_count: usize) -> Self {
        Self {
            counts: vec![0; window_count],
            touched: Vec::new(),
            ranked: Vec::new(),
            candidate_ids: Vec::new(),
            minimizer_hashes: Vec::new(),
            minimizers: Vec::new(),
            minimizer_deque: Vec::new(),
        }
    }

    fn reset_counts(&mut self) {
        for idx in self.touched.drain(..) {
            if idx < self.counts.len() {
                self.counts[idx] = 0;
            }
        }
        self.ranked.clear();
        self.candidate_ids.clear();
    }

    fn resize_counts(&mut self, window_count: usize) {
        self.reset_counts();
        self.counts.resize(window_count, 0);
    }

    fn bump(&mut self, idx: usize, weight: u32) -> Result<(), String> {
        let Some(count) = self.counts.get_mut(idx) else {
            return Err("candidate posting id out of bounds".to_owned());
        };
        if *count == 0 {
            self.touched.push(idx);
        }
        *count = count.saturating_add(weight);
        Ok(())
    }

    fn ranked_candidate_ids(&mut self, limit: usize) -> &[usize] {
        self.ranked.clear();
        self.ranked.reserve(self.touched.len());
        for idx in self.touched.drain(..) {
            let count = self.counts[idx];
            self.counts[idx] = 0;
            self.ranked.push((idx, count));
        }
        self.ranked
            .sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
        self.candidate_ids.clear();
        self.candidate_ids
            .extend(self.ranked.iter().take(limit).map(|(idx, _)| *idx));
        &self.candidate_ids
    }
}

#[derive(Clone, Debug)]
struct CandidateWindow {
    target_id: usize,
    target_start: usize,
    sketch: Vec<f32>,
    quantized: QuantizedVector,
    simhash64: u64,
    simhash128: [u64; 2],
}

#[derive(Clone, Debug)]
struct ExperimentalCandidateIndex {
    config: ReferenceIndexConfig,
    retrieval: RetrievalOptions,
    target_names: Vec<String>,
    target_offsets: Vec<usize>,
    windows: Vec<CandidateWindow>,
    minimizer_postings: HashMap<u64, Vec<usize>>,
    simhash_band_postings: Vec<HashMap<u16, Vec<usize>>>,
    mih_postings: Vec<HashMap<u16, Vec<usize>>>,
    centroids: Vec<Vec<f32>>,
    centroid_lists: Vec<Vec<usize>>,
    centroid_graph: Vec<Vec<usize>>,
}

impl ExperimentalCandidateIndex {
    fn build_records(
        records: &[dino_quant::SequenceRecord],
        config: ReferenceIndexConfig,
        retrieval: RetrievalOptions,
    ) -> Result<Self, String> {
        config.validate()?;
        validate_retrieval_options(retrieval)?;
        if records.is_empty() {
            return Err("at least one sequence record is required".to_owned());
        }

        let mut windows = Vec::new();
        let mut target_names = Vec::with_capacity(records.len());
        let mut target_offsets = Vec::with_capacity(records.len());
        let mut linear_offset = 0_usize;
        for (record_idx, record) in records.iter().enumerate() {
            target_names.push(record.name.clone());
            target_offsets.push(linear_offset);
            if record.bases.len() >= config.window_len {
                let window_count = ((record.bases.len() - config.window_len) / config.stride) + 1;
                windows.reserve(window_count);
                for target_start in
                    (0..=record.bases.len() - config.window_len).step_by(config.stride)
                {
                    let target_end = target_start + config.window_len;
                    let sketch = dna_kmer_sketch(
                        &record.bases[target_start..target_end],
                        config.k,
                        config.dim,
                    )?;
                    let quantized = config.quantizer.encode(&sketch)?;
                    let simhash64 = simhash_code64(&sketch, 0x51d0_a11e_0000_0064);
                    let simhash128 = [simhash64, simhash_code64(&sketch, 0x51d0_a11e_0000_0128)];
                    windows.push(CandidateWindow {
                        target_id: record_idx,
                        target_start,
                        sketch,
                        quantized,
                        simhash64,
                        simhash128,
                    });
                }
            }
            linear_offset += record.bases.len();
            if record_idx + 1 < records.len() {
                linear_offset += 1;
            }
        }

        if windows.is_empty() {
            return Err("no reference record is long enough for the configured window".to_owned());
        }

        let mut index = Self {
            config,
            retrieval,
            target_names,
            target_offsets,
            windows,
            minimizer_postings: HashMap::new(),
            simhash_band_postings: Vec::new(),
            mih_postings: Vec::new(),
            centroids: Vec::new(),
            centroid_lists: Vec::new(),
            centroid_graph: Vec::new(),
        };
        index.build_auxiliary(records)?;
        Ok(index)
    }

    fn build_auxiliary(&mut self, records: &[dino_quant::SequenceRecord]) -> Result<(), String> {
        match self.retrieval.mode {
            RetrievalMode::Scan => {}
            RetrievalMode::Minimizer => self.build_minimizer_postings(records)?,
            RetrievalMode::Simhash => self.build_simhash_band_postings(),
            RetrievalMode::Mih => self.build_mih_postings(),
            RetrievalMode::Ivf => self.build_centroid_lists(),
            RetrievalMode::Hnsw => {
                self.build_centroid_lists();
                self.build_centroid_graph();
            }
        }
        Ok(())
    }

    fn build_minimizer_postings(
        &mut self,
        records: &[dino_quant::SequenceRecord],
    ) -> Result<(), String> {
        let mut window_idx = 0_usize;
        let mut hashes = Vec::new();
        let mut minimizers = Vec::new();
        let mut deque = Vec::new();
        for record in records {
            if record.bases.len() >= self.config.window_len {
                for target_start in
                    (0..=record.bases.len() - self.config.window_len).step_by(self.config.stride)
                {
                    let target_end = target_start + self.config.window_len;
                    minimizer_hashes_into(
                        &record.bases[target_start..target_end],
                        self.retrieval.minimizer_k,
                        self.retrieval.minimizer_window,
                        &mut hashes,
                        &mut minimizers,
                        &mut deque,
                    )?;
                    for hash in &minimizers {
                        self.minimizer_postings
                            .entry(*hash)
                            .or_default()
                            .push(window_idx);
                    }
                    window_idx += 1;
                }
            }
        }
        Ok(())
    }

    fn build_simhash_band_postings(&mut self) {
        let bands = self.retrieval.simhash_bands;
        self.simhash_band_postings = (0..bands).map(|_| HashMap::new()).collect();
        for (idx, window) in self.windows.iter().enumerate() {
            for band in 0..bands {
                let key = simhash_band_key(window.simhash64, band, bands);
                self.simhash_band_postings[band]
                    .entry(key)
                    .or_default()
                    .push(idx);
            }
        }
    }

    fn build_mih_postings(&mut self) {
        self.mih_postings = (0..8).map(|_| HashMap::new()).collect();
        for (idx, window) in self.windows.iter().enumerate() {
            for part in 0..8 {
                let key = simhash128_part(window.simhash128, part);
                self.mih_postings[part].entry(key).or_default().push(idx);
            }
        }
    }

    fn build_centroid_lists(&mut self) {
        let centroid_count = self.retrieval.centroids.min(self.windows.len()).max(1);
        self.centroids.clear();
        for centroid_idx in 0..centroid_count {
            let window_idx = if centroid_count == 1 {
                0
            } else {
                centroid_idx * (self.windows.len() - 1) / (centroid_count - 1)
            };
            self.centroids.push(self.windows[window_idx].sketch.clone());
        }
        self.centroid_lists = (0..self.centroids.len()).map(|_| Vec::new()).collect();
        for (idx, window) in self.windows.iter().enumerate() {
            let centroid_idx = nearest_centroid(&window.sketch, &self.centroids);
            self.centroid_lists[centroid_idx].push(idx);
        }
    }

    fn build_centroid_graph(&mut self) {
        let neighbors = 8_usize.min(self.centroids.len().saturating_sub(1));
        self.centroid_graph = Vec::with_capacity(self.centroids.len());
        for idx in 0..self.centroids.len() {
            let mut scores = Vec::with_capacity(self.centroids.len().saturating_sub(1));
            for other in 0..self.centroids.len() {
                if idx == other {
                    continue;
                }
                scores.push((
                    other,
                    dot_unchecked(&self.centroids[idx], &self.centroids[other]),
                ));
            }
            scores.sort_by(|left, right| right.1.total_cmp(&left.1));
            self.centroid_graph.push(
                scores
                    .into_iter()
                    .take(neighbors)
                    .map(|(other, _)| other)
                    .collect(),
            );
        }
    }

    #[cfg(test)]
    fn search_sequence(&self, query: &[u8], top_k: usize) -> Result<IndexedSearchResult, String> {
        let mut scratch = CandidateScratch::new(self.windows.len());
        self.search_sequence_with_scratch(query, top_k, &mut scratch)
    }

    fn search_sequence_with_scratch(
        &self,
        query: &[u8],
        top_k: usize,
        scratch: &mut CandidateScratch,
    ) -> Result<IndexedSearchResult, String> {
        let query_sketch = dna_kmer_sketch(query, self.config.k, self.config.dim)?;
        let prepared_query = QuantizedVector::prepare_approximate_query(
            &query_sketch,
            self.config.quantizer.rotation_seed,
            self.config
                .quantizer
                .use_qjl_residual
                .then_some(self.config.quantizer.qjl_seed),
        )?;
        let ids = self.retrieve_candidate_ids(query, &query_sketch, scratch)?;
        let considered = ids.len();
        let mut top = Vec::with_capacity(top_k.min(considered));
        for &idx in ids {
            let window = &self.windows[idx];
            let score = window
                .quantized
                .approximate_dot_prepared_query(&prepared_query)?;
            let target_end = window.target_start + self.config.window_len;
            let linear_start = self.target_offsets[window.target_id] + window.target_start;
            push_top_hit(
                &mut top,
                top_k,
                SearchHit {
                    target_name: self.target_names[window.target_id].clone(),
                    target_start: window.target_start,
                    target_end,
                    start: linear_start,
                    end: linear_start + self.config.window_len,
                    score,
                },
            );
        }
        top.sort_by(|left, right| right.score.total_cmp(&left.score));
        Ok(IndexedSearchResult {
            hits: top,
            considered,
        })
    }

    fn retrieve_candidate_ids<'a>(
        &self,
        query: &[u8],
        sketch: &[f32],
        scratch: &'a mut CandidateScratch,
    ) -> Result<&'a [usize], String> {
        match self.retrieval.mode {
            RetrievalMode::Scan => {
                scratch.candidate_ids.clear();
                scratch.candidate_ids.extend(0..self.windows.len());
                Ok(&scratch.candidate_ids)
            }
            RetrievalMode::Minimizer => self.retrieve_minimizer(query, scratch),
            RetrievalMode::Simhash => {
                scratch.candidate_ids = self.retrieve_simhash(sketch);
                Ok(&scratch.candidate_ids)
            }
            RetrievalMode::Mih => {
                scratch.candidate_ids = self.retrieve_mih(sketch);
                Ok(&scratch.candidate_ids)
            }
            RetrievalMode::Ivf => {
                scratch.candidate_ids = self.retrieve_ivf(sketch);
                Ok(&scratch.candidate_ids)
            }
            RetrievalMode::Hnsw => {
                scratch.candidate_ids = self.retrieve_hnsw(sketch);
                Ok(&scratch.candidate_ids)
            }
        }
    }

    fn retrieve_minimizer<'a>(
        &self,
        query: &[u8],
        scratch: &'a mut CandidateScratch,
    ) -> Result<&'a [usize], String> {
        if scratch.counts.len() != self.windows.len() {
            scratch.resize_counts(self.windows.len());
        } else {
            scratch.reset_counts();
        }
        minimizer_hashes_into(
            query,
            self.retrieval.minimizer_k,
            self.retrieval.minimizer_window,
            &mut scratch.minimizer_hashes,
            &mut scratch.minimizers,
            &mut scratch.minimizer_deque,
        )?;
        for pos in 0..scratch.minimizers.len() {
            let hash = scratch.minimizers[pos];
            if let Some(posting) = self.minimizer_postings.get(&hash) {
                for &idx in posting {
                    scratch.bump(idx, 1)?;
                }
            }
        }
        Ok(scratch.ranked_candidate_ids(self.retrieval.candidate_limit))
    }

    fn retrieve_simhash(&self, sketch: &[f32]) -> Vec<usize> {
        let code = simhash_code64(sketch, 0x51d0_a11e_0000_0064);
        let mut counts = HashMap::new();
        let bands = self.retrieval.simhash_bands;
        for band in 0..bands {
            let key = simhash_band_key(code, band, bands);
            if let Some(posting) = self.simhash_band_postings[band].get(&key) {
                for &idx in posting {
                    bump_count(&mut counts, idx, 1);
                }
            }
        }
        top_counted_ids(counts, self.retrieval.candidate_limit)
    }

    fn retrieve_mih(&self, sketch: &[f32]) -> Vec<usize> {
        let code = [
            simhash_code64(sketch, 0x51d0_a11e_0000_0064),
            simhash_code64(sketch, 0x51d0_a11e_0000_0128),
        ];
        let mut counts = HashMap::new();
        for part in 0..8 {
            let key = simhash128_part(code, part);
            self.bump_mih_part(part, key, 4, &mut counts);
            for bit in 0..16 {
                self.bump_mih_part(part, key ^ (1_u16 << bit), 1, &mut counts);
            }
        }
        top_counted_ids(counts, self.retrieval.candidate_limit)
    }

    fn bump_mih_part(&self, part: usize, key: u16, weight: u32, counts: &mut HashMap<usize, u32>) {
        if let Some(posting) = self.mih_postings[part].get(&key) {
            for &idx in posting {
                bump_count(counts, idx, weight);
            }
        }
    }

    fn retrieve_ivf(&self, sketch: &[f32]) -> Vec<usize> {
        let centroids = top_centroids(sketch, &self.centroids, self.retrieval.probes);
        let mut ids = Vec::new();
        for centroid_idx in centroids {
            ids.extend(self.centroid_lists[centroid_idx].iter().copied());
        }
        self.limit_by_raw_dot(ids, sketch)
    }

    fn retrieve_hnsw(&self, sketch: &[f32]) -> Vec<usize> {
        let centroids = graph_centroid_search(
            sketch,
            &self.centroids,
            &self.centroid_graph,
            self.retrieval.probes,
        );
        let mut ids = Vec::new();
        for centroid_idx in centroids {
            ids.extend(self.centroid_lists[centroid_idx].iter().copied());
        }
        self.limit_by_raw_dot(ids, sketch)
    }

    fn limit_by_raw_dot(&self, ids: Vec<usize>, sketch: &[f32]) -> Vec<usize> {
        if ids.len() <= self.retrieval.candidate_limit {
            return ids;
        }
        let mut scored = ids
            .into_iter()
            .map(|idx| (idx, dot_unchecked(sketch, &self.windows[idx].sketch)))
            .collect::<Vec<_>>();
        scored.sort_by(|left, right| right.1.total_cmp(&left.1));
        scored
            .into_iter()
            .take(self.retrieval.candidate_limit)
            .map(|(idx, _)| idx)
            .collect()
    }

    fn window_count(&self) -> usize {
        self.windows.len()
    }

    fn compression_ratio(&self) -> f32 {
        let compressed = self
            .windows
            .iter()
            .map(|window| window.quantized.compressed_bytes())
            .sum::<usize>();
        if compressed == 0 {
            return 0.0;
        }
        (self.windows.len() * self.config.dim * std::mem::size_of::<f32>()) as f32
            / compressed as f32
    }
}

fn validate_retrieval_options(options: RetrievalOptions) -> Result<(), String> {
    if options.candidate_limit == 0 {
        return Err("candidate-limit must be non-zero".to_owned());
    }
    if options.minimizer_k == 0 || options.minimizer_k > 31 {
        return Err("minimizer-k must be in 1..=31".to_owned());
    }
    if options.minimizer_window == 0 {
        return Err("minimizer-window must be non-zero".to_owned());
    }
    if options.simhash_bands < 4 || options.simhash_bands > 16 || 64 % options.simhash_bands != 0 {
        return Err("simhash-bands must divide 64 and be in 4..=16".to_owned());
    }
    if options.centroids == 0 {
        return Err("centroids must be non-zero".to_owned());
    }
    if options.probes == 0 {
        return Err("probes must be non-zero".to_owned());
    }
    Ok(())
}

const MINIMIZER_CACHE_MAGIC: &[u8; 8] = b"DQMINI3\n";

fn load_minimizer_reference_cache(
    path: &str,
    reference_path: &str,
    config: ReferenceIndexConfig,
    retrieval: RetrievalOptions,
) -> Result<Option<(ExperimentalCandidateIndex, usize)>, String> {
    if retrieval.mode != RetrievalMode::Minimizer {
        return Err("--reference-cache currently supports only --retrieval minimizer".to_owned());
    }
    let Ok(mut file) = File::open(path) else {
        return Ok(None);
    };
    let mut magic = [0_u8; 8];
    file.read_exact(&mut magic)
        .map_err(|err| format!("failed to read reference cache magic {path}: {err}"))?;
    if &magic != MINIMIZER_CACHE_MAGIC {
        return Ok(None);
    }
    let cached_reference = read_string(&mut file)?;
    let cached_config = read_index_config(&mut file)?;
    let cached_retrieval = read_retrieval_options(&mut file)?;
    let reference_bases = read_usize(&mut file)?;
    if cached_reference != reference_path
        || !index_config_matches(cached_config, config)
        || cached_retrieval != retrieval
    {
        return Ok(None);
    }
    let target_count = read_usize(&mut file)?;
    let mut target_names = Vec::with_capacity(target_count);
    let mut target_offsets = Vec::with_capacity(target_count);
    for _ in 0..target_count {
        target_names.push(read_string(&mut file)?);
        target_offsets.push(read_usize(&mut file)?);
    }
    let window_count = read_usize(&mut file)?;
    let mut windows = Vec::with_capacity(window_count);
    for _ in 0..window_count {
        let target_id = read_u32(&mut file)? as usize;
        if target_id >= target_names.len() {
            return Err("reference cache window target id out of bounds".to_owned());
        }
        windows.push(CandidateWindow {
            target_id,
            target_start: read_u32(&mut file)? as usize,
            sketch: Vec::new(),
            quantized: read_cached_quantized_vector(&mut file, config)?,
            simhash64: 0,
            simhash128: [0, 0],
        });
    }
    let posting_count = read_usize(&mut file)?;
    let mut minimizer_postings = HashMap::with_capacity(posting_count);
    for _ in 0..posting_count {
        let hash = read_u64(&mut file)?;
        let ids_len = read_u32(&mut file)? as usize;
        let mut ids = Vec::with_capacity(ids_len);
        for _ in 0..ids_len {
            ids.push(read_u32(&mut file)? as usize);
        }
        minimizer_postings.insert(hash, ids);
    }
    Ok(Some((
        ExperimentalCandidateIndex {
            config,
            retrieval,
            target_names,
            target_offsets,
            windows,
            minimizer_postings,
            simhash_band_postings: Vec::new(),
            mih_postings: Vec::new(),
            centroids: Vec::new(),
            centroid_lists: Vec::new(),
            centroid_graph: Vec::new(),
        },
        reference_bases,
    )))
}

fn save_minimizer_reference_cache(
    path: &str,
    reference_path: &str,
    reference_bases: usize,
    index: &ExperimentalCandidateIndex,
) -> Result<(), String> {
    if index.retrieval.mode != RetrievalMode::Minimizer {
        return Ok(());
    }
    if let Some(parent) = std::path::Path::new(path).parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .map_err(|err| format!("failed to create reference cache directory: {err}"))?;
    }
    let tmp = format!("{path}.tmp");
    let mut file = File::create(&tmp)
        .map_err(|err| format!("failed to create reference cache {tmp}: {err}"))?;
    file.write_all(MINIMIZER_CACHE_MAGIC)
        .map_err(|err| format!("failed to write reference cache magic {tmp}: {err}"))?;
    write_string(&mut file, reference_path)?;
    write_index_config(&mut file, index.config)?;
    write_retrieval_options(&mut file, index.retrieval)?;
    write_usize(&mut file, reference_bases)?;
    write_usize(&mut file, index.target_names.len())?;
    for (name, offset) in index.target_names.iter().zip(&index.target_offsets) {
        write_string(&mut file, name)?;
        write_usize(&mut file, *offset)?;
    }
    write_usize(&mut file, index.windows.len())?;
    for window in &index.windows {
        write_u32_from_usize(&mut file, window.target_id)?;
        write_u32_from_usize(&mut file, window.target_start)?;
        write_cached_quantized_vector(&mut file, &window.quantized, index.config)?;
    }
    write_usize(&mut file, index.minimizer_postings.len())?;
    for (hash, ids) in &index.minimizer_postings {
        write_u64(&mut file, *hash)?;
        write_u32_from_usize(&mut file, ids.len())?;
        for &id in ids {
            write_u32_from_usize(&mut file, id)?;
        }
    }
    file.flush()
        .map_err(|err| format!("failed to flush reference cache {tmp}: {err}"))?;
    std::fs::rename(&tmp, path)
        .map_err(|err| format!("failed to install reference cache {path}: {err}"))?;
    Ok(())
}

fn index_config_matches(left: ReferenceIndexConfig, right: ReferenceIndexConfig) -> bool {
    left.k == right.k
        && left.dim == right.dim
        && left.window_len == right.window_len
        && left.stride == right.stride
        && left.quantizer.bits == right.quantizer.bits
        && left.quantizer.clip_sigma.to_bits() == right.quantizer.clip_sigma.to_bits()
        && left.quantizer.rotation_seed == right.quantizer.rotation_seed
        && left.quantizer.qjl_seed == right.quantizer.qjl_seed
        && left.quantizer.use_qjl_residual == right.quantizer.use_qjl_residual
}

fn write_index_config(writer: &mut impl Write, config: ReferenceIndexConfig) -> Result<(), String> {
    write_usize(writer, config.k)?;
    write_usize(writer, config.dim)?;
    write_usize(writer, config.window_len)?;
    write_usize(writer, config.stride)?;
    write_u8(writer, config.quantizer.bits)?;
    write_f32(writer, config.quantizer.clip_sigma)?;
    write_u64(writer, config.quantizer.rotation_seed)?;
    write_u64(writer, config.quantizer.qjl_seed)?;
    write_u8(writer, u8::from(config.quantizer.use_qjl_residual))
}

fn read_index_config(reader: &mut impl Read) -> Result<ReferenceIndexConfig, String> {
    Ok(ReferenceIndexConfig {
        k: read_usize(reader)?,
        dim: read_usize(reader)?,
        window_len: read_usize(reader)?,
        stride: read_usize(reader)?,
        quantizer: QuantizerConfig {
            bits: read_u8(reader)?,
            clip_sigma: read_f32(reader)?,
            rotation_seed: read_u64(reader)?,
            qjl_seed: read_u64(reader)?,
            use_qjl_residual: read_u8(reader)? != 0,
        },
    })
}

fn write_retrieval_options(
    writer: &mut impl Write,
    retrieval: RetrievalOptions,
) -> Result<(), String> {
    write_u8(writer, retrieval_mode_code(retrieval.mode))?;
    write_usize(writer, retrieval.candidate_limit)?;
    write_usize(writer, retrieval.minimizer_k)?;
    write_usize(writer, retrieval.minimizer_window)?;
    write_usize(writer, retrieval.simhash_bands)?;
    write_usize(writer, retrieval.centroids)?;
    write_usize(writer, retrieval.probes)
}

fn read_retrieval_options(reader: &mut impl Read) -> Result<RetrievalOptions, String> {
    Ok(RetrievalOptions {
        mode: retrieval_mode_from_code(read_u8(reader)?)?,
        candidate_limit: read_usize(reader)?,
        minimizer_k: read_usize(reader)?,
        minimizer_window: read_usize(reader)?,
        simhash_bands: read_usize(reader)?,
        centroids: read_usize(reader)?,
        probes: read_usize(reader)?,
    })
}

fn retrieval_mode_code(mode: RetrievalMode) -> u8 {
    match mode {
        RetrievalMode::Scan => 0,
        RetrievalMode::Minimizer => 1,
        RetrievalMode::Simhash => 2,
        RetrievalMode::Mih => 3,
        RetrievalMode::Ivf => 4,
        RetrievalMode::Hnsw => 5,
    }
}

fn retrieval_mode_from_code(code: u8) -> Result<RetrievalMode, String> {
    match code {
        0 => Ok(RetrievalMode::Scan),
        1 => Ok(RetrievalMode::Minimizer),
        2 => Ok(RetrievalMode::Simhash),
        3 => Ok(RetrievalMode::Mih),
        4 => Ok(RetrievalMode::Ivf),
        5 => Ok(RetrievalMode::Hnsw),
        _ => Err("unknown retrieval mode in reference cache".to_owned()),
    }
}

fn write_cached_quantized_vector(
    writer: &mut impl Write,
    vector: &QuantizedVector,
    config: ReferenceIndexConfig,
) -> Result<(), String> {
    let snapshot = vector.snapshot();
    if snapshot.dim != config.dim
        || snapshot.bits != config.quantizer.bits
        || snapshot.codes_len != config.dim
        || snapshot.codes_bits != config.quantizer.bits
        || snapshot.rotation_seed != config.quantizer.rotation_seed
    {
        return Err("quantized vector does not match reference cache config".to_owned());
    }
    let expected_code_bytes = (config.dim * usize::from(config.quantizer.bits)).div_ceil(8);
    if snapshot.codes_data.len() != expected_code_bytes {
        return Err(
            "quantized vector code byte length does not match reference cache config".to_owned(),
        );
    }
    writer
        .write_all(&snapshot.codes_data)
        .map_err(|err| format!("failed to write reference cache quantized codes: {err}"))?;
    if config.quantizer.use_qjl_residual {
        match snapshot.qjl_residual {
            Some(residual) => {
                write_u8(writer, 1)?;
                write_usize(writer, residual.signs.len())?;
                for sign in residual.signs {
                    write_u64(writer, sign)?;
                }
                write_f32(writer, residual.norm)?;
            }
            None => write_u8(writer, 0)?,
        }
    }
    Ok(())
}

fn read_cached_quantized_vector(
    reader: &mut impl Read,
    config: ReferenceIndexConfig,
) -> Result<QuantizedVector, String> {
    let code_bytes = (config.dim * usize::from(config.quantizer.bits)).div_ceil(8);
    let mut codes_data = vec![0_u8; code_bytes];
    reader
        .read_exact(&mut codes_data)
        .map_err(|err| format!("failed to read reference cache quantized codes: {err}"))?;
    let qjl_residual = if config.quantizer.use_qjl_residual {
        if read_u8(reader)? == 0 {
            None
        } else {
            let signs_len = read_usize(reader)?;
            let mut signs = Vec::with_capacity(signs_len);
            for _ in 0..signs_len {
                signs.push(read_u64(reader)?);
            }
            Some(QjlResidualSnapshot {
                signs,
                dim: config.dim,
                norm: read_f32(reader)?,
                seed: config.quantizer.qjl_seed,
            })
        }
    } else {
        None
    };
    QuantizedVector::from_snapshot(QuantizedVectorSnapshot {
        dim: config.dim,
        bits: config.quantizer.bits,
        clip: config.quantizer.clip_sigma / (config.dim as f32).sqrt(),
        rotation_seed: config.quantizer.rotation_seed,
        codes_data,
        codes_len: config.dim,
        codes_bits: config.quantizer.bits,
        qjl_residual,
    })
}

fn write_string(writer: &mut impl Write, value: &str) -> Result<(), String> {
    write_bytes(writer, value.as_bytes())
}

fn read_string(reader: &mut impl Read) -> Result<String, String> {
    String::from_utf8(read_bytes(reader)?)
        .map_err(|err| format!("reference cache contained invalid UTF-8: {err}"))
}

fn write_bytes(writer: &mut impl Write, value: &[u8]) -> Result<(), String> {
    write_usize(writer, value.len())?;
    writer
        .write_all(value)
        .map_err(|err| format!("failed to write reference cache bytes: {err}"))
}

fn read_bytes(reader: &mut impl Read) -> Result<Vec<u8>, String> {
    let len = read_usize(reader)?;
    let mut value = vec![0_u8; len];
    reader
        .read_exact(&mut value)
        .map_err(|err| format!("failed to read reference cache bytes: {err}"))?;
    Ok(value)
}

fn write_usize(writer: &mut impl Write, value: usize) -> Result<(), String> {
    write_u64(writer, value as u64)
}

fn read_usize(reader: &mut impl Read) -> Result<usize, String> {
    usize::try_from(read_u64(reader)?)
        .map_err(|_| "reference cache integer does not fit usize".to_owned())
}

fn write_u32_from_usize(writer: &mut impl Write, value: usize) -> Result<(), String> {
    let value =
        u32::try_from(value).map_err(|_| "reference cache integer does not fit u32".to_owned())?;
    write_u32(writer, value)
}

fn write_u32(writer: &mut impl Write, value: u32) -> Result<(), String> {
    writer
        .write_all(&value.to_le_bytes())
        .map_err(|err| format!("failed to write reference cache u32: {err}"))
}

fn read_u32(reader: &mut impl Read) -> Result<u32, String> {
    let mut bytes = [0_u8; 4];
    reader
        .read_exact(&mut bytes)
        .map_err(|err| format!("failed to read reference cache u32: {err}"))?;
    Ok(u32::from_le_bytes(bytes))
}

fn write_u64(writer: &mut impl Write, value: u64) -> Result<(), String> {
    writer
        .write_all(&value.to_le_bytes())
        .map_err(|err| format!("failed to write reference cache u64: {err}"))
}

fn read_u64(reader: &mut impl Read) -> Result<u64, String> {
    let mut bytes = [0_u8; 8];
    reader
        .read_exact(&mut bytes)
        .map_err(|err| format!("failed to read reference cache u64: {err}"))?;
    Ok(u64::from_le_bytes(bytes))
}

fn write_f32(writer: &mut impl Write, value: f32) -> Result<(), String> {
    writer
        .write_all(&value.to_le_bytes())
        .map_err(|err| format!("failed to write reference cache f32: {err}"))
}

fn read_f32(reader: &mut impl Read) -> Result<f32, String> {
    let mut bytes = [0_u8; 4];
    reader
        .read_exact(&mut bytes)
        .map_err(|err| format!("failed to read reference cache f32: {err}"))?;
    Ok(f32::from_le_bytes(bytes))
}

fn write_u8(writer: &mut impl Write, value: u8) -> Result<(), String> {
    writer
        .write_all(&[value])
        .map_err(|err| format!("failed to write reference cache u8: {err}"))
}

fn read_u8(reader: &mut impl Read) -> Result<u8, String> {
    let mut byte = [0_u8; 1];
    reader
        .read_exact(&mut byte)
        .map_err(|err| format!("failed to read reference cache u8: {err}"))?;
    Ok(byte[0])
}

#[derive(Clone, Debug)]
struct ProteinEmitOptions {
    database_path: String,
    query_path: String,
    top_k: usize,
    k: usize,
    dim: usize,
    bits: u8,
    max_bases: Option<usize>,
    use_qjl_residual: bool,
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
    let config = build_index_config(options.bench)?;
    if options.retrieval.mode != RetrievalMode::Scan
        && let Some(cache_path) = options.reference_cache.as_deref()
        && let Some((index, reference_bases)) = load_minimizer_reference_cache(
            cache_path,
            &options.reference_path,
            config,
            options.retrieval,
        )?
    {
        return run_emit_candidates_with_index(
            options,
            index,
            reference_bases,
            Duration::ZERO,
            true,
        );
    }

    let mut records = read_fasta_file(&options.reference_path)?;
    truncate_records(&mut records, options.bench.max_bases);
    let reference_bases = records
        .iter()
        .map(|record| record.bases.len())
        .sum::<usize>();
    if options.retrieval.mode == RetrievalMode::Scan {
        return run_emit_candidates_scan(options, records, config, reference_bases);
    }

    run_emit_candidates_indexed(options, records, config, reference_bases)
}

fn run_emit_candidates_scan(
    options: EmitOptions,
    records: Vec<dino_quant::SequenceRecord>,
    config: ReferenceIndexConfig,
    reference_bases: usize,
) -> Result<(), String> {
    let build_start = Instant::now();
    let index = ReferenceWindowIndex::build_records(&records, config)?;
    let build_elapsed = build_start.elapsed();

    println!(
        "query_name\tquery_len\trank\ttarget_name\ttarget_start\ttarget_end\tlinear_start\tlinear_end\tscore"
    );
    let mut reads = 0_usize;
    let mut candidates = 0_usize;
    let search_start = Instant::now();
    visit_fastq_slices_file_limit(&options.reads_path, options.max_reads, |record| {
        let hits = index.search_sequence(record.bases, options.bench.top_k)?;
        reads += 1;
        let read_name = String::from_utf8_lossy(record.name);
        for (rank, hit) in hits.iter().enumerate() {
            candidates += 1;
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
                read_name,
                record.bases.len(),
                rank + 1,
                hit.target_name,
                hit.target_start,
                hit.target_end,
                hit.start,
                hit.end,
                hit.score
            );
        }
        Ok(())
    })?;
    let search_elapsed = search_start.elapsed();

    eprintln!(
        "dino-quant emit-candidates: retrieval=scan reference_bases={} windows={} reads={} candidates={} considered={} avg_considered={:.2} build_ms={:.2} read_qps={:.2} compression_ratio={:.2}",
        reference_bases,
        index.window_count(),
        reads,
        candidates,
        index.window_count() * reads,
        index.window_count() as f32,
        ms(build_elapsed),
        qps(reads, search_elapsed),
        index.compression_ratio()
    );
    Ok(())
}

fn run_emit_candidates_indexed(
    options: EmitOptions,
    records: Vec<dino_quant::SequenceRecord>,
    config: ReferenceIndexConfig,
    reference_bases: usize,
) -> Result<(), String> {
    let build_start = Instant::now();
    let index = ExperimentalCandidateIndex::build_records(&records, config, options.retrieval)?;
    let build_elapsed = build_start.elapsed();
    if let Some(cache_path) = options.reference_cache.as_deref() {
        save_minimizer_reference_cache(
            cache_path,
            &options.reference_path,
            reference_bases,
            &index,
        )?;
    }

    run_emit_candidates_with_index(options, index, reference_bases, build_elapsed, false)
}

fn run_emit_candidates_with_index(
    options: EmitOptions,
    index: ExperimentalCandidateIndex,
    reference_bases: usize,
    build_elapsed: Duration,
    cache_hit: bool,
) -> Result<(), String> {
    println!(
        "query_name\tquery_len\trank\ttarget_name\ttarget_start\ttarget_end\tlinear_start\tlinear_end\tscore"
    );
    let mut reads = 0_usize;
    let mut candidates = 0_usize;
    let mut considered = 0_usize;
    let mut scratch = CandidateScratch::new(index.window_count());
    let search_start = Instant::now();
    visit_fastq_slices_file_limit(&options.reads_path, options.max_reads, |record| {
        let result =
            index.search_sequence_with_scratch(record.bases, options.bench.top_k, &mut scratch)?;
        considered += result.considered;
        reads += 1;
        let read_name = String::from_utf8_lossy(record.name);
        for (rank, hit) in result.hits.iter().enumerate() {
            candidates += 1;
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
                read_name,
                record.bases.len(),
                rank + 1,
                hit.target_name,
                hit.target_start,
                hit.target_end,
                hit.start,
                hit.end,
                hit.score
            );
        }
        Ok(())
    })?;
    let search_elapsed = search_start.elapsed();

    eprintln!(
        "dino-quant emit-candidates: retrieval={} reference_bases={} windows={} reads={} candidates={} considered={} avg_considered={:.2} build_ms={:.2} cache_hit={} read_qps={:.2} compression_ratio={:.2}",
        options.retrieval.mode.as_str(),
        reference_bases,
        index.window_count(),
        reads,
        candidates,
        considered,
        considered as f32 / reads.max(1) as f32,
        ms(build_elapsed),
        cache_hit,
        qps(reads, search_elapsed),
        index.compression_ratio()
    );
    Ok(())
}

fn run_emit_candidate_reference(args: Vec<String>) -> Result<(), String> {
    let options = parse_candidate_reference_args(&args)?;
    let records = read_fasta_file(&options.reference_path)?;
    let intervals = read_candidate_intervals(&options.candidates_path, options.padding)?;
    let mut emitted_records = 0_usize;
    let mut emitted_bases = 0_usize;
    let mut input_intervals = 0_usize;
    let mut merged_intervals = 0_usize;

    for record in &records {
        let Some(record_intervals) = intervals.get(&record.name) else {
            if options.mask_reference {
                emit_masked_fasta_record(&record.name, &record.bases, &[]);
            }
            continue;
        };
        input_intervals += record_intervals.len();
        let bounded = record_intervals
            .iter()
            .filter_map(|interval| bound_interval(*interval, record.bases.len()))
            .collect::<Vec<_>>();
        let merged = merge_intervals(bounded, options.merge_gap);
        merged_intervals += merged.len();
        if options.mask_reference {
            emit_masked_fasta_record(&record.name, &record.bases, &merged);
            emitted_records += 1;
            emitted_bases += record.bases.len();
            continue;
        }
        for interval in merged {
            let header = format!("{}:{}-{}", record.name, interval.start, interval.end);
            emit_fasta_record(&header, &record.bases[interval.start..interval.end]);
            emitted_records += 1;
            emitted_bases += interval.end - interval.start;
        }
    }

    eprintln!(
        "dino-quant emit-candidate-reference: input_intervals={} merged_intervals={} records={} bases={}",
        input_intervals, merged_intervals, emitted_records, emitted_bases
    );
    Ok(())
}

fn read_candidate_intervals(
    path: &str,
    padding: usize,
) -> Result<HashMap<String, Vec<Interval>>, String> {
    let input = std::fs::read_to_string(path)
        .map_err(|err| format!("failed to read candidate TSV {path}: {err}"))?;
    let mut intervals: HashMap<String, Vec<Interval>> = HashMap::new();
    for (line_idx, line) in input.lines().enumerate() {
        if line_idx == 0 && line.starts_with("query_name\t") {
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 6 {
            return Err(format!(
                "candidate TSV line {} has fewer than 6 fields",
                line_idx + 1
            ));
        }
        let target_name = fields[3].to_owned();
        let start = fields[4]
            .parse::<usize>()
            .map_err(|err| format!("invalid target_start on line {}: {err}", line_idx + 1))?;
        let end = fields[5]
            .parse::<usize>()
            .map_err(|err| format!("invalid target_end on line {}: {err}", line_idx + 1))?;
        if start >= end {
            continue;
        }
        intervals.entry(target_name).or_default().push(Interval {
            start: start.saturating_sub(padding),
            end: end.saturating_add(padding),
        });
    }
    Ok(intervals)
}

fn bound_interval(interval: Interval, len: usize) -> Option<Interval> {
    let start = interval.start.min(len);
    let end = interval.end.min(len);
    (start < end).then_some(Interval { start, end })
}

fn merge_intervals(mut intervals: Vec<Interval>, merge_gap: usize) -> Vec<Interval> {
    if intervals.is_empty() {
        return intervals;
    }
    intervals.sort_by_key(|interval| (interval.start, interval.end));
    let mut merged = Vec::with_capacity(intervals.len());
    let mut current = intervals[0];
    for interval in intervals.into_iter().skip(1) {
        if interval.start <= current.end.saturating_add(merge_gap) {
            current.end = current.end.max(interval.end);
        } else {
            merged.push(current);
            current = interval;
        }
    }
    merged.push(current);
    merged
}

fn emit_fasta_record(header: &str, bases: &[u8]) {
    println!(">{header}");
    for chunk in bases.chunks(80) {
        println!("{}", String::from_utf8_lossy(chunk));
    }
}

fn emit_masked_fasta_record(header: &str, bases: &[u8], intervals: &[Interval]) {
    println!(">{header}");
    let mut line = Vec::with_capacity(80);
    let mut interval_idx = 0_usize;
    for (pos, &base) in bases.iter().enumerate() {
        while interval_idx < intervals.len() && intervals[interval_idx].end <= pos {
            interval_idx += 1;
        }
        let keep = interval_idx < intervals.len()
            && intervals[interval_idx].start <= pos
            && pos < intervals[interval_idx].end;
        line.push(if keep { base } else { b'N' });
        if line.len() == 80 {
            println!("{}", String::from_utf8_lossy(&line));
            line.clear();
        }
    }
    if !line.is_empty() {
        println!("{}", String::from_utf8_lossy(&line));
    }
}

fn run_emit_protein_candidates(args: Vec<String>) -> Result<(), String> {
    let options = parse_protein_emit_args(&args)?;
    let mut database = read_protein_fasta_file(&options.database_path)?;
    truncate_records(&mut database, options.max_bases);
    let queries = read_protein_fasta_file(&options.query_path)?;
    let quantizer = QuantizerConfig {
        bits: options.bits,
        use_qjl_residual: options.use_qjl_residual,
        ..QuantizerConfig::default()
    };
    quantizer.validate_for_cli(options.dim)?;

    let build_start = Instant::now();
    let mut targets = Vec::with_capacity(database.len());
    for record in &database {
        let sketch = protein_kmer_sketch(&record.bases, options.k, options.dim)?;
        let sketch = quantizer.encode(&sketch)?;
        targets.push((&record.name, record.bases.len(), sketch));
    }
    let build_elapsed = build_start.elapsed();

    println!("query_name\tquery_len\trank\ttarget_name\ttarget_len\tscore");
    let search_start = Instant::now();
    let mut candidates = 0_usize;
    for query in &queries {
        let query_sketch = protein_kmer_sketch(&query.bases, options.k, options.dim)?;
        let mut top = Vec::with_capacity(options.top_k.min(targets.len()));
        for (name, len, sketch) in &targets {
            push_protein_hit(
                &mut top,
                options.top_k,
                ProteinHit {
                    target_name: (*name).clone(),
                    target_len: *len,
                    score: sketch.approximate_dot_query(&query_sketch)?,
                },
            );
        }
        top.sort_by(|left, right| right.score.total_cmp(&left.score));
        for (rank, hit) in top.iter().enumerate() {
            candidates += 1;
            println!(
                "{}\t{}\t{}\t{}\t{}\t{:.6}",
                query.name,
                query.bases.len(),
                rank + 1,
                hit.target_name,
                hit.target_len,
                hit.score
            );
        }
    }
    let search_elapsed = search_start.elapsed();
    eprintln!(
        "dino-quant emit-protein-candidates: targets={} queries={} candidates={} build_ms={:.2} query_qps={:.2}",
        targets.len(),
        queries.len(),
        candidates,
        ms(build_elapsed),
        qps(queries.len(), search_elapsed)
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
    let mut max_reads = None;
    let mut retrieval = RetrievalOptions::default();
    let mut reference_cache = None;
    let mut idx = 2;
    while idx < args.len() {
        match args[idx].as_str() {
            "--retrieval" => {
                let value = args
                    .get(idx + 1)
                    .ok_or_else(|| "--retrieval requires a value".to_owned())?;
                retrieval.mode = RetrievalMode::parse(value)?;
                idx += 2;
            }
            "--candidate-limit" => {
                retrieval.candidate_limit = parse_next_usize(args, &mut idx, "--candidate-limit")?;
            }
            "--reference-cache" => {
                reference_cache = Some(
                    args.get(idx + 1)
                        .ok_or_else(|| "--reference-cache requires a path".to_owned())?
                        .clone(),
                );
                idx += 2;
            }
            "--minimizer-k" => {
                retrieval.minimizer_k = parse_next_usize(args, &mut idx, "--minimizer-k")?;
            }
            "--minimizer-window" => {
                retrieval.minimizer_window =
                    parse_next_usize(args, &mut idx, "--minimizer-window")?;
            }
            "--simhash-bands" => {
                retrieval.simhash_bands = parse_next_usize(args, &mut idx, "--simhash-bands")?;
            }
            "--centroids" => {
                retrieval.centroids = parse_next_usize(args, &mut idx, "--centroids")?;
            }
            "--probes" => {
                retrieval.probes = parse_next_usize(args, &mut idx, "--probes")?;
            }
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
            "--max-reads" => {
                max_reads = Some(parse_next_usize(args, &mut idx, "--max-reads")?);
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
        max_reads,
        bench,
        retrieval,
        reference_cache,
    })
}

fn parse_candidate_reference_args(args: &[String]) -> Result<CandidateReferenceOptions, String> {
    let reference_path = args
        .first()
        .ok_or_else(|| "emit-candidate-reference requires a reference FASTA path".to_owned())?
        .clone();
    let candidates_path = args
        .get(1)
        .ok_or_else(|| "emit-candidate-reference requires a candidate TSV path".to_owned())?
        .clone();
    let mut options = CandidateReferenceOptions {
        reference_path,
        candidates_path,
        padding: 64,
        merge_gap: 0,
        mask_reference: false,
    };
    let mut idx = 2;
    while idx < args.len() {
        match args[idx].as_str() {
            "--padding" => {
                options.padding = parse_next_usize(args, &mut idx, "--padding")?;
            }
            "--merge-gap" => {
                options.merge_gap = parse_next_usize(args, &mut idx, "--merge-gap")?;
            }
            "--mask-reference" => {
                options.mask_reference = true;
                idx += 1;
            }
            other => return Err(format!("unknown emit-candidate-reference option: {other}")),
        }
    }
    Ok(options)
}

fn parse_protein_emit_args(args: &[String]) -> Result<ProteinEmitOptions, String> {
    let database_path = args
        .first()
        .ok_or_else(|| "emit-protein-candidates requires a protein FASTA database".to_owned())?
        .clone();
    let query_path = args
        .get(1)
        .ok_or_else(|| "emit-protein-candidates requires a protein FASTA query file".to_owned())?
        .clone();
    let mut options = ProteinEmitOptions {
        database_path,
        query_path,
        top_k: 10,
        k: 5,
        dim: 256,
        bits: 4,
        max_bases: None,
        use_qjl_residual: false,
    };
    let mut idx = 2;
    while idx < args.len() {
        match args[idx].as_str() {
            "--top-k" => {
                options.top_k = parse_next_usize(args, &mut idx, "--top-k")?;
            }
            "--k" => {
                options.k = parse_next_usize(args, &mut idx, "--k")?;
            }
            "--dim" => {
                options.dim = parse_next_usize(args, &mut idx, "--dim")?;
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
            other => return Err(format!("unknown emit-protein-candidates option: {other}")),
        }
    }
    if options.top_k == 0 {
        return Err("top-k must be non-zero".to_owned());
    }
    Ok(options)
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

fn truncate_records(records: &mut Vec<dino_quant::SequenceRecord>, max_bases: Option<usize>) {
    let Some(mut remaining) = max_bases else {
        return;
    };
    let mut keep = 0_usize;
    for record in records.iter_mut() {
        if remaining == 0 {
            break;
        }
        if record.bases.len() > remaining {
            record.bases.truncate(remaining);
        }
        remaining = remaining.saturating_sub(record.bases.len());
        keep += 1;
    }
    records.truncate(keep);
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
                target_name: "reference".to_owned(),
                target_start: window.start,
                target_end: window.end,
                start: window.start,
                end: window.end,
                score: dot(&window.sketch, query_sketch)?,
            },
        );
    }
    top.sort_by(|left, right| right.score.total_cmp(&left.score));
    Ok(top)
}

#[derive(Clone, Debug)]
struct ProteinHit {
    target_name: String,
    target_len: usize,
    score: f32,
}

fn push_protein_hit(top: &mut Vec<ProteinHit>, top_k: usize, hit: ProteinHit) {
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

#[cfg(test)]
fn minimizer_hashes(seq: &[u8], k: usize, window: usize) -> Result<Vec<u64>, String> {
    let mut hashes = Vec::new();
    let mut minimizers = Vec::new();
    let mut deque = Vec::new();
    minimizer_hashes_into(seq, k, window, &mut hashes, &mut minimizers, &mut deque)?;
    Ok(minimizers)
}

fn minimizer_hashes_into(
    seq: &[u8],
    k: usize,
    window: usize,
    hashes: &mut Vec<u64>,
    minimizers: &mut Vec<u64>,
    deque: &mut Vec<usize>,
) -> Result<(), String> {
    hashes.clear();
    minimizers.clear();
    if k == 0 || k > 31 {
        return Err("minimizer k must be in 1..=31".to_owned());
    }
    if window == 0 {
        return Err("minimizer window must be non-zero".to_owned());
    }
    if seq.len() < k {
        return Ok(());
    }

    canonical_kmer_hashes_into(seq, k, hashes);
    if hashes.is_empty() {
        return Ok(());
    }

    let span = window.min(hashes.len());
    let mut last = None;
    if span <= 16 {
        for start in 0..=hashes.len() - span {
            let mut best = hashes[start];
            for &hash in &hashes[start + 1..start + span] {
                if hash < best {
                    best = hash;
                }
            }
            if last != Some(best) {
                minimizers.push(best);
                last = Some(best);
            }
        }
        return Ok(());
    }

    deque.clear();
    let mut head = 0_usize;
    for (idx, &hash) in hashes.iter().enumerate() {
        while deque.len() > head && hashes[*deque.last().unwrap_or(&idx)] > hash {
            deque.pop();
        }
        deque.push(idx);
        if idx + 1 < span {
            continue;
        }
        let start = idx + 1 - span;
        while head < deque.len() && deque[head] < start {
            head += 1;
        }
        let Some(&best_idx) = deque.get(head) else {
            return Err("minimizer deque unexpectedly emptied".to_owned());
        };
        let best = hashes[best_idx];
        if last != Some(best) {
            minimizers.push(best);
            last = Some(best);
        }
    }
    Ok(())
}

#[cfg(test)]
fn canonical_kmer_hash(kmer: &[u8]) -> Option<u64> {
    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    for (idx, &base) in kmer.iter().enumerate() {
        let code = u64::from(dna_base_code(base)?);
        forward = (forward << 2) | code;
        let rc = 3_u64 - code;
        reverse |= rc << (idx * 2);
    }
    Some(mix_u64_local(
        forward.min(reverse) ^ ((kmer.len() as u64) << 56),
    ))
}

fn canonical_kmer_hashes_into(seq: &[u8], k: usize, hashes: &mut Vec<u64>) {
    hashes.reserve(seq.len() + 1 - k);
    let mask = (1_u64 << (2 * k)) - 1;
    let rc_shift = 2 * (k - 1);
    let mut valid = 0_usize;
    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    for &base in seq {
        let Some(code) = dna_base_code(base).map(u64::from) else {
            valid = 0;
            forward = 0;
            reverse = 0;
            continue;
        };
        forward = ((forward << 2) | code) & mask;
        reverse = (reverse >> 2) | ((3 - code) << rc_shift);
        valid += 1;
        if valid >= k {
            hashes.push(mix_u64_local(forward.min(reverse) ^ ((k as u64) << 56)));
        }
    }
}

fn dna_base_code(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn simhash_code64(sketch: &[f32], seed: u64) -> u64 {
    let mut accum = [0.0_f32; 64];
    for (dim, &value) in sketch.iter().enumerate() {
        if value == 0.0 {
            continue;
        }
        let hash = mix_u64_local(seed ^ dim as u64);
        let bucket = (hash as usize) & 63;
        let sign = if (hash >> 63) == 0 { 1.0 } else { -1.0 };
        accum[bucket] += sign * value;
    }
    let mut code = 0_u64;
    for (bit, &value) in accum.iter().enumerate() {
        if value >= 0.0 {
            code |= 1_u64 << bit;
        }
    }
    code
}

fn simhash_band_key(code: u64, band: usize, bands: usize) -> u16 {
    let bits = 64 / bands;
    let mask = if bits == 16 {
        u16::MAX
    } else {
        ((1_u32 << bits) - 1) as u16
    };
    ((code >> (band * bits)) as u16) & mask
}

fn simhash128_part(code: [u64; 2], part: usize) -> u16 {
    if part < 4 {
        ((code[0] >> (part * 16)) & 0xffff) as u16
    } else {
        ((code[1] >> ((part - 4) * 16)) & 0xffff) as u16
    }
}

fn bump_count(counts: &mut HashMap<usize, u32>, idx: usize, weight: u32) {
    counts
        .entry(idx)
        .and_modify(|count| *count = count.saturating_add(weight))
        .or_insert(weight);
}

fn top_counted_ids(counts: HashMap<usize, u32>, limit: usize) -> Vec<usize> {
    let mut ranked = counts.into_iter().collect::<Vec<_>>();
    ranked.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    ranked.into_iter().take(limit).map(|(idx, _)| idx).collect()
}

fn nearest_centroid(query: &[f32], centroids: &[Vec<f32>]) -> usize {
    top_centroids(query, centroids, 1)
        .into_iter()
        .next()
        .unwrap_or(0)
}

fn top_centroids(query: &[f32], centroids: &[Vec<f32>], probes: usize) -> Vec<usize> {
    let mut scores = centroids
        .iter()
        .enumerate()
        .map(|(idx, centroid)| (idx, dot_unchecked(query, centroid)))
        .collect::<Vec<_>>();
    scores.sort_by(|left, right| right.1.total_cmp(&left.1));
    scores
        .into_iter()
        .take(probes.min(centroids.len()))
        .map(|(idx, _)| idx)
        .collect()
}

fn graph_centroid_search(
    query: &[f32],
    centroids: &[Vec<f32>],
    graph: &[Vec<usize>],
    probes: usize,
) -> Vec<usize> {
    if centroids.is_empty() {
        return Vec::new();
    }
    let mut seen = vec![false; centroids.len()];
    let mut frontier = vec![0_usize];
    seen[0] = true;
    let mut scored = Vec::new();
    while let Some(idx) = frontier.pop() {
        scored.push((idx, dot_unchecked(query, &centroids[idx])));
        for &next in &graph[idx] {
            if !seen[next] {
                seen[next] = true;
                frontier.push(next);
            }
        }
        if scored.len() >= probes.saturating_mul(8).max(probes) {
            break;
        }
    }
    scored.sort_by(|left, right| right.1.total_cmp(&left.1));
    scored
        .into_iter()
        .take(probes.min(centroids.len()))
        .map(|(idx, _)| idx)
        .collect()
}

fn dot_unchecked(left: &[f32], right: &[f32]) -> f32 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn mix_u64_local(mut value: u64) -> u64 {
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
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

trait QuantizerCliValidation {
    fn validate_for_cli(&self, dim: usize) -> Result<(), String>;
}

impl QuantizerCliValidation for QuantizerConfig {
    fn validate_for_cli(&self, dim: usize) -> Result<(), String> {
        self.encode(&vec![0.0; dim]).map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_config(window_len: usize, stride: usize) -> ReferenceIndexConfig {
        ReferenceIndexConfig {
            k: 3,
            dim: 32,
            window_len,
            stride,
            quantizer: QuantizerConfig {
                use_qjl_residual: false,
                ..QuantizerConfig::default()
            },
        }
    }

    #[test]
    fn counted_ids_prefer_higher_counts_then_lower_index() {
        let mut counts = HashMap::new();
        bump_count(&mut counts, 9, 1);
        bump_count(&mut counts, 4, 2);
        bump_count(&mut counts, 2, 2);

        assert_eq!(top_counted_ids(counts, 3), vec![2, 4, 9]);
    }

    #[test]
    fn monotonic_minimizers_match_naive_window_scan() -> Result<(), String> {
        let seq = b"ACGTACGTNNNNGGGGTTTTAAAACCCCACGTACGT";
        for window in [1, 2, 3, 8, 100] {
            let fast = minimizer_hashes(seq, 3, window)?;
            let naive = naive_minimizer_hashes(seq, 3, window)?;
            assert_eq!(fast, naive, "window={window}");
        }
        Ok(())
    }

    #[test]
    fn rolling_kmer_hashes_match_window_hashes() {
        let seq = b"ACGTACGTNNNNGGGGTTTTAAAACCCCACGTACGT";
        for k in [1, 3, 15, 31] {
            let mut rolling = Vec::new();
            canonical_kmer_hashes_into(seq, k, &mut rolling);
            let windowed = seq
                .windows(k)
                .filter_map(canonical_kmer_hash)
                .collect::<Vec<_>>();
            assert_eq!(rolling, windowed, "k={k}");
        }
    }

    fn naive_minimizer_hashes(seq: &[u8], k: usize, window: usize) -> Result<Vec<u64>, String> {
        if k == 0 || k > 31 {
            return Err("minimizer k must be in 1..=31".to_owned());
        }
        if window == 0 {
            return Err("minimizer window must be non-zero".to_owned());
        }
        if seq.len() < k {
            return Ok(Vec::new());
        }

        let hashes = seq
            .windows(k)
            .filter_map(canonical_kmer_hash)
            .collect::<Vec<_>>();
        if hashes.is_empty() {
            return Ok(Vec::new());
        }

        let span = window.min(hashes.len());
        let mut minimizers = Vec::new();
        let mut last = None;
        for start in 0..=hashes.len() - span {
            let mut best = hashes[start];
            for &hash in &hashes[start + 1..start + span] {
                if hash < best {
                    best = hash;
                }
            }
            if last != Some(best) {
                minimizers.push(best);
                last = Some(best);
            }
        }
        Ok(minimizers)
    }

    #[test]
    fn candidate_scratch_reuses_dense_counts_with_same_ranking() -> Result<(), String> {
        let mut scratch = CandidateScratch::new(10);
        scratch.bump(9, 1)?;
        scratch.bump(4, 2)?;
        scratch.bump(2, 2)?;
        assert_eq!(scratch.ranked_candidate_ids(3), &[2, 4, 9]);

        scratch.bump(1, 3)?;
        scratch.bump(9, 1)?;
        assert_eq!(scratch.ranked_candidate_ids(2), &[1, 9]);
        Ok(())
    }

    #[test]
    fn minimizer_retrieval_finds_exact_window() -> Result<(), String> {
        let records = vec![dino_quant::SequenceRecord {
            name: "ref".to_owned(),
            bases: b"ACGTACGTGGGGTTTT".to_vec(),
        }];
        let retrieval = RetrievalOptions {
            mode: RetrievalMode::Minimizer,
            candidate_limit: 16,
            minimizer_k: 3,
            minimizer_window: 2,
            ..RetrievalOptions::default()
        };
        let index =
            ExperimentalCandidateIndex::build_records(&records, test_config(8, 8), retrieval)?;

        let result = index.search_sequence(b"ACGTACGT", 4)?;

        assert!(result.considered > 0);
        assert_eq!(result.hits.first().map(|hit| hit.target_start), Some(0));
        Ok(())
    }

    #[test]
    fn centroid_retrieval_honors_candidate_limit() -> Result<(), String> {
        let records = vec![dino_quant::SequenceRecord {
            name: "ref".to_owned(),
            bases: b"ACGTACGTGGGGTTTTAAAACCCC".to_vec(),
        }];
        let retrieval = RetrievalOptions {
            mode: RetrievalMode::Ivf,
            candidate_limit: 2,
            centroids: 1,
            probes: 1,
            ..RetrievalOptions::default()
        };
        let index =
            ExperimentalCandidateIndex::build_records(&records, test_config(8, 4), retrieval)?;

        let result = index.search_sequence(b"ACGTACGT", 4)?;

        assert!(result.considered <= 2);
        assert!(result.hits.len() <= 2);
        Ok(())
    }

    #[test]
    fn candidate_reference_intervals_are_bounded_and_merged() {
        let intervals = vec![
            Interval { start: 8, end: 12 },
            Interval { start: 0, end: 4 },
            Interval { start: 5, end: 7 },
            Interval { start: 30, end: 50 },
        ];
        let bounded = intervals
            .into_iter()
            .filter_map(|interval| bound_interval(interval, 32))
            .collect::<Vec<_>>();

        assert_eq!(
            merge_intervals(bounded, 1),
            vec![
                Interval { start: 0, end: 12 },
                Interval { start: 30, end: 32 }
            ]
        );
    }
}
