use dino_quant::{QuantizedVector, SearchHit};
use std::ops::Deref;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum RetrievalMode {
    Scan,
    Minimizer,
    Simhash,
    Mih,
    Ivf,
    Hnsw,
}

impl RetrievalMode {
    pub(crate) fn parse(value: &str) -> Result<Self, String> {
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

    pub(crate) fn as_str(self) -> &'static str {
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
pub(crate) struct RetrievalOptions {
    pub(crate) mode: RetrievalMode,
    pub(crate) candidate_limit: usize,
    pub(crate) minimizer_k: usize,
    pub(crate) minimizer_window: usize,
    pub(crate) simhash_bands: usize,
    pub(crate) centroids: usize,
    pub(crate) probes: usize,
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
pub(crate) struct IndexedSearchHit {
    pub(crate) target_id: usize,
    pub(crate) hit: SearchHit,
}

impl Deref for IndexedSearchHit {
    type Target = SearchHit;

    fn deref(&self) -> &Self::Target {
        &self.hit
    }
}

#[derive(Clone, Debug)]
pub(crate) struct IndexedSearchResult {
    pub(crate) hits: Vec<IndexedSearchHit>,
    pub(crate) considered: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct CandidateScratch {
    pub(crate) counts: Vec<u32>,
    pub(crate) touched: Vec<usize>,
    pub(crate) ranked: Vec<(usize, u32)>,
    pub(crate) candidate_ids: Vec<usize>,
    pub(crate) minimizer_hashes: Vec<u64>,
    pub(crate) minimizers: Vec<u64>,
    pub(crate) minimizer_deque: Vec<usize>,
}

impl CandidateScratch {
    pub(crate) fn new(window_count: usize) -> Self {
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

    pub(crate) fn reset_counts(&mut self) {
        for idx in self.touched.drain(..) {
            if idx < self.counts.len() {
                self.counts[idx] = 0;
            }
        }
        self.ranked.clear();
        self.candidate_ids.clear();
    }

    pub(crate) fn resize_counts(&mut self, window_count: usize) {
        self.reset_counts();
        self.counts.resize(window_count, 0);
    }

    pub(crate) fn bump(&mut self, idx: usize, weight: u32) -> Result<(), String> {
        let Some(count) = self.counts.get_mut(idx) else {
            return Err("candidate posting id out of bounds".to_owned());
        };
        if *count == 0 {
            self.touched.push(idx);
        }
        *count = count.saturating_add(weight);
        Ok(())
    }

    pub(crate) fn ranked_candidate_ids(&mut self, limit: usize) -> &[usize] {
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
pub(crate) struct CandidateWindow {
    pub(crate) target_id: usize,
    pub(crate) target_start: usize,
    pub(crate) sketch: Option<Vec<f32>>,
    pub(crate) quantized: QuantizedVector,
    pub(crate) simhash64: u64,
    pub(crate) simhash128: [u64; 2],
}

#[derive(Clone, Debug)]
pub(crate) struct ExperimentalCandidateIndex {
    pub(crate) config: dino_quant::ReferenceIndexConfig,
    pub(crate) retrieval: RetrievalOptions,
    pub(crate) target_names: Vec<String>,
    pub(crate) target_offsets: Vec<usize>,
    pub(crate) windows: Vec<CandidateWindow>,
    pub(crate) minimizer_postings: std::collections::HashMap<u64, Vec<usize>>,
    pub(crate) simhash_band_postings: Vec<std::collections::HashMap<u16, Vec<usize>>>,
    pub(crate) mih_postings: Vec<std::collections::HashMap<u16, Vec<usize>>>,
    pub(crate) centroids: Vec<Vec<f32>>,
    pub(crate) centroid_lists: Vec<Vec<usize>>,
    pub(crate) centroid_graph: Vec<Vec<usize>>,
}
