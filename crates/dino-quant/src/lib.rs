use std::f32::consts::PI;
use std::{io::Read, path::Path};

const BASE_A: u8 = 0;
const BASE_C: u8 = 1;
const BASE_G: u8 = 2;
const BASE_T: u8 = 3;
const CONTIG_SEPARATOR: u8 = b'N';

#[derive(Clone, Copy, Debug)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        mix_u64(self.state)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct QuantizerConfig {
    pub bits: u8,
    pub clip_sigma: f32,
    pub rotation_seed: u64,
    pub qjl_seed: u64,
    pub use_qjl_residual: bool,
}

impl Default for QuantizerConfig {
    fn default() -> Self {
        Self {
            bits: 4,
            clip_sigma: 4.0,
            rotation_seed: 0x51d0_51d0_51d0_51d0,
            qjl_seed: 0x71b0_7171_b071_7171,
            use_qjl_residual: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct QjlResidual {
    signs: Vec<u64>,
    dim: usize,
    norm: f32,
    seed: u64,
}

#[derive(Clone, Debug)]
struct PackedCodes {
    data: Vec<u8>,
    len: usize,
    bits: u8,
}

#[derive(Clone, Debug)]
pub struct QuantizedVector {
    dim: usize,
    bits: u8,
    clip: f32,
    rotation_seed: u64,
    codes: PackedCodes,
    qjl_residual: Option<QjlResidual>,
}

#[derive(Clone, Copy, Debug)]
pub struct ReconstructionMetrics {
    pub mse: f32,
    pub cosine: f32,
    pub dot_error: f32,
}

#[derive(Clone, Debug)]
pub struct SequenceRecord {
    pub name: String,
    pub bases: Vec<u8>,
}

#[derive(Clone, Copy, Debug)]
pub struct SequenceRecordRef<'a> {
    pub name: &'a [u8],
    pub bases: &'a [u8],
}

#[derive(Clone, Copy, Debug)]
pub struct ReferenceIndexConfig {
    pub k: usize,
    pub dim: usize,
    pub window_len: usize,
    pub stride: usize,
    pub quantizer: QuantizerConfig,
}

#[derive(Clone, Debug)]
pub struct SearchHit {
    pub start: usize,
    pub end: usize,
    pub score: f32,
}

#[derive(Clone, Debug)]
pub struct ReferenceWindowIndex {
    config: ReferenceIndexConfig,
    windows: Vec<ReferenceWindow>,
}

#[derive(Clone, Debug)]
struct ReferenceWindow {
    start: usize,
    end: usize,
    sketch: QuantizedVector,
}

impl QuantizerConfig {
    pub fn encode(&self, vector: &[f32]) -> Result<QuantizedVector, String> {
        self.validate(vector.len())?;

        let dim = vector.len();
        let clip = self.clip_sigma / (dim as f32).sqrt();
        let rotated = rotate(vector, self.rotation_seed)?;
        let codes = encode_scalar_codes(&rotated, self.bits, clip)?;

        let mut quantized = QuantizedVector {
            dim,
            bits: self.bits,
            clip,
            rotation_seed: self.rotation_seed,
            codes,
            qjl_residual: None,
        };

        if self.use_qjl_residual {
            let base = quantized.decode()?;
            let residual = subtract(vector, &base)?;
            quantized.qjl_residual = QjlResidual::encode(&residual, self.qjl_seed)?;
        }

        Ok(quantized)
    }

    fn validate(&self, dim: usize) -> Result<(), String> {
        if dim == 0 {
            return Err("vector dimension must be non-zero".to_owned());
        }
        if !dim.is_power_of_two() {
            return Err(
                "vector dimension must be a power of two for the Hadamard rotation".to_owned(),
            );
        }
        if !(1..=8).contains(&self.bits) {
            return Err("quantizer bits must be in 1..=8".to_owned());
        }
        if !self.clip_sigma.is_finite() || self.clip_sigma <= 0.0 {
            return Err("clip_sigma must be finite and positive".to_owned());
        }
        Ok(())
    }
}

impl QuantizedVector {
    pub fn decode(&self) -> Result<Vec<f32>, String> {
        if self.codes.len() != self.dim {
            return Err("quantized code length does not match dimension".to_owned());
        }

        let rotated = decode_scalar_codes(&self.codes, self.bits, self.clip)?;
        let mut decoded = inverse_rotate(&rotated, self.rotation_seed)?;
        if let Some(residual) = &self.qjl_residual {
            let correction = residual.decode()?;
            if correction.len() != decoded.len() {
                return Err("QJL correction length does not match decoded vector".to_owned());
            }
            for (dst, corr) in decoded.iter_mut().zip(correction) {
                *dst += corr;
            }
        }
        Ok(decoded)
    }

    pub fn compressed_bits(&self) -> usize {
        let scalar_bits = self.codes.byte_len() * 8;
        let qjl_bits = self
            .qjl_residual
            .as_ref()
            .map_or(0, |residual| residual.dim + 32);
        scalar_bits + qjl_bits
    }

    pub fn compressed_bytes(&self) -> usize {
        self.compressed_bits().div_ceil(8)
    }

    pub fn scalar_dot_rotated_query(&self, rotated_query: &[f32]) -> Result<f32, String> {
        if rotated_query.len() != self.dim {
            return Err("rotated query dimension does not match quantized vector".to_owned());
        }
        let rotated_query_sum = rotated_query.iter().sum();
        self.scalar_dot_rotated_query_with_sum(rotated_query, rotated_query_sum)
    }

    fn scalar_dot_rotated_query_with_sum(
        &self,
        rotated_query: &[f32],
        rotated_query_sum: f32,
    ) -> Result<f32, String> {
        if rotated_query.len() != self.dim {
            return Err("rotated query dimension does not match quantized vector".to_owned());
        }
        self.codes
            .decoded_dot(rotated_query, rotated_query_sum, self.clip)
    }

    fn approximate_dot_4bit_lookup(
        &self,
        lookup: &[[f32; 256]],
        rotated_query_sum: f32,
        rotated_qjl_query: Option<&[f32]>,
    ) -> Result<f32, String> {
        if self.bits != 4 {
            return Err("4-bit lookup scoring requires a 4-bit quantized vector".to_owned());
        }
        let scale = 2.0 * self.clip / 15.0;
        let scalar_score = (-self.clip * rotated_query_sum)
            + scale * self.codes.weighted_code_sum_4_lookup(lookup)?;
        Ok(scalar_score + self.residual_dot_rotated_query(rotated_qjl_query)?)
    }

    fn residual_dot_rotated_query(&self, rotated_qjl_query: Option<&[f32]>) -> Result<f32, String> {
        match (&self.qjl_residual, rotated_qjl_query) {
            (Some(residual), Some(query)) => residual.dot_rotated_query(query),
            (Some(_), None) => Err("QJL residual scoring requires a rotated QJL query".to_owned()),
            (None, _) => Ok(0.0),
        }
    }

    pub fn scalar_dot_query(&self, query: &[f32]) -> Result<f32, String> {
        let rotated_query = rotate(query, self.rotation_seed)?;
        self.scalar_dot_rotated_query(&rotated_query)
    }

    pub fn approximate_dot_query(&self, query: &[f32]) -> Result<f32, String> {
        let rotated_query = rotate(query, self.rotation_seed)?;
        let rotated_query_sum = rotated_query.iter().sum();
        let mut score =
            self.scalar_dot_rotated_query_with_sum(&rotated_query, rotated_query_sum)?;
        if let Some(residual) = &self.qjl_residual {
            let rotated_residual_query = rotate(query, residual.seed)?;
            score += residual.dot_rotated_query(&rotated_residual_query)?;
        }
        Ok(score)
    }
}

impl PackedCodes {
    fn encode(values: &[u8], bits: u8) -> Result<Self, String> {
        if !(1..=8).contains(&bits) {
            return Err("packed code bits must be in 1..=8".to_owned());
        }
        let mask = (1_u16 << bits) - 1;
        let mut data = vec![0_u8; (values.len() * usize::from(bits)).div_ceil(8)];
        for (idx, code) in values.iter().enumerate() {
            if u16::from(*code) > mask {
                return Err("scalar code exceeds bit width".to_owned());
            }
            let bit_pos = idx * usize::from(bits);
            let byte_idx = bit_pos / 8;
            let bit_offset = bit_pos % 8;
            let shifted = u16::from(*code) << bit_offset;
            data[byte_idx] |= shifted as u8;
            if bit_offset + usize::from(bits) > 8 {
                data[byte_idx + 1] |= (shifted >> 8) as u8;
            }
        }

        Ok(Self {
            data,
            len: values.len(),
            bits,
        })
    }

    fn len(&self) -> usize {
        self.len
    }

    fn byte_len(&self) -> usize {
        self.data.len()
    }

    fn get(&self, idx: usize) -> Result<u8, String> {
        if idx >= self.len {
            return Err("packed code index out of bounds".to_owned());
        }
        let bit_pos = idx * usize::from(self.bits);
        let byte_idx = bit_pos / 8;
        let bit_offset = bit_pos % 8;
        let mut value = u16::from(self.data[byte_idx] >> bit_offset);
        if bit_offset + usize::from(self.bits) > 8 {
            value |= u16::from(self.data[byte_idx + 1]) << (8 - bit_offset);
        }
        Ok((value & ((1_u16 << self.bits) - 1)) as u8)
    }

    fn decoded_dot(
        &self,
        rotated_query: &[f32],
        rotated_query_sum: f32,
        clip: f32,
    ) -> Result<f32, String> {
        if rotated_query.len() != self.len {
            return Err("rotated query length does not match packed code length".to_owned());
        }
        if !clip.is_finite() || clip <= 0.0 {
            return Err("clip must be finite and positive".to_owned());
        }

        let levels = (1_u16 << self.bits) - 1;
        let scale = 2.0 * clip / levels as f32;
        let offset_sum = -clip * rotated_query_sum;
        Ok(offset_sum + scale * self.weighted_code_sum(rotated_query)?)
    }

    fn weighted_code_sum(&self, rotated_query: &[f32]) -> Result<f32, String> {
        match self.bits {
            2 => Ok(self.weighted_code_sum_2(rotated_query)),
            4 => Ok(self.weighted_code_sum_4(rotated_query)),
            8 => Ok(self.weighted_code_sum_8(rotated_query)),
            _ => self.weighted_code_sum_generic(rotated_query),
        }
    }

    fn weighted_code_sum_2(&self, rotated_query: &[f32]) -> f32 {
        let mut sum = 0.0_f32;
        for (byte_idx, byte) in self.data.iter().enumerate() {
            let idx = byte_idx * 4;
            if idx >= self.len {
                break;
            }
            sum += f32::from(byte & 0b0000_0011) * rotated_query[idx];
            if idx + 1 < self.len {
                sum += f32::from((byte >> 2) & 0b0000_0011) * rotated_query[idx + 1];
            }
            if idx + 2 < self.len {
                sum += f32::from((byte >> 4) & 0b0000_0011) * rotated_query[idx + 2];
            }
            if idx + 3 < self.len {
                sum += f32::from(byte >> 6) * rotated_query[idx + 3];
            }
        }
        sum
    }

    fn weighted_code_sum_4(&self, rotated_query: &[f32]) -> f32 {
        let mut sum = 0.0_f32;
        for (byte_idx, byte) in self.data.iter().enumerate() {
            let idx = byte_idx * 2;
            if idx >= self.len {
                break;
            }
            sum += f32::from(byte & 0x0f) * rotated_query[idx];
            if idx + 1 < self.len {
                sum += f32::from(byte >> 4) * rotated_query[idx + 1];
            }
        }
        sum
    }

    fn weighted_code_sum_4_lookup(&self, lookup: &[[f32; 256]]) -> Result<f32, String> {
        if lookup.len() != self.data.len() {
            return Err("4-bit lookup table length does not match packed code bytes".to_owned());
        }
        let mut sum = 0.0_f32;
        for (idx, byte) in self.data.iter().enumerate() {
            sum += lookup[idx][usize::from(*byte)];
        }
        Ok(sum)
    }

    fn weighted_code_sum_8(&self, rotated_query: &[f32]) -> f32 {
        let mut sum = 0.0_f32;
        for (code, query_value) in self.data.iter().zip(rotated_query) {
            sum += f32::from(*code) * query_value;
        }
        sum
    }

    fn weighted_code_sum_generic(&self, rotated_query: &[f32]) -> Result<f32, String> {
        let mut sum = 0.0_f32;
        let mask = (1_u32 << self.bits) - 1;
        let mut byte_idx = 0;
        let mut bit_buffer = 0_u32;
        let mut bits_in_buffer = 0_u8;
        for query_value in rotated_query {
            while bits_in_buffer < self.bits {
                if byte_idx >= self.data.len() {
                    return Err("packed code buffer ended early".to_owned());
                }
                bit_buffer |= u32::from(self.data[byte_idx]) << bits_in_buffer;
                bits_in_buffer += 8;
                byte_idx += 1;
            }

            let code = bit_buffer & mask;
            bit_buffer >>= self.bits;
            bits_in_buffer -= self.bits;
            sum += code as f32 * query_value;
        }
        Ok(sum)
    }

    fn decode_all(&self, clip: f32) -> Result<Vec<f32>, String> {
        if !clip.is_finite() || clip <= 0.0 {
            return Err("clip must be finite and positive".to_owned());
        }
        let levels = (1_u16 << self.bits) - 1;
        let scale = 2.0 * clip / levels as f32;
        let offset = -clip;
        let mut decoded = Vec::with_capacity(self.len);
        for idx in 0..self.len {
            decoded.push(offset + f32::from(self.get(idx)?) * scale);
        }
        Ok(decoded)
    }
}

impl ReferenceIndexConfig {
    pub fn validate(&self) -> Result<(), String> {
        if self.window_len == 0 {
            return Err("window_len must be non-zero".to_owned());
        }
        if self.stride == 0 {
            return Err("stride must be non-zero".to_owned());
        }
        if self.k == 0 || self.k > 31 {
            return Err("k must be in 1..=31".to_owned());
        }
        self.quantizer.validate(self.dim)
    }
}

impl ReferenceWindowIndex {
    pub fn build(reference: &[u8], config: ReferenceIndexConfig) -> Result<Self, String> {
        config.validate()?;
        if reference.len() < config.window_len {
            return Err("reference is shorter than the configured window length".to_owned());
        }

        let window_count = ((reference.len() - config.window_len) / config.stride) + 1;
        let mut windows = Vec::with_capacity(window_count);
        for start in (0..=reference.len() - config.window_len).step_by(config.stride) {
            let end = start + config.window_len;
            let sketch = dna_kmer_sketch(&reference[start..end], config.k, config.dim)?;
            let sketch = config.quantizer.encode(&sketch)?;
            windows.push(ReferenceWindow { start, end, sketch });
        }

        Ok(Self { config, windows })
    }

    pub fn search_sequence(&self, query: &[u8], top_k: usize) -> Result<Vec<SearchHit>, String> {
        let sketch = dna_kmer_sketch(query, self.config.k, self.config.dim)?;
        self.search_sketch(&sketch, top_k)
    }

    pub fn search_sketch(&self, query: &[f32], top_k: usize) -> Result<Vec<SearchHit>, String> {
        if top_k == 0 {
            return Ok(Vec::new());
        }
        if query.len() != self.config.dim {
            return Err("query sketch dimension does not match index dimension".to_owned());
        }

        let rotated_query = rotate(query, self.config.quantizer.rotation_seed)?;
        let rotated_query_sum = rotated_query.iter().sum();
        let rotated_qjl_query = if self.config.quantizer.use_qjl_residual {
            Some(rotate(query, self.config.quantizer.qjl_seed)?)
        } else {
            None
        };
        if self.config.quantizer.bits == 4 {
            return self.search_rotated_4bit_lookup(
                &rotated_query,
                rotated_query_sum,
                rotated_qjl_query.as_deref(),
                top_k,
            );
        }

        let mut top = Vec::with_capacity(top_k.min(self.windows.len()));
        for window in &self.windows {
            let scalar_score = window
                .sketch
                .scalar_dot_rotated_query_with_sum(&rotated_query, rotated_query_sum)?;
            let score = scalar_score
                + window
                    .sketch
                    .residual_dot_rotated_query(rotated_qjl_query.as_deref())?;
            push_top_hit(
                &mut top,
                top_k,
                SearchHit {
                    start: window.start,
                    end: window.end,
                    score,
                },
            );
        }
        top.sort_by(|left, right| right.score.total_cmp(&left.score));
        Ok(top)
    }

    fn search_rotated_4bit_lookup(
        &self,
        rotated_query: &[f32],
        rotated_query_sum: f32,
        rotated_qjl_query: Option<&[f32]>,
        top_k: usize,
    ) -> Result<Vec<SearchHit>, String> {
        let lookup = build_4bit_lookup(rotated_query);
        let mut top = Vec::with_capacity(top_k.min(self.windows.len()));
        for window in &self.windows {
            let score = window.sketch.approximate_dot_4bit_lookup(
                &lookup,
                rotated_query_sum,
                rotated_qjl_query,
            )?;
            push_top_hit(
                &mut top,
                top_k,
                SearchHit {
                    start: window.start,
                    end: window.end,
                    score,
                },
            );
        }
        top.sort_by(|left, right| right.score.total_cmp(&left.score));
        Ok(top)
    }

    pub fn window_count(&self) -> usize {
        self.windows.len()
    }

    pub fn compressed_bytes(&self) -> usize {
        self.windows
            .iter()
            .map(|window| window.sketch.compressed_bytes())
            .sum()
    }

    pub fn raw_vector_bytes(&self) -> usize {
        self.windows.len() * self.config.dim * std::mem::size_of::<f32>()
    }

    pub fn compression_ratio(&self) -> f32 {
        if self.compressed_bytes() == 0 {
            return 0.0;
        }
        self.raw_vector_bytes() as f32 / self.compressed_bytes() as f32
    }

    pub fn config(&self) -> ReferenceIndexConfig {
        self.config
    }
}

impl QjlResidual {
    fn encode(residual: &[f32], seed: u64) -> Result<Option<Self>, String> {
        let norm = l2_norm(residual);
        if norm <= f32::EPSILON {
            return Ok(None);
        }

        let projected = rotate(residual, seed)?;
        let mut signs = vec![0_u64; projected.len().div_ceil(64)];
        for (idx, value) in projected.iter().enumerate() {
            if *value >= 0.0 {
                signs[idx / 64] |= 1_u64 << (idx % 64);
            }
        }

        Ok(Some(Self {
            signs,
            dim: projected.len(),
            norm,
            seed,
        }))
    }

    fn decode(&self) -> Result<Vec<f32>, String> {
        if self.dim == 0 || !self.dim.is_power_of_two() {
            return Err("QJL residual dimension must be a non-zero power of two".to_owned());
        }
        let scale = (PI / 2.0).sqrt() * self.norm / self.dim as f32;
        let mut projected = vec![0.0_f32; self.dim];
        for (idx, slot) in projected.iter_mut().enumerate() {
            let bit = (self.signs[idx / 64] >> (idx % 64)) & 1;
            *slot = if bit == 1 { scale } else { -scale };
        }
        inverse_rotate(&projected, self.seed)
    }

    fn dot_rotated_query(&self, rotated_query: &[f32]) -> Result<f32, String> {
        if rotated_query.len() != self.dim {
            return Err("QJL rotated query dimension does not match residual dimension".to_owned());
        }
        let scale = (PI / 2.0).sqrt() * self.norm / self.dim as f32;
        let mut sum = 0.0_f32;
        for (idx, query_value) in rotated_query.iter().enumerate() {
            let bit = (self.signs[idx / 64] >> (idx % 64)) & 1;
            let sign = if bit == 1 { 1.0 } else { -1.0 };
            sum += sign * query_value;
        }
        Ok(scale * sum)
    }
}

pub fn rotate(vector: &[f32], seed: u64) -> Result<Vec<f32>, String> {
    validate_hadamard_dim(vector.len())?;
    let mut rotated = Vec::with_capacity(vector.len());
    for (idx, value) in vector.iter().enumerate() {
        rotated.push(*value * coordinate_sign(seed, idx));
    }
    hadamard_in_place(&mut rotated);
    let scale = 1.0 / (vector.len() as f32).sqrt();
    for value in &mut rotated {
        *value *= scale;
    }
    Ok(rotated)
}

pub fn inverse_rotate(rotated: &[f32], seed: u64) -> Result<Vec<f32>, String> {
    validate_hadamard_dim(rotated.len())?;
    let mut vector = rotated.to_vec();
    hadamard_in_place(&mut vector);
    let scale = 1.0 / (rotated.len() as f32).sqrt();
    for (idx, value) in vector.iter_mut().enumerate() {
        *value *= scale * coordinate_sign(seed, idx);
    }
    Ok(vector)
}

pub fn dna_kmer_sketch(seq: &[u8], k: usize, dim: usize) -> Result<Vec<f32>, String> {
    if k == 0 || k > 31 {
        return Err("k must be in 1..=31".to_owned());
    }
    if dim == 0 {
        return Err("sketch dimension must be non-zero".to_owned());
    }
    if seq.len() < k {
        return Ok(vec![0.0; dim]);
    }

    let mut sketch = vec![0.0_f32; dim];
    for window in seq.windows(k) {
        if let Some(code) = canonical_kmer_code(window) {
            let hash = mix_u64(code ^ ((k as u64) << 56));
            let bucket = (hash as usize) % dim;
            let sign = if (hash >> 63) == 0 { 1.0 } else { -1.0 };
            sketch[bucket] += sign;
        }
    }
    l2_normalize(&mut sketch);
    Ok(sketch)
}

#[cfg(test)]
fn parse_fasta_bytes(input: &[u8]) -> Result<Vec<SequenceRecord>, String> {
    let mut records = Vec::new();
    dino_seq::visit_fasta_bytes(input, |record: dino_seq::FastaVisitRecord<'_>| {
        let mut bases = Vec::with_capacity(record.seq().len());
        append_sequence_line(record.seq(), &mut bases).map_err(dino_seq_format_error)?;
        records.push(SequenceRecord {
            name: parse_record_name(record.name_without_gt()).map_err(dino_seq_format_error)?,
            bases,
        });
        Ok(())
    })
    .map_err(|err| format!("failed to parse FASTA bytes with dino-seq: {err}"))?;
    if records.is_empty() {
        return Err("FASTA input did not contain any records".to_owned());
    }
    Ok(records)
}

pub fn read_fasta_file(path: impl AsRef<Path>) -> Result<Vec<SequenceRecord>, String> {
    let path = path.as_ref();
    let mut reader = dino_seq::open_fasta_for_reference(path).map_err(|err| {
        format!(
            "failed to open FASTA with dino-seq {}: {err}",
            path.display()
        )
    })?;
    let mut records = Vec::new();
    reader
        .visit_records(|record| {
            let mut bases = Vec::with_capacity(record.seq().len());
            append_sequence_line(record.seq(), &mut bases).map_err(dino_seq_format_error)?;
            records.push(SequenceRecord {
                name: parse_record_name(record.name_without_gt()).map_err(dino_seq_format_error)?,
                bases,
            });
            Ok(())
        })
        .map_err(|err| {
            format!(
                "failed to parse FASTA with dino-seq {}: {err}",
                path.display()
            )
        })?;
    if records.is_empty() {
        return Err("FASTA input did not contain any records".to_owned());
    }
    Ok(records)
}

#[cfg(test)]
fn parse_fastq_bytes(input: &[u8]) -> Result<Vec<SequenceRecord>, String> {
    let mut records = Vec::new();
    dino_seq::visit_fastq_bytes(
        input,
        dino_seq::FastqConfig::default(),
        |record: dino_seq::FastqVisitRecord<'_>| {
            let header = record.name();
            let name = header.strip_prefix(b"@").unwrap_or(header);
            let mut bases = Vec::with_capacity(record.seq().len());
            append_sequence_line(record.seq(), &mut bases).map_err(dino_seq_format_error)?;
            records.push(SequenceRecord {
                name: parse_record_name(name).map_err(dino_seq_format_error)?,
                bases,
            });
            Ok(())
        },
    )
    .map_err(|err| format!("failed to parse FASTQ bytes with dino-seq: {err}"))?;
    if records.is_empty() {
        return Err("FASTQ input did not contain any records".to_owned());
    }
    Ok(records)
}

pub fn visit_fastq_slices_file(
    path: impl AsRef<Path>,
    visitor: impl FnMut(SequenceRecordRef<'_>) -> Result<(), String>,
) -> Result<(), String> {
    let path = path.as_ref();
    let mut reader = dino_seq::open_fastq(path).map_err(|err| {
        format!(
            "failed to open FASTQ with dino-seq {}: {err}",
            path.display()
        )
    })?;
    visit_fastq_slices_with_reader(&mut reader, visitor).map_err(|err| {
        format!(
            "failed to parse FASTQ with dino-seq {}: {err}",
            path.display()
        )
    })
}

fn visit_fastq_slices_with_reader<R: Read>(
    reader: &mut dino_seq::FastqReader<R>,
    mut visitor: impl FnMut(SequenceRecordRef<'_>) -> Result<(), String>,
) -> Result<(), String> {
    let chunk_config = dino_seq::FastqChunkConfig::new(4 * 1024 * 1024).min_records(1);
    let mut records = 0_usize;
    while reader
        .next_chunk_with_sink(chunk_config, &mut |record: dino_seq::FastqVisitRecord<
            '_,
        >| {
            let header = record.name();
            let name = header.strip_prefix(b"@").unwrap_or(header);
            let name = parse_record_name_bytes(name).map_err(dino_seq_format_error)?;
            validate_sequence_bases(record.seq()).map_err(dino_seq_format_error)?;
            visitor(SequenceRecordRef {
                name,
                bases: record.seq(),
            })
            .map_err(dino_seq_format_error)?;
            records += 1;
            Ok(())
        })
        .map_err(|err| err.to_string())?
        .is_some()
    {}
    if records == 0 {
        return Err("FASTQ input did not contain any records".to_owned());
    }
    Ok(())
}

pub fn concatenate_records(records: &[SequenceRecord]) -> Result<Vec<u8>, String> {
    if records.is_empty() {
        return Err("at least one sequence record is required".to_owned());
    }
    let total_bases = records
        .iter()
        .map(|record| record.bases.len())
        .sum::<usize>();
    let mut concatenated = Vec::with_capacity(total_bases + records.len().saturating_sub(1));
    for (idx, record) in records.iter().enumerate() {
        if idx > 0 {
            concatenated.push(CONTIG_SEPARATOR);
        }
        concatenated.extend_from_slice(&record.bases);
    }
    Ok(concatenated)
}

pub fn reconstruction_metrics(
    original: &[f32],
    decoded: &[f32],
) -> Result<ReconstructionMetrics, String> {
    Ok(ReconstructionMetrics {
        mse: mse(original, decoded)?,
        cosine: cosine_similarity(original, decoded)?,
        dot_error: (dot(original, original)? - dot(original, decoded)?).abs(),
    })
}

pub fn mutate_dna(seq: &[u8], every: usize) -> Vec<u8> {
    if every == 0 {
        return seq.to_vec();
    }
    let mut mutated = seq.to_vec();
    for idx in (every - 1..mutated.len()).step_by(every) {
        mutated[idx] = match mutated[idx].to_ascii_uppercase() {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            b'T' => b'A',
            other => other,
        };
    }
    mutated
}

pub fn synthetic_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = SplitMix64::new(seed);
    let mut seq = Vec::with_capacity(len);
    for _ in 0..len {
        let base = match rng.next_u64() & 3 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        seq.push(base);
    }
    seq
}

pub fn intervals_overlap(
    left_start: usize,
    left_end: usize,
    right_start: usize,
    right_end: usize,
) -> bool {
    left_start < right_end && right_start < left_end
}

pub fn dot(left: &[f32], right: &[f32]) -> Result<f32, String> {
    if left.len() != right.len() {
        return Err("vectors must have equal length".to_owned());
    }
    Ok(left.iter().zip(right).map(|(a, b)| a * b).sum())
}

pub fn mse(left: &[f32], right: &[f32]) -> Result<f32, String> {
    if left.len() != right.len() {
        return Err("vectors must have equal length".to_owned());
    }
    if left.is_empty() {
        return Err("vectors must be non-empty".to_owned());
    }
    let sum: f32 = left
        .iter()
        .zip(right)
        .map(|(a, b)| {
            let delta = a - b;
            delta * delta
        })
        .sum();
    Ok(sum / left.len() as f32)
}

pub fn cosine_similarity(left: &[f32], right: &[f32]) -> Result<f32, String> {
    let denom = l2_norm(left) * l2_norm(right);
    if denom <= f32::EPSILON {
        return Err("cosine similarity is undefined for zero vectors".to_owned());
    }
    Ok(dot(left, right)? / denom)
}

pub fn l2_normalize(vector: &mut [f32]) {
    let norm = l2_norm(vector);
    if norm <= f32::EPSILON {
        return;
    }
    for value in vector {
        *value /= norm;
    }
}

fn encode_scalar_codes(rotated: &[f32], bits: u8, clip: f32) -> Result<PackedCodes, String> {
    let levels = (1_u16 << bits) - 1;
    let inv_width = 1.0 / (2.0 * clip);
    let mut codes = Vec::with_capacity(rotated.len());
    for value in rotated {
        let clipped = value.clamp(-clip, clip);
        let normalized = (clipped + clip) * inv_width;
        codes.push((normalized * levels as f32).round() as u8);
    }
    PackedCodes::encode(&codes, bits)
}

fn decode_scalar_codes(codes: &PackedCodes, bits: u8, clip: f32) -> Result<Vec<f32>, String> {
    if !(1..=8).contains(&bits) {
        return Err("quantizer bits must be in 1..=8".to_owned());
    }
    if !clip.is_finite() || clip <= 0.0 {
        return Err("clip must be finite and positive".to_owned());
    }
    if codes.bits != bits {
        return Err("packed code bit width does not match quantizer bit width".to_owned());
    }
    codes.decode_all(clip)
}

fn subtract(left: &[f32], right: &[f32]) -> Result<Vec<f32>, String> {
    if left.len() != right.len() {
        return Err("vectors must have equal length".to_owned());
    }
    Ok(left.iter().zip(right).map(|(a, b)| a - b).collect())
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

fn build_4bit_lookup(rotated_query: &[f32]) -> Vec<[f32; 256]> {
    let mut lookup = Vec::with_capacity(rotated_query.len().div_ceil(2));
    for chunk in rotated_query.chunks(2) {
        let first = chunk[0];
        let second = if chunk.len() == 2 { chunk[1] } else { 0.0 };
        let mut table = [0.0_f32; 256];
        for byte in 0_u16..=255 {
            let low = f32::from((byte & 0x0f) as u8);
            let high = f32::from((byte >> 4) as u8);
            table[usize::from(byte)] = low * first + high * second;
        }
        lookup.push(table);
    }
    lookup
}

fn trim_ascii(mut bytes: &[u8]) -> &[u8] {
    while matches!(bytes.first(), Some(b' ' | b'\t' | b'\r' | b'\n')) {
        bytes = &bytes[1..];
    }
    while matches!(bytes.last(), Some(b' ' | b'\t' | b'\r' | b'\n')) {
        bytes = &bytes[..bytes.len() - 1];
    }
    bytes
}

fn parse_record_name(header: &[u8]) -> Result<String, String> {
    let name = parse_record_name_bytes(header)?;
    String::from_utf8(name.to_vec()).map_err(|_| "sequence record name must be UTF-8".to_owned())
}

fn parse_record_name_bytes(header: &[u8]) -> Result<&[u8], String> {
    let trimmed = trim_ascii(header);
    if trimmed.is_empty() {
        return Err("sequence record header must contain a name".to_owned());
    }
    trimmed
        .split(|byte| matches!(*byte, b' ' | b'\t'))
        .next()
        .filter(|name| !name.is_empty())
        .ok_or_else(|| "sequence record header must contain a name".to_owned())
}

fn append_sequence_line(line: &[u8], dst: &mut Vec<u8>) -> Result<(), String> {
    for base in line {
        match base.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' | b'N' => dst.push(base.to_ascii_uppercase()),
            b' ' | b'\t' | b'\r' => {}
            other => {
                return Err(format!(
                    "unsupported sequence character '{}' in input",
                    char::from(other)
                ));
            }
        }
    }
    Ok(())
}

fn validate_sequence_bases(line: &[u8]) -> Result<(), String> {
    for base in line {
        match base.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' | b'N' => {}
            other => {
                return Err(format!(
                    "unsupported sequence character '{}' in input",
                    char::from(other)
                ));
            }
        }
    }
    Ok(())
}

fn dino_seq_format_error(message: String) -> dino_seq::FastqError {
    dino_seq::FastqError::Format(message)
}

fn l2_norm(vector: &[f32]) -> f32 {
    vector.iter().map(|value| value * value).sum::<f32>().sqrt()
}

fn validate_hadamard_dim(dim: usize) -> Result<(), String> {
    if dim == 0 {
        return Err("Hadamard dimension must be non-zero".to_owned());
    }
    if !dim.is_power_of_two() {
        return Err("Hadamard dimension must be a power of two".to_owned());
    }
    Ok(())
}

fn hadamard_in_place(values: &mut [f32]) {
    let mut stride = 1;
    while stride < values.len() {
        let step = stride * 2;
        for start in (0..values.len()).step_by(step) {
            for idx in start..start + stride {
                let left = values[idx];
                let right = values[idx + stride];
                values[idx] = left + right;
                values[idx + stride] = left - right;
            }
        }
        stride = step;
    }
}

fn coordinate_sign(seed: u64, idx: usize) -> f32 {
    if mix_u64(seed ^ idx as u64) & 1 == 0 {
        1.0
    } else {
        -1.0
    }
}

fn canonical_kmer_code(window: &[u8]) -> Option<u64> {
    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    for (idx, base) in window.iter().enumerate() {
        let code = base_code(*base)?;
        forward = (forward << 2) | u64::from(code);
        let rc = u64::from(BASE_T - code);
        reverse |= rc << (idx * 2);
    }
    Some(forward.min(reverse))
}

fn base_code(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(BASE_A),
        b'C' => Some(BASE_C),
        b'G' => Some(BASE_G),
        b'T' => Some(BASE_T),
        _ => None,
    }
}

fn mix_u64(mut value: u64) -> u64 {
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotation_round_trips() -> Result<(), String> {
        let mut vector = (0..128)
            .map(|idx| ((idx as f32 + 1.0) * 0.17).sin())
            .collect::<Vec<_>>();
        l2_normalize(&mut vector);

        let rotated = rotate(&vector, 17)?;
        let decoded = inverse_rotate(&rotated, 17)?;
        let err = mse(&vector, &decoded)?;
        assert!(err < 1.0e-12, "round-trip MSE was {err}");
        Ok(())
    }

    #[test]
    fn more_scalar_bits_reduce_error() -> Result<(), String> {
        let mut vector = (0..256)
            .map(|idx| ((idx as f32 + 3.0) * 0.11).cos())
            .collect::<Vec<_>>();
        l2_normalize(&mut vector);

        let low = QuantizerConfig {
            bits: 2,
            use_qjl_residual: false,
            ..QuantizerConfig::default()
        }
        .encode(&vector)?
        .decode()?;
        let high = QuantizerConfig {
            bits: 5,
            use_qjl_residual: false,
            ..QuantizerConfig::default()
        }
        .encode(&vector)?
        .decode()?;

        assert!(mse(&vector, &high)? < mse(&vector, &low)?);
        Ok(())
    }

    #[test]
    fn packed_codes_round_trip_values() -> Result<(), String> {
        let values = (0..97).map(|idx| (idx % 7) as u8).collect::<Vec<_>>();
        let packed = PackedCodes::encode(&values, 3)?;

        assert!(packed.byte_len() < values.len());
        for (idx, expected) in values.iter().enumerate() {
            assert_eq!(packed.get(idx)?, *expected);
        }
        Ok(())
    }

    #[test]
    fn dna_sketches_are_normalized() -> Result<(), String> {
        let seq = synthetic_dna(512, 9);
        let sketch = dna_kmer_sketch(&seq, 15, 128)?;
        let norm = l2_norm(&sketch);
        assert!((norm - 1.0).abs() < 1.0e-5, "norm was {norm}");
        Ok(())
    }

    #[test]
    fn parses_fasta_records_and_concatenates_with_separator() -> Result<(), String> {
        let records = parse_fasta_bytes(b">chr1 description\nacgt\nNN\n>chr2\nTTA\n")?;
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "chr1");
        assert_eq!(records[0].bases, b"ACGTNN");
        assert_eq!(records[1].name, "chr2");
        assert_eq!(concatenate_records(&records)?, b"ACGTNNNTTA");
        Ok(())
    }

    #[test]
    fn parses_fastq_records() -> Result<(), String> {
        let records = parse_fastq_bytes(b"@read1 comment\nacgtn\n+\nIIIII\n@read2\nTTA\n+\n###\n")?;
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "read1");
        assert_eq!(records[0].bases, b"ACGTN");
        assert_eq!(records[1].name, "read2");
        assert_eq!(records[1].bases, b"TTA");
        Ok(())
    }

    #[test]
    fn rejects_invalid_sequence_character() {
        let err = parse_fasta_bytes(b">chr1\nACGTX\n").expect_err("invalid base should fail");
        assert!(err.contains("unsupported sequence character"));
    }

    #[test]
    fn qjl_residual_decodes_finite_vector() -> Result<(), String> {
        let seq = synthetic_dna(2048, 42);
        let vector = dna_kmer_sketch(&seq, 17, 256)?;
        let quantized = QuantizerConfig {
            bits: 3,
            use_qjl_residual: true,
            ..QuantizerConfig::default()
        }
        .encode(&vector)?;
        let decoded = quantized.decode()?;

        assert_eq!(decoded.len(), vector.len());
        assert!(decoded.iter().all(|value| value.is_finite()));
        assert!(quantized.compressed_bits() < vector.len() * 32);
        Ok(())
    }

    #[test]
    fn approximate_dot_matches_decoded_qjl_dot() -> Result<(), String> {
        let reference = dna_kmer_sketch(&synthetic_dna(2048, 42), 17, 256)?;
        let query = dna_kmer_sketch(&mutate_dna(&synthetic_dna(2048, 77), 19), 17, 256)?;
        let quantized = QuantizerConfig {
            bits: 4,
            use_qjl_residual: true,
            ..QuantizerConfig::default()
        }
        .encode(&reference)?;

        let decoded = quantized.decode()?;
        let decoded_dot = dot(&decoded, &query)?;
        let approximate_dot = quantized.approximate_dot_query(&query)?;

        assert!(
            (decoded_dot - approximate_dot).abs() < 1.0e-5,
            "decoded_dot={decoded_dot} approximate_dot={approximate_dot}"
        );
        Ok(())
    }

    #[test]
    fn compressed_index_finds_mutated_source_window() -> Result<(), String> {
        let reference = synthetic_dna(8192, 0xabc);
        let config = ReferenceIndexConfig {
            k: 15,
            dim: 256,
            window_len: 512,
            stride: 128,
            quantizer: QuantizerConfig {
                bits: 4,
                use_qjl_residual: false,
                ..QuantizerConfig::default()
            },
        };
        let index = ReferenceWindowIndex::build(&reference, config)?;
        let source_start = 2816;
        let source_end = source_start + config.window_len;
        let query = mutate_dna(&reference[source_start..source_end], 41);
        let hits = index.search_sequence(&query, 5)?;

        assert!(
            hits.iter()
                .any(|hit| { intervals_overlap(hit.start, hit.end, source_start, source_end) })
        );
        assert!(index.compression_ratio() > 6.0);
        Ok(())
    }

    #[test]
    fn canonical_kmers_match_reverse_complements() -> Result<(), String> {
        let left = canonical_kmer_code(b"ACGTT").ok_or("left k-mer should be valid")?;
        let right = canonical_kmer_code(b"AACGT").ok_or("right k-mer should be valid")?;
        assert_eq!(left, right);
        Ok(())
    }
}
