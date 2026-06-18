use crate::retrieval::{
    CandidateWindow, ExperimentalCandidateIndex, RetrievalMode, RetrievalOptions,
};
use dino_quant::{
    QjlResidualSnapshot, QuantizedVector, QuantizedVectorSnapshot, QuantizerConfig,
    ReferenceIndexConfig,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};

const MINIMIZER_CACHE_MAGIC: &[u8; 8] = b"DQMINI4\n";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct ReferenceFingerprint {
    records: usize,
    pub(crate) bases: usize,
    digest: u64,
}

pub(crate) fn reference_fingerprint(
    records: &[dino_quant::SequenceRecord],
) -> ReferenceFingerprint {
    let mut digest = 0xcbf2_9ce4_8422_2325_u64;
    let mut bases = 0_usize;
    for record in records {
        hash_fingerprint_bytes(&mut digest, record.name.as_bytes());
        hash_fingerprint_u64(&mut digest, record.name.len() as u64);
        hash_fingerprint_bytes(&mut digest, &record.bases);
        hash_fingerprint_u64(&mut digest, record.bases.len() as u64);
        bases += record.bases.len();
    }
    hash_fingerprint_u64(&mut digest, records.len() as u64);
    ReferenceFingerprint {
        records: records.len(),
        bases,
        digest,
    }
}

pub(crate) fn load_minimizer_reference_cache(
    path: &str,
    reference_path: &str,
    fingerprint: ReferenceFingerprint,
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
    let cached_fingerprint = read_reference_fingerprint(&mut file)?;
    let cached_config = read_index_config(&mut file)?;
    let cached_retrieval = read_retrieval_options(&mut file)?;
    let reference_bases = read_usize(&mut file)?;
    if cached_reference != reference_path
        || cached_fingerprint != fingerprint
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
            sketch: None,
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

pub(crate) fn save_minimizer_reference_cache(
    path: &str,
    reference_path: &str,
    fingerprint: ReferenceFingerprint,
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
    write_reference_fingerprint(&mut file, fingerprint)?;
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

fn hash_fingerprint_bytes(digest: &mut u64, bytes: &[u8]) {
    for byte in bytes {
        *digest ^= u64::from(*byte);
        *digest = digest.wrapping_mul(0x0000_0100_0000_01b3);
    }
}

fn hash_fingerprint_u64(digest: &mut u64, value: u64) {
    hash_fingerprint_bytes(digest, &value.to_le_bytes());
}

fn write_reference_fingerprint(
    writer: &mut impl Write,
    fingerprint: ReferenceFingerprint,
) -> Result<(), String> {
    write_usize(writer, fingerprint.records)?;
    write_usize(writer, fingerprint.bases)?;
    write_u64(writer, fingerprint.digest)
}

fn read_reference_fingerprint(reader: &mut impl Read) -> Result<ReferenceFingerprint, String> {
    Ok(ReferenceFingerprint {
        records: read_usize(reader)?,
        bases: read_usize(reader)?,
        digest: read_u64(reader)?,
    })
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
