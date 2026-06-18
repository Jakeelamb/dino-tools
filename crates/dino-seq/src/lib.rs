#![warn(missing_docs)]
//! Low-allocation FASTQ and FASTA parsing.
//!
//! Dino Seq is a library core for downstream scientific tools that need
//! low-allocation streams of ordinary FASTQ records and multiline FASTA
//! records. The default crate is the raw parser core; gzip and BGZF transports
//! are explicit feature opt-ins.
//!
//! The default feature set builds on stable Rust with no compression
//! dependencies. SIMD acceleration is available through the explicit `simd`
//! feature on stable Rust targets that expose supported `std::arch` intrinsics.
//!
//! # Choosing an entry point
//!
//! Use [`FastqReader`] or [`FastaReader`] directly over a concrete
//! [`std::fs::File`] for raw-file hot paths. Use [`count_fastq_read`] for
//! count/stat workloads that do not need record views. Use [`visit_fastq_bytes`]
//! when a complete FASTQ byte buffer is already resident in memory. Use
//! [`open_fastq`], [`open_fastq_with_config`], [`open_fasta`], or
//! [`open_fasta_with_config`] when file-path convenience and optional
//! gzip/BGZF magic detection are worth the boxed transport. Use
//! [`PairedFastqReader`] or [`open_paired_fastq`] for ordered R1/R2 streams. Use
//! [`visit_fasta_bytes`] for resident multiline FASTA, and
//! [`visit_two_line_fasta_bytes`] or [`visit_two_line_fasta_read`] for strict
//! canonical two-line FASTA fast paths.
//!
//! # Scope
//!
//! Dino Seq parses four-line FASTQ records. It validates ordered paired-end
//! reads in separate R1/R2 streams or adjacent interleaved records, but it does
//! not synchronize reordered mates. It does not trim adapters, filter reads,
//! align reads, or generate quality-control reports.
//!
//! # Lifetimes and allocation
//!
//! Batches borrow from the reader's reusable storage. A [`FastqBatch`] or
//! [`FastaBatch`] is valid until the next mutable reader call. Clone or copy
//! record data if it must outlive the batch. This design keeps the parser
//! low-allocation, but it means callers should process each batch before
//! advancing the reader.
//!
//! # Feature flags
//!
//! - `gzip` enables ordinary gzip auto-detection and streaming decode.
//! - `bgzf` enables BGZF readers, writers, indexing, and adaptive parallel
//!   decoding.
//! - `libdeflate` enables optional libdeflate BGZF backends and an explicit
//!   buffered gzip opener.
//! - `mmap` enables resident file visitors backed by memory maps.
//! - `transport` enables both `gzip` and `bgzf` for callers that want the old
//!   all-transport convenience build.
//! - `pure-rust-compression` is a compatibility alias for `transport`.
//! - `simd` enables stable `std::arch` scanner and packing paths where
//!   supported.
//! - `asm-scan` (x86-64 only) uses a hand-written AVX2 newline scanner instead of
//!   LLVM-generated intrinsics for the internal newline scan, for experiments.
//!   When both `simd` and `asm-scan` are enabled on x86-64, `asm-scan` takes
//!   precedence for newline scanning only.
//!
//! # Example
//!
//! ```
//! use dino_seq::FastqReader;
//!
//! let data = b"@r1\nACGT\n+\nIIII\n";
//! let mut reader = FastqReader::new(&data[..]);
//! let mut records = 0;
//!
//! while let Some(batch) = reader.next_batch()? {
//!     for record in batch.records() {
//!         assert_eq!(record.seq(), b"ACGT");
//!         records += 1;
//!     }
//! }
//!
//! assert_eq!(records, 1);
//! # Ok::<(), dino_seq::FastqError>(())
//! ```
//!
//! # Paired reads
//!
//! ```
//! use dino_seq::PairedFastqReader;
//!
//! let r1 = b"@frag/1\nACGT\n+\nIIII\n";
//! let r2 = b"@frag/2\nTGCA\n+\nJJJJ\n";
//! let mut reader = PairedFastqReader::new(&r1[..], &r2[..]);
//! let batch = reader.next_pair_batch()?.expect("one paired batch");
//! let pair = batch.pairs().next().expect("one read pair");
//!
//! assert_eq!(pair.pair_id(), b"frag");
//! assert_eq!(pair.first().seq(), b"ACGT");
//! assert_eq!(pair.second().seq(), b"TGCA");
//! # Ok::<(), dino_seq::FastqError>(())
//! ```

#[cfg(feature = "bgzf")]
mod bgzf;
mod error;
mod fasta;
mod fastq;
mod fastq_frame;
#[cfg(feature = "mmap")]
mod mmap;
/// Base/quality packing and trusted four-line FASTQ pack paths.
///
/// The high-level FASTQ readers expose borrowed records. This module provides
/// the lower-level side-channel representation used when downstream code wants
/// packed two-bit bases, ambiguity masks, and Phred+33 summaries without
/// allocating an owned record per read.
pub mod pack;
mod scan;
mod source;

#[cfg(feature = "bgzf")]
pub use bgzf::{
    BGZF_EOF_BLOCK, BgzfAutoReader, BgzfDecodedBlock, BgzfDecodedBlockReader, BgzfDeflateBackend,
    BgzfIndex, BgzfIndexEntry, BgzfInflateBackend, BgzfParallelConfig, BgzfParallelReader,
    BgzfPipelineMetrics, BgzfPipelineMetricsSnapshot, BgzfReader, BgzfSeekReader,
    BgzfVirtualOffset, BgzfWriter, build_bgzf_index, build_bgzf_index_strict,
    compress_bgzf_parallel, compress_bgzf_parallel_with_deflate_backend, decompress_bgzf_parallel,
    decompress_bgzf_parallel_with_inflate_backend,
};
pub use error::{FastqError, FastqPosition, Result};
#[cfg(feature = "bgzf")]
pub use fasta::{BgzfFastaReferenceChunks, BgzfIndexedFastaReader, build_fasta_index_bgzf};
pub use fasta::{
    FastaBatch, FastaConfig, FastaIndex, FastaIndexEntry, FastaPartition, FastaPartitionConfig,
    FastaReader, FastaRecord, FastaRecordRef, FastaRecordSink, FastaReferenceChunk,
    FastaReferenceChunkRef, FastaReferenceChunkSink, FastaReferenceChunks, FastaShape, FastaStats,
    FastaVisitRecord, IndexedFastaReader, OwnedFastaBatch, OwnedFastaRecord, build_fasta_index,
    count_fasta_bytes, count_fasta_read, count_two_line_fasta_bytes, count_two_line_fasta_read,
    detect_fasta_shape, plan_fasta_partitions, visit_fasta_bytes, visit_fasta_bytes_auto,
    visit_two_line_fasta_bytes, visit_two_line_fasta_read,
};
pub use fastq::{
    FastqBatch, FastqChunkConfig, FastqChunkSinkExt, FastqChunkStats, FastqConfig, FastqPair,
    FastqReader, FastqRecord, FastqRecordSink, FastqStats, FastqVisitRecord, InterleavedPairs,
    PairValidation, PairedFastqBatch, PairedFastqPairs, PairedFastqReader, PairedRecords,
    PairingMode, RecordRef, count_fastq_bytes, count_fastq_read, count_fastq_read_with_config,
    paired_records, strip_pair_suffix, visit_fastq_bytes,
};
#[cfg(feature = "mmap")]
pub use mmap::{count_fasta_mmap, count_fastq_mmap, visit_fasta_mmap, visit_fastq_mmap};
pub use source::{
    DetectedInputKind, detect_file_input_kind, open_fasta, open_fasta_for_reference,
    open_fasta_with_config, open_fastq, open_fastq_with_config, open_paired_fastq,
    open_paired_fastq_with_config, open_paired_fastq_with_configs,
};
#[cfg(all(feature = "gzip", feature = "libdeflate"))]
pub use source::{
    LibdeflateGzipLimits, open_fasta_gzip_libdeflate, open_fasta_gzip_libdeflate_with_config,
    open_fasta_gzip_libdeflate_with_limits, open_fastq_gzip_libdeflate,
    open_fastq_gzip_libdeflate_with_config, open_fastq_gzip_libdeflate_with_limits,
};
#[cfg(feature = "bgzf")]
pub use source::{
    open_fastq_bgzf_adaptive, open_fastq_bgzf_flate2, open_fastq_bgzf_parallel,
    open_fastq_bgzf_parallel_with_backend, open_fastq_bgzf_parallel_with_config,
    open_fastq_bgzf_parallel_with_options, open_fastq_bgzf_with_backend,
};
