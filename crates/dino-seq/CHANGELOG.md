# Changelog

All notable changes to dino_seq are documented here. The project is pre-1.0;
public APIs may still change between minor releases.

## Unreleased

### Added

- FASTQ streaming through borrowed batches, one-pass visitors, and bounded chunk
  sinks for caller-owned output buffers.
- FASTA streaming through borrowed batches, visitors, stats helpers, strict
  two-line fast paths, and `.fai`-style indexing helpers.
- Ordered paired-end FASTQ validation for separate R1/R2 streams and adjacent
  interleaved records.
- Optional gzip, BGZF, libdeflate, mmap, SIMD, and pack-path features.
- A focused adoption example in `examples/fastq_chunk_sink.rs`.

### Changed

- Trimmed the crate to keep benchmark and publication artifacts out of the
  package surface.
- Kept only the nightly microbenchmark file as an in-tree benchmark smoke.

### Removed

- Removed checked benchmark snapshots, publication figures, replication docs,
  peer harnesses, and benchmark orchestration scripts from `crates/dino-seq`.
- Removed the `dino-seq-bench` and `dino-seq-fixture` package binaries.
- Removed the public hidden `benchutil` module from the library surface.
