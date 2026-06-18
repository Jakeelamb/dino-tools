# dino_seq

Dino Seq is a small Rust parser crate for streaming FASTQ and FASTA inputs. It
focuses on low-allocation record access for downstream tools rather than
preprocessing workflows.

The default crate is the raw parser core. Compression is a transport layer:
enable `gzip`, `bgzf`, or `transport` only when file-magic convenience matters.

## Scope

Dino Seq provides:

- borrowed FASTQ batches over any `Read`
- one-pass FASTQ visitors for callers that do not need batch side tables
- FASTQ count/stat fast paths that do not build record views
- bounded FASTQ chunk sinks for callers that fill their own output buffers
- streaming multiline FASTA batches and visitors
- ordered paired-end validation for separate R1/R2 or adjacent interleaved reads
- optional BGZF, gzip, mmap, libdeflate, SIMD, and pack-path features

Dino Seq does not trim adapters, filter reads, align reads, produce QC reports,
or synchronize reordered mates.

## FASTQ Entry Points

- `FastqReader::next_batch`: use when the caller needs `RecordRef` side tables.
- `FastqReader::new(std::fs::File::open(path)?)`: raw-file hot path with a
  concrete reader type.
- `FastqReader::count_records`, `count_fastq_read`, and `count_fastq_bytes`:
  use for count/stat workloads that do not need record views.
- `FastqReader::visit_records`: use for whole-stream one-pass consumers.
- `FastqReader::next_chunk_with_sink`: use for aligners or pipelines that need
  bounded chunks in caller-owned structs.
- `visit_fastq_bytes`: use when the complete FASTQ input is already resident.
- `open_fastq` and `open_fastq_with_config`: convenience boxed file openers.
  They open raw files by default and detect gzip/BGZF when those features are
  enabled.

```rust
use dino_seq::FastqReader;

let data = b"@r1\nACGT\n+\nIIII\n";
let mut reader = FastqReader::new(&data[..]);

while let Some(batch) = reader.next_batch()? {
    for record in batch.records() {
        assert_eq!(record.seq(), b"ACGT");
    }
}
# Ok::<(), dino_seq::FastqError>(())
```

Use `FastqReader::next_chunk_with_sink` when downstream code wants to fill
owned buffers without building a batch side table.

## FASTA Entry Points

- `FastaReader::next_batch`: stream multiline FASTA batches.
- `FastaReader::new(std::fs::File::open(path)?)`: raw-file hot path with a
  concrete reader type.
- `FastaReader::visit_records`: visit records without retaining a batch.
- `visit_fasta_bytes`: parse resident multiline FASTA bytes.
- `visit_two_line_fasta_bytes` and `visit_two_line_fasta_read`: strict fast
  paths for canonical `>header` / `sequence` FASTA.
- `open_fasta` and `open_fasta_with_config`: convenience boxed file openers.
  They open raw files by default and detect gzip/BGZF when those features are
  enabled.

```rust
use dino_seq::FastaReader;

let data = b">chr1\nAC\nGT\n";
let mut reader = FastaReader::new(&data[..]);
let batch = reader.next_batch()?.unwrap();
let record = batch.records().next().unwrap();

assert_eq!(record.name_without_gt(), b"chr1");
assert_eq!(record.seq(), b"ACGT");
# Ok::<(), dino_seq::FastqError>(())
```

## Features

- default: raw FASTQ/FASTA parser core only
- `gzip`: ordinary gzip input by gzip magic
- `bgzf`: BGZF reader, writer, detection, indexing, and parallel block helpers
- `transport`: gzip + BGZF convenience transports
- `pure-rust-compression`: compatibility alias for `transport`
- `libdeflate`: optional libdeflate BGZF backends and explicit gzip openers
- `mmap`: resident file visitors and counters backed by read-only memory maps
- `simd`: stable `std::arch` acceleration where supported

## CLI

```bash
cargo run --release --bin dino-seq -- stats reads.fastq
cargo run --release --bin dino-seq -- stats --format fasta reference.fa
cargo run --release --bin dino-seq -- checksum --format fastq reads.fastq
cargo run --release --bin dino-seq -- fasta-index reference.fa
cargo run --release --bin dino-seq -- fasta-fetch reference.fa --fai reference.fa.fai --name chr1 --start 0 --end 100
```

## Checks

```bash
cargo fmt --all --check
cargo test -p dino-seq --no-default-features
cargo test -p dino-seq --all-features
cargo clippy -p dino-seq --all-targets --all-features -- -D warnings
cargo bench --all-features
cargo package --allow-dirty --list
```
