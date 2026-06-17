# dino_seq

Dino Seq is a small Rust parser crate for streaming FASTQ and FASTA inputs. It
focuses on low-allocation record access for downstream tools rather than
preprocessing workflows.

The core invariant is simple: decompression is a transport layer. Raw, gzip, and
BGZF inputs feed the same format parser after input detection.

## Scope

Dino Seq provides:

- borrowed FASTQ batches over any `Read`
- one-pass FASTQ visitors for callers that do not need batch side tables
- bounded FASTQ chunk sinks for callers that fill their own output buffers
- streaming multiline FASTA batches and visitors
- ordered paired-end validation for separate R1/R2 or adjacent interleaved reads
- optional BGZF, gzip, mmap, libdeflate, SIMD, and pack-path features

Dino Seq does not trim adapters, filter reads, align reads, produce QC reports,
or synchronize reordered mates.

## FASTQ Entry Points

- `FastqReader::next_batch`: use when the caller needs `RecordRef` side tables.
- `FastqReader::visit_records`: use for whole-stream one-pass consumers.
- `FastqReader::next_chunk_with_sink`: use for aligners or pipelines that need
  bounded chunks in caller-owned structs.
- `visit_fastq_bytes`: use when the complete FASTQ input is already resident.
- `open_fastq` and `open_fastq_with_config`: use for raw, gzip, or BGZF files.

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

For owned downstream buffers, see
[`examples/fastq_chunk_sink.rs`](examples/fastq_chunk_sink.rs).

## FASTA Entry Points

- `FastaReader::next_batch`: stream multiline FASTA batches.
- `FastaReader::visit_records`: visit records without retaining a batch.
- `visit_fasta_bytes`: parse resident multiline FASTA bytes.
- `visit_two_line_fasta_bytes` and `visit_two_line_fasta_read`: strict fast
  paths for canonical `>header` / `sequence` FASTA.
- `open_fasta` and `open_fasta_with_config`: use for raw, gzip, or BGZF files.

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

- default: raw, gzip, and BGZF streaming
- `gzip`: ordinary gzip input by gzip magic
- `bgzf`: BGZF reader, writer, detection, indexing, and parallel block helpers
- `libdeflate`: optional libdeflate BGZF backends and explicit gzip openers
- `mmap`: resident file visitors backed by read-only memory maps
- `simd`: stable `std::arch` acceleration where supported
- `asm-scan`: x86-64 newline scanner experiment

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
cargo test --all
cargo clippy --all-targets --all-features -- -D warnings
cargo bench --all-features
cargo package --allow-dirty --list
```

## Documentation

- [`docs/API_SURFACE.md`](docs/API_SURFACE.md): intended public API tiers.
- [`BENCHMARKING.md`](BENCHMARKING.md): small in-tree benchmark guidance.
