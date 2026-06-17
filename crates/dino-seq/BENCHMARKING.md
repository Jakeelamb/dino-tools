# Benchmarking

Dino Seq keeps only a small in-tree benchmark surface. The crate should stay
focused on parser adoption, not publication artifact generation.

## In-Tree Smoke

Run the stable smoke benchmark when changing parser hot paths:

```bash
cargo bench --all-features
```

The benchmark file covers:

- validated FASTQ batch parsing
- no-validation FASTQ batch parsing
- validated FASTQ visitor parsing without batch side tables
- parse-plus-pack side-channel work
- trusted direct parse-plus-pack side-channel work without batch side tables
- BGZF decode-then-parse when the `bgzf` feature is enabled

Run one case when collecting perf counters or flamegraphs:

```bash
cargo bench --bench throughput --all-features -- --list
cargo bench --bench throughput --all-features -- --case pack --records 500000 --iters 20
bench_bin=$(find ../../target/release/deps -maxdepth 1 -type f -executable -name 'throughput-*' | head -n1)
perf stat -d "$bench_bin" --case pack --records 500000 --iters 20
```

Optional local tooling lives in `tools/perf/`. It is a removable sidecar around
the same stable bench target and writes generated output to `target/perf/`:

```bash
tools/perf/run.sh list
DINO_SEQ_PERF_CASE=pack tools/perf/run.sh perf
DINO_SEQ_PERF_CASE=pack tools/perf/run.sh flamegraph
```

The same knobs are available as `DINO_SEQ_BENCH_CASE`,
`DINO_SEQ_BENCH_RECORDS`, `DINO_SEQ_BENCH_READ_LEN`, and
`DINO_SEQ_BENCH_ITERS`.
Use `--case trusted-pack` only to compare the lower-level trusted direct pack
surface against the batch-backed fused pack path on a specific workload.

These timings are local smoke checks. They are not publication evidence and
should not be used for broad performance claims without a separate artifact
pipeline outside this crate.

## Adoption-Shaped Manual Checks

For API-level checks, prefer tiny examples over checked benchmark snapshots:

```bash
cargo run --example fastq_chunk_sink < reads.fastq
cargo run --release --bin dino-seq -- stats reads.fastq
cargo run --release --bin dino-seq -- stats --format fasta reference.fa
```

Chunk consumers should tune `FastqConfig::slab_size`,
`FastqChunkConfig::new(target_bases)`, and
`FastqChunkConfig::min_records` against their own workload. Larger slabs are not
automatically faster for chunked consumers because tail carry and cache pressure
can dominate.

## What Stays Out Of This Crate

Do not check generated benchmark runs, figures, publication notebooks,
replication kits, peer harnesses, or biological datasets into `dino-seq`.
Keep those in an external benchmark/artifact repository or under ignored local
output directories.
