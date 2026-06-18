# Contributing

Dino Seq is intended for scientific and systems use. Contributions should
preserve correctness, reproducibility, and honest benchmark boundaries.

## Development Checks

Run focused checks while editing:

```bash
cargo fmt --all --check
cargo test -p dino-seq --no-default-features
cargo test -p dino-seq --all-features
cargo clippy -p dino-seq --all-targets --all-features -- -D warnings
cargo package -p dino-seq --allow-dirty --list
```

Run the in-tree smoke benchmark when changing parser hot paths:

```bash
cargo bench -p dino-seq --bench throughput --all-features -- --list
```

## Benchmark Claims

Do not add or strengthen performance claims without benchmark evidence from the
same commit. A useful benchmark contribution records:

- exact command lines and environment variables;
- hardware, OS, Rust toolchain, and feature flags;
- comparator versions;
- raw output from the benchmark tool used;
- checksum, record-count, or base-count parity where relevant.

Machine-local benchmark artifacts are evidence, not universal claims. Keep
parser-only, compression, and workflow-tool comparisons separate.

Use the GitHub benchmark issue template for new public performance wording or
new comparison rows. It asks for the same raw artifacts, command lines, and
claim-boundary checks expected during review.

## Data

Do not commit large biological FASTQ datasets. Use scripts and manifests to
point at local or public datasets, and keep generated large intermediates under
`target/` or another ignored output directory.

## Security And Correctness Reports

Use `SECURITY.md` for vulnerabilities or denial-of-service behavior on
untrusted FASTQ, gzip, or BGZF input. Use the parser bug issue template for
ordinary correctness reports and include a minimal synthetic reproducer instead
of large biological data.
