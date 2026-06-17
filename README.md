# Dino Tools

Dino Tools is a small workspace for standalone Rust tools.

There is no umbrella `dino` command. Each tool owns its own package, binary,
README, tests, and release surface. The root workspace exists only so local
development can build and test the checked-in tools together.

## Current Shape

- `crates/dino-seq`: FASTQ/FASTA streaming parser and ingest CLI.
- `crates/dino-quant`: TurboQuant-inspired DNA sketch compression experiments.
- `dna/`: nucleotide terminal animations (`DNA` CLI).
- `CHANGELOG.md`: reviewable history for workspace-level changes only.

## Near-Term Rule

Keep tools isolated unless shared code proves it is worth extracting. Prefer a
little duplication over a framework layer that couples unrelated tools.

## Install

```sh
cargo install --path crates/dino-seq
cargo install --path dna
```

## Commands

```sh
cargo run -p dino-seq -- stats --help
cargo run -p dino-quant -- demo
cargo run -p dna --bin DNA -- --help
cargo test --workspace
```
