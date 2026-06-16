# Dino Tools

Dino Tools is the umbrella workspace for a Rust-first bioinformatics tool suite.

The suite command is `dino`. Individual tools such as `trex` remain independent
repositories until they are stable enough to promote into the suite.

## Current Shape

- `crates/dino-cli`: umbrella command and discovery surface.
- `crates/dino-core`: shared suite metadata and tool registry types.
- `crates/dino-io`: shared lightweight bioinformatics IO helpers.
- `crates/dino-seq`: FASTQ/FASTA streaming parser and ingest CLI.
- `docs/`: naming, workspace, and promotion rules.
- `CHANGELOG.md`: reviewable history for workspace-level changes.

## Organized Tools

- `dna/`: nucleotide terminal animations (`DNA` CLI). Standalone Rust crate;
  see [dna/README.md](dna/README.md).

- `dino-seq`: workspace ingest/parser for FASTQ/FASTA streaming. See
  [docs/tools/dino-seq.md](docs/tools/dino-seq.md) and
  [crates/dino-seq/README.md](crates/dino-seq/README.md).

## Near-Term Rule

Do not merge existing tools into this workspace by default. Develop each tool in
its own repository, then promote only stable shared contracts into `dino-core` or
`dino-io`.

## Commands

```sh
cargo run -p dino-cli -- list
cargo run -p dino-cli -- status
cargo run -p dino-cli -- dino-seq
cargo run -p dino-seq -- stats --help
cargo test --workspace
```

For the full local readiness gate, see [docs/WORKSPACE.md](docs/WORKSPACE.md).
