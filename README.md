# Dino Tools

Dino Tools is the umbrella workspace for a Rust-first bioinformatics tool suite.

The suite command is `dino`. Individual tools such as `trex` and `microraptor`
remain independent repositories until they are stable enough to promote into the
suite.

## Current Shape

- `crates/dino-cli`: umbrella command and discovery surface.
- `crates/dino-core`: shared suite metadata and tool registry types.
- `crates/dino-io`: shared lightweight bioinformatics IO helpers.
- `docs/`: naming, workspace, and promotion rules.
- `CHANGELOG.md`: reviewable history for workspace-level changes.

## Organized Tools

- `dna/`: nucleotide terminal animations (`DNA` CLI). Standalone Rust crate;
  see [dna/README.md](dna/README.md).

- `microraptor`: external ingest/parser tool for FASTQ/FASTA streaming. See
  [docs/tools/microraptor.md](docs/tools/microraptor.md).

## Near-Term Rule

Do not merge existing tools into this workspace by default. Develop each tool in
its own repository, then promote only stable shared contracts into `dino-core` or
`dino-io`.

## Commands

```sh
cargo run -p dino-cli -- list
cargo run -p dino-cli -- status
cargo run -p dino-cli -- microraptor
cargo test --workspace
```

For the full local readiness gate, see [docs/WORKSPACE.md](docs/WORKSPACE.md).
