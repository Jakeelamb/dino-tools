# Workspace Boundaries

This workspace is the suite layer, not a dumping ground for every bioinformatics
experiment.

## Promote Into Dino Tools When

- The behavior is stable enough to document.
- Two or more tools need the same contract.
- The code has tests and does not depend on a tool-specific workflow.
- The crate boundary can stay small.

## Keep External When

- The tool is still changing rapidly.
- The code is deeply tied to one executable.
- The API is not yet proven by a second caller.
- Moving it would slow down focused development.

## Initial Promotion Order

1. Shared names and registry metadata in `dino-core`.
2. Streaming file format helpers in `dino-io`.
3. Stable shared data structures, one at a time.
4. Tool dispatch only after the standalone tool interface is stable.

## External Tool Slots

External tools can be organized here before they are promoted. Start by adding
their source location, command name, and suite-facing contracts to `dino-core`,
then add a focused page under `docs/tools/`.

This is the current pattern for Dino Seq: the FASTQ/FASTA ingest implementation
lives in `crates/dino-seq`, and `dino dino-seq` reports workspace metadata.

## Readiness Gates

Use these checks before committing workspace changes:

```sh
cargo fmt --all --check
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo run -p dino-cli -- list
cargo run -p dino-cli -- dino-seq
```

`dino-core` and `dino-io` should also pass `cargo package --allow-dirty` when
their public metadata changes. `dino-cli` is an unpublished umbrella binary and
depends on local workspace crates, so its readiness gate is the workspace build,
test, clippy, and command smoke test rather than crates.io package verification.
