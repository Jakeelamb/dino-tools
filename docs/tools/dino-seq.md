# Dino Seq

Dino Seq is the Dino Tools ingest/parser for FASTQ/FASTA streaming. It lives in
`crates/dino-seq` as a workspace crate.

## Source

- Workspace crate: `crates/dino-seq`
- Library: `dino_seq`
- Command: `dino-seq`

## Suite Contract

Use Dino Seq for sequence-input transport and parsing:

- raw, gzip, and BGZF FASTQ opening
- raw, gzip, and BGZF FASTA opening
- streaming FASTQ and FASTA batch readers
- ordered paired FASTQ validation
- FASTA `.fai` index construction
- reference range fetches and reference chunk streaming

Do not duplicate FASTQ/FASTA parsing in other shared crates while this contract
is owned by Dino Seq.

## Commands

```sh
cargo run -p dino-seq -- stats --help
cargo run -p dino-cli -- dino-seq
```

See `crates/dino-seq/README.md` for the full CLI and library surface.
