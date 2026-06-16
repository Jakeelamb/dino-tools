# Naming

## Suite

- Suite name: `Dino Tools`
- Umbrella command: `dino`
- Repository: `dino-tools`

## Crates

Use the `dino-*` prefix for suite-owned shared crates:

- `dino-core`
- `dino-io`
- `dino-kmers`
- `dino-graph`
- `dino-align`
- `dino-bench`

## Tools

Use dinosaur names for user-facing tools. Keep each tool name tied to a durable
job instead of a temporary implementation detail.

| Tool | Intended Role | Current Status |
| --- | --- | --- |
| `trex` | assembler / heavy de novo engine | external |
| `microraptor` | FASTQ/FASTA parsing and ingest | external |
| `velociraptor` | fast search, sketching, and prefiltering | planned |
| `ankylosaurus` | QC, validation, and defensive filtering | planned |
| `stegosaurus` | graph layout and scaffolding | planned |
| `triceratops` | comparison, reconciliation, and consensus | planned |
| `brachiosaurus` | long-read and large-reference workflows | planned |
| `archaeopteryx` | import/export and compatibility bridges | planned |

## Command Pattern

Prefer:

```sh
dino list
dino status
dino trex assemble ...
dino microraptor stats ...
```

The last two forms are future promotion targets. Until a tool is promoted, keep
the standalone repository command as the source of truth.

Use focused docs under `docs/tools/` for external tools that have a clear
suite-facing role before they are promoted into workspace crates.
