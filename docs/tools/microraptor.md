# Microraptor

Microraptor is the Dino Tools ingest/parser tool slot. It remains a standalone
repository while Dino Tools records the suite-facing contract.

## Source

- Local repo: `/home/jake/Projects/microraptor`
- Public repo: `https://github.com/Jakeelamb/microraptor`
- Command: `microraptor`
- Dino Tools status: `external`

## Suite Contract

Use Microraptor for sequence-input transport and parsing:

- raw, gzip, and BGZF FASTQ opening
- raw, gzip, and BGZF FASTA opening
- streaming FASTQ and FASTA batch readers
- ordered paired FASTQ validation
- FASTA `.fai` index construction
- reference range fetches and reference chunk streaming

Do not duplicate FASTQ/FASTA parsing in Dino Tools shared crates while this
contract is owned by Microraptor.

## Current Boundary

`dino-core` records Microraptor metadata so `dino microraptor` can explain where
the tool lives and what contract it owns. Dino Tools does not depend on the
Microraptor crate yet.

Promote code or types into Dino Tools only after another suite tool needs the
same stable contract and the shared boundary is small enough for `dino-io` or a
new focused `dino-*` crate.
