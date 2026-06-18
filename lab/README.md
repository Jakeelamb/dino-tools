# Quant Lab

Tracked here:

- `tools.tsv`: one representative tool per quant application lane.
- `datasets.tsv`: public input data and checksums.
- `matrix.tsv`: first-wave benchmark rows.
- `runs.tsv`: baseline and quant benchmark runs.

Untracked by default:

- modified tool checkouts
- references, reads, indexes, and other input data
- candidate TSVs, benchmark logs, profiles, and run outputs

Default external lab root:

```text
../dino-quant-lab/
```

Use `scripts/quant-lab init` to create the external directories.
Use `scripts/quant-lab acquire all` to fetch the first-wave public inputs.
Use `scripts/quant-lab prepare all` to build first-wave indexes and signatures.
Use `scripts/quant-lab prepare-assisted-minimap2` to cache the current
minimizer-gated masked reference, minimap2 index, and baseline PAF.
Use `scripts/quant-lab run-suite ready` to run missing ready rows once.
Use `scripts/quant-lab run quant_minimap2_assisted_cached_scale_ecoli` to
compare cached full-reference and masked-assisted minimap2 across read counts.
Use `scripts/quant-lab summarize` to print external run metrics.
Use `scripts/quant-lab score-minimap2-candidates CANDIDATES_TSV MINIMAP2_PAF`
to score candidate-window recall against minimap2 mappings.
Use `scripts/quant-lab score-translated-paf BASELINE_PAF ASSISTED_PAF` to
score reduced-reference mapper output with translated coordinates.
Use `scripts/quant-lab optimize-retrieval MODE` to sweep one retrieval mode
across parameter settings and read counts.
Use `scripts/quant-lab new-run TOOL APP VARIANT` to allocate a run directory.

First-wave read rows use deterministic 10k-read slices. Set
`DINO_QUANT_FORCE=1` only when a successful matrix row should be rerun.
Set `DINO_QUANT_OPT_SIZES="100 1000 3000"` to control optimization read counts.
Set `DINO_QUANT_REBUILD_CACHE=1` to rebuild assisted-minimap2 cached artifacts.
Parity rows compare dino-quant candidate windows against baseline mapper output;
they are recall gates, not accelerated mapper implementations.
