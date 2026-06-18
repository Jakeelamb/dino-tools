# dino-quant

TurboQuant-style DNA sketch compression experiments.

```bash
cargo run -p dino-quant --release -- demo
cargo run -p dino-quant --release -- bench-fasta reference.fa.gz
cargo run -p dino-quant --release -- emit-candidates reference.fa reads.fastq
cargo test -p dino-quant
```

`dino-seq` is the only FASTA/FASTQ input backend. Quantized candidates are a
prefilter only; exact alignment stays downstream.

`emit-candidates` writes `target_id` with each row so
`emit-candidate-reference` can reconstruct multi-record references even when
FASTA headers share the same first token. Minimizer reference caches are tied to
the parsed reference content and are rebuilt automatically when that content
changes.
