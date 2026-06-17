# Gungraun Sidecar

Instruction-count benchmarks for Dino Seq hot paths. This is a standalone Cargo
package so removing `tools/perf/gungraun/` removes the dependency and harness.

Prerequisites:

```bash
valgrind --version
tools/perf/run.sh gungraun-install
```

Run:

```bash
tools/perf/run.sh gungraun
```

Gungraun reports Callgrind-style instruction and cycle estimates. Use it for
small before/after comparisons, not wall-clock throughput. Keep
`benches/throughput.rs` and `tools/perf/run.sh bench` for MiB/s numbers.
