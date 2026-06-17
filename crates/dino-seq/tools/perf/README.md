# Perf Tools

Optional local profiling wrapper for Dino Seq hot paths. This folder is a
removable sidecar: deleting `tools/perf/` removes the extra tooling surface.
Generated output goes under `target/perf/`, which stays out of package artifacts.

```bash
tools/perf/run.sh list
tools/perf/run.sh bench
tools/perf/run.sh perf
tools/perf/run.sh hyperfine
tools/perf/run.sh flamegraph
tools/perf/run.sh gungraun-install
tools/perf/run.sh gungraun
```

Knobs:

```bash
DINO_SEQ_PERF_CASE=pack
DINO_SEQ_PERF_RECORDS=500000
DINO_SEQ_PERF_ITERS=20
DINO_SEQ_PERF_READ_LEN=150
DINO_SEQ_PERF_OUT=target/perf
CARGO_TARGET_DIR=target
```

The wrapper uses the stable `throughput` bench target for wall-clock runs.
`gungraun-install` installs the matching `gungraun-runner` binary under
`tools/perf/gungraun/.bin/`. `gungraun` runs the standalone sidecar package in
`tools/perf/gungraun/`, writes `target/perf/gungraun-hotpaths.txt`, and reports
Callgrind-style instruction counts for small before/after checks. `perf`,
`hyperfine`, `cargo flamegraph`, and `gungraun` are optional host/tooling
surfaces, not normal crate dependencies. It defaults `CARGO_TARGET_DIR` to the
crate-local `target/` directory so the build cache is also removable and works
from restricted checkouts.
