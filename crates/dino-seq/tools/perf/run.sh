#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
usage: tools/perf/run.sh <list|bench|perf|hyperfine|flamegraph|gungraun-install|gungraun>

env:
  DINO_SEQ_PERF_CASE       benchmark case, default: pack
  DINO_SEQ_PERF_RECORDS    synthetic record count, default: 500000
  DINO_SEQ_PERF_ITERS      iterations per run, default: 20
  DINO_SEQ_PERF_READ_LEN   read length, default: 150
  DINO_SEQ_PERF_OUT        output directory, default: target/perf
  CARGO_TARGET_DIR         cargo output directory, default: target
EOF
}

command_exists() {
  command -v "$1" >/dev/null 2>&1
}

cargo_target_dir() {
  cargo metadata --no-deps --format-version 1 \
    | sed -n 's/.*"target_directory":"\([^"]*\)".*/\1/p'
}

bench_bin() {
  local target_dir
  target_dir="$(cargo_target_dir)"
  find "$target_dir/release/deps" -maxdepth 1 -type f -executable -name 'throughput-*' \
    | sort \
    | tail -n 1
}

case_name="${DINO_SEQ_PERF_CASE:-pack}"
records="${DINO_SEQ_PERF_RECORDS:-500000}"
iters="${DINO_SEQ_PERF_ITERS:-20}"
read_len="${DINO_SEQ_PERF_READ_LEN:-150}"
out_dir="${DINO_SEQ_PERF_OUT:-target/perf}"
mode="${1:-}"

export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-target}"

if [ -z "$mode" ]; then
  usage
  exit 2
fi

mkdir -p "$out_dir"

case "$mode" in
  list)
    cargo bench --bench throughput --all-features -- --list
    ;;
  bench)
    cargo bench --bench throughput --all-features -- \
      --case "$case_name" \
      --records "$records" \
      --read-len "$read_len" \
      --iters "$iters" \
      | tee "$out_dir/bench-$case_name.txt"
    ;;
  perf)
    command_exists perf || {
      echo "perf is not installed" >&2
      exit 127
    }
    cargo bench --bench throughput --all-features --no-run >/dev/null
    bin="$(bench_bin)"
    [ -n "$bin" ] || {
      echo "throughput bench binary not found" >&2
      exit 1
    }
    perf stat -d -o "$out_dir/perf-stat-$case_name.txt" "$bin" \
      --case "$case_name" \
      --records "$records" \
      --read-len "$read_len" \
      --iters "$iters"
    cat "$out_dir/perf-stat-$case_name.txt"
    ;;
  hyperfine)
    command_exists hyperfine || {
      echo "hyperfine is not installed" >&2
      exit 127
    }
    cargo bench --bench throughput --all-features --no-run >/dev/null
    bin="$(bench_bin)"
    [ -n "$bin" ] || {
      echo "throughput bench binary not found" >&2
      exit 1
    }
    hyperfine --warmup 1 --export-json "$out_dir/hyperfine-$case_name.json" \
      "$bin --case $case_name --records $records --read-len $read_len --iters $iters"
    ;;
  flamegraph)
    command_exists cargo-flamegraph || command_exists flamegraph || {
      echo "cargo-flamegraph is not installed" >&2
      exit 127
    }
    cargo flamegraph --bench throughput --features bgzf,gzip \
      -o "$out_dir/flamegraph-$case_name.svg" -- \
      --case "$case_name" \
      --records "$records" \
      --read-len "$read_len" \
      --iters "$iters"
    ;;
  gungraun-install)
    cargo install --version 0.19.2 --root tools/perf/gungraun/.bin gungraun-runner
    ;;
  gungraun)
    command_exists valgrind || {
      echo "valgrind is not installed" >&2
      exit 127
    }
    runner="$PWD/tools/perf/gungraun/.bin/bin/gungraun-runner"
    if [ -x "$runner" ] && [ -z "${GUNGRAUN_RUNNER:-}" ]; then
      export GUNGRAUN_RUNNER="$runner"
    fi
    cargo bench --manifest-path tools/perf/gungraun/Cargo.toml --bench hotpaths 2>&1 \
      | tee "$out_dir/gungraun-hotpaths.txt"
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac
