#!/usr/bin/env bash
set -euo pipefail

mkdir -p target/profiles

args=()
if [[ -n "${DINO_SEQ_INPUT:-}" ]]; then
  args+=(--input "${DINO_SEQ_INPUT}")
fi
if [[ -n "${DINO_SEQ_MODE:-}" ]]; then
  args+=(--mode "${DINO_SEQ_MODE}")
fi

cargo build --release --bin dino-seq-bench
perf stat -d -r "${DINO_SEQ_PERF_REPEATS:-3}" \
  target/release/dino-seq-bench \
  "${args[@]}" \
  --records "${DINO_SEQ_RECORDS:-500000}" \
  --read-len "${DINO_SEQ_READ_LEN:-150}" \
  --iters "${DINO_SEQ_ITERS:-3}" \
  --workers "${DINO_SEQ_WORKERS:-$(nproc)}"

perf record -F 999 -g -o target/profiles/dino-seq-bench.perf.data -- \
  target/release/dino-seq-bench \
  "${args[@]}" \
  --records "${DINO_SEQ_RECORDS:-500000}" \
  --read-len "${DINO_SEQ_READ_LEN:-150}" \
  --iters 1 \
  --workers "${DINO_SEQ_WORKERS:-$(nproc)}"

perf report -i target/profiles/dino-seq-bench.perf.data --stdio \
  > target/profiles/dino-seq-bench.perf.txt

echo "wrote target/profiles/dino-seq-bench.perf.data"
echo "wrote target/profiles/dino-seq-bench.perf.txt"
