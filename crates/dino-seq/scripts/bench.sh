#!/usr/bin/env bash
set -euo pipefail

args=()
if [[ -n "${DINO_SEQ_INPUT:-}" ]]; then
  args+=(--input "${DINO_SEQ_INPUT}")
fi
if [[ -n "${DINO_SEQ_MODE:-}" ]]; then
  args+=(--mode "${DINO_SEQ_MODE}")
fi

cargo run --release --bin dino-seq-bench -- \
  "${args[@]}" \
  --records "${DINO_SEQ_RECORDS:-200000}" \
  --read-len "${DINO_SEQ_READ_LEN:-150}" \
  --iters "${DINO_SEQ_ITERS:-7}" \
  --workers "${DINO_SEQ_WORKERS:-$(nproc)}"
