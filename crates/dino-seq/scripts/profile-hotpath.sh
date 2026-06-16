#!/usr/bin/env bash
set -euo pipefail

records="${DINO_SEQ_PROFILE_RECORDS:-500000}"
read_len="${DINO_SEQ_PROFILE_READ_LEN:-150}"
iters="${DINO_SEQ_PROFILE_ITERS:-3}"
workers="${DINO_SEQ_WORKERS:-$(nproc)}"
input_dir="${DINO_SEQ_PROFILE_INPUT_DIR:-target/profile-inputs}"
profile_dir="${DINO_SEQ_PROFILE_DIR:-target/profiles}"
pattern="${DINO_SEQ_PROFILE_PATTERN:-cyclic}"
profile_bgzf_parallel="${DINO_SEQ_PROFILE_BGZF_PARALLEL:-0}"

mkdir -p "${input_dir}" "${profile_dir}"

cargo +nightly build --release --all-features --bin dino-seq-bench --bin dino-seq-fixture
target/release/dino-seq-fixture \
  --out-dir "${input_dir}" \
  --records "${records}" \
  --read-len "${read_len}" \
  --pattern "${pattern}"

input="${DINO_SEQ_PROFILE_INPUT:-${input_dir}/single.fastq}"
jsonl="${profile_dir}/dino_seq-hotpath.jsonl"
: > "${jsonl}"

for mode in parse pack; do
  profile_args=()
  if [[ "${mode}" == "pack" && "${profile_bgzf_parallel}" != "0" ]]; then
    profile_args+=(--profile-bgzf-parallel)
  fi

  printf 'benchmarking mode=%s input=%s\n' "${mode}" "${input}"
  target/release/dino-seq-bench \
    --input "${input}" \
    --mode "${mode}" \
    --iters "${iters}" \
    --workers "${workers}" \
    "${profile_args[@]}" \
    --json >> "${jsonl}"

  if command -v perf >/dev/null 2>&1; then
    perf stat -d \
      -o "${profile_dir}/dino_seq-${mode}.perf-stat.txt" \
      target/release/dino-seq-bench \
        --input "${input}" \
        --mode "${mode}" \
        --iters "${iters}" \
        --workers "${workers}" \
        "${profile_args[@]}" >/dev/null || true

    perf record -F 999 -g \
      -o "${profile_dir}/dino_seq-${mode}.perf.data" -- \
      target/release/dino-seq-bench \
        --input "${input}" \
        --mode "${mode}" \
        --iters 1 \
        --workers "${workers}" \
        "${profile_args[@]}" >/dev/null 2>&1 || true

    if [[ -f "${profile_dir}/dino_seq-${mode}.perf.data" ]]; then
      perf report \
        -i "${profile_dir}/dino_seq-${mode}.perf.data" \
        --stdio > "${profile_dir}/dino_seq-${mode}.perf-report.txt" || true
    fi
  fi
done

printf 'wrote %s\n' "${jsonl}"
printf 'wrote %s/dino_seq-parse.perf-stat.txt if perf was permitted\n' "${profile_dir}"
printf 'wrote %s/dino_seq-pack.perf-stat.txt if perf was permitted\n' "${profile_dir}"
