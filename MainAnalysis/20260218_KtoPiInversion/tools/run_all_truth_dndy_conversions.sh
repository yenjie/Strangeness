#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
RUNNER="${ROOT_DIR}/tools/run_truth_dndy_from_tree.sh"
MAX_PARALLEL="${MAX_PARALLEL:-6}"

declare -a JOBS=(
  "result/20260314/pythia8_truth/pythia8_zpole_truth_400k.root result/20260314/pythia8_truth/pythia8_zpole_truth_400k_dndy.root"
  "result/20260314/pythia8_rope_truth/pythia8_zpole_truth_rope_400k.root result/20260314/pythia8_rope_truth/pythia8_zpole_truth_rope_400k_dndy.root"
  "result/20260314/pythia8_dire_truth/pythia8_zpole_truth_dire_400k.root result/20260314/pythia8_dire_truth/pythia8_zpole_truth_dire_400k_dndy.root"
  "result/20260314/herwig_truth/herwig_ee_zpole_truth_400k.root result/20260314/herwig_truth/herwig_ee_zpole_truth_400k_dndy.root"
  "result/20260314/sherpa_truth/sherpa_ee_zpole_truth_400k.root result/20260314/sherpa_truth/sherpa_ee_zpole_truth_400k_dndy.root"
  "result/20260314/xscape_colorless/xscape_epem_colorless_zpole_truth_20k.root result/20260314/xscape_colorless/xscape_epem_colorless_zpole_truth_20k_dndy.root"
  "result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth_clean20k.root result/20260314/xscape_hybrid/xscape_epem_hybrid_zpole_truth_clean20k_dndy.root"
)

mkdir -p "${ROOT_DIR}/logs/dndy_generators"

running=0
for job in "${JOBS[@]}"; do
  infile="${job%% *}"
  outfile="${job##* }"
  name="$(basename "${outfile}" .root)"
  (
    cd "${ROOT_DIR}"
    bash "${RUNNER}" "${infile}" "${outfile}" > "logs/dndy_generators/${name}.log" 2>&1
  ) &
  running=$((running + 1))
  if [ "${running}" -ge "${MAX_PARALLEL}" ]; then
    wait -n
    running=$((running - 1))
  fi
done
wait
