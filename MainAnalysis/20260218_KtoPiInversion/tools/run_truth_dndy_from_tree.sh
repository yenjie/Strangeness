#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <input_truth.root> <output_dndy.root> [template_dndy.root]" >&2
  exit 1
fi

INPUT_ROOT="$1"
OUTPUT_ROOT="$2"
TEMPLATE_ROOT="${3:-output/systematics_20260314_dndy/nominal_unfold_dndy.root}"

source /raid5/root/root-v6.34.04/root/bin/thisroot.sh

PYTHIA8_PREFIX="${PWD}/external/pythia8317-src/install"
PYTHIA8_CONFIG="${PYTHIA8_PREFIX}/bin/pythia8-config"
BUILD_DIR="${PWD}/tools/.build"
mkdir -p "${BUILD_DIR}"

EXE="${BUILD_DIR}/build_truth_kpi_vs_dndy_from_tree"
SRC="${PWD}/tools/build_truth_kpi_vs_dndy_from_tree.cc"

if [ ! -x "${EXE}" ] || [ "${SRC}" -nt "${EXE}" ]; then
  g++ -O2 -std=c++17 "${SRC}" -o "${EXE}" \
    $(root-config --cflags --libs) \
    $( ${PYTHIA8_CONFIG} --cxxflags --libs | sed 's/-std=c++11//g' ) \
    -std=c++17
fi

"${EXE}" "${INPUT_ROOT}" "${TEMPLATE_ROOT}" "${OUTPUT_ROOT}"
