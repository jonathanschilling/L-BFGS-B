#!/usr/bin/env bash
# Regenerate the committed JSON reference outputs in tests/expected/.
# Run this after an intentional algorithm change, then review the diff
# before committing.
#
# Usage: bash tests/regenerate_expected.sh [build_dir]
#   build_dir defaults to <repo>/build
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${1:-${HERE}/../build}"

if [[ ! -d "${BUILD}" ]]; then
   echo "build directory not found: ${BUILD}" >&2
   echo "build the project first: cmake -S . -B build && cmake --build build -j" >&2
   exit 1
fi

mkdir -p "${HERE}/expected"

for drv in driver1_f77 driver1_f90 driver2_f77 driver2_f90 driver3_f77 driver3_f90; do
   if [[ ! -x "${BUILD}/${drv}" ]]; then
      echo "missing executable: ${BUILD}/${drv}" >&2
      exit 1
   fi
   wd="$(mktemp -d)"
   trap 'rm -rf "${wd}"' EXIT
   LBFGSB_JSON_OUTPUT="${wd}/${drv}.json" \
   LBFGSB_TLIMIT=86400 \
   LBFGSB_NFG_LIMIT=100000 \
   OPENBLAS_NUM_THREADS=1 \
   OMP_NUM_THREADS=1 \
      "${BUILD}/${drv}" >/dev/null
   cp "${wd}/${drv}.json" "${HERE}/expected/${drv}.json"
   rm -rf "${wd}"
   trap - EXIT
   echo "regenerated tests/expected/${drv}.json"
done
