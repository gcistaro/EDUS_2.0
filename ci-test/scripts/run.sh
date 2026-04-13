#!/bin/bash
set -e

# ─────────────────────────────────────────────
# Locate project root from script location
# ─────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}/../.."

# ─────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────
BUILD_DIR="${PROJECT_ROOT}/build"
INPUT="${PROJECT_ROOT}/ci-test/inputs/hBN_IPA.json"
REF_DIR="${PROJECT_ROOT}/ci-test/outputs/hBN_IPA"
OUT_DIR="${PROJECT_ROOT}/CTEST/hBN_IPA"

# ─────────────────────────────────────────────
# Clean output directory
# ─────────────────────────────────────────────
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# ─────────────────────────────────────────────
# Run simulation
# ─────────────────────────────────────────────
"${BUILD_DIR}/EDUS" "${INPUT}"

# ─────────────────────────────────────────────
# Compare outputs (numerical regression)
# ─────────────────────────────────────────────
python3 "${PROJECT_ROOT}/ci-test/compare.py" \
    Output/Population.txt "${REF_DIR}/Population.txt"

python3 "${PROJECT_ROOT}/ci-test/compare.py" \
    Output/Velocity.txt "${REF_DIR}/Velocity.txt"


python3 "${PROJECT_ROOT}/ci-test/compare.py" \
    Output/Population_wannier.txt "${REF_DIR}/Population_wannier.txt"

ls "${PROJECT_ROOT}/PostProces/"

python3 "${PROJECT_ROOT}/PostProces/Absorbance.py" \
         "--smearing=0.6"

python3 "${PROJECT_ROOT}/ci-test/compare.py" \
    absorbance.txt "${REF_DIR}/absorbance.txt"
