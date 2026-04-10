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
INPUT="${PROJECT_ROOT}/ci-test/inputs/hBN.json"
REF_DIR="${PROJECT_ROOT}/ci-test/outputs/hBN"
OUT_DIR="${BUILD_DIR}/CTEST/hBN"

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
    Output/DM0.txt "${REF_DIR}/DM0.txt"
