#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

python3 05_coloc/01_prepare_inputs.py "$@"
bash 05_coloc/02_run_easycoloc.sh
