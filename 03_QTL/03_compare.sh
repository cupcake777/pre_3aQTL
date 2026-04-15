#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

R_BIN="${R_BIN:-Rscript}"

mkdir -p compare
"$R_BIN" compare/01_combined_analysis.R
