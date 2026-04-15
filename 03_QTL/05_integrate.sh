#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

R_BIN="${R_BIN:-Rscript}"

mkdir -p integrate
(cd integrate && "$R_BIN" 01_stage_vs_interaction.R)
