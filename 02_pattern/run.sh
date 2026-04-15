#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$script_dir"

echo "[1/3] trajectory clustering and ranked enrichment"
Rscript "$script_dir/01_trajectory_clustering.R" "$script_dir"

echo "[2/3] functional enrichment"
kegg_input="$(ls -1 "$script_dir"/05b_GSEA_KEGG_K*_all.csv 2>/dev/null | sort -V | tail -n 1 || true)"
if [[ -z "$kegg_input" ]]; then
  kegg_input="$script_dir/05b_GSEA_KEGG_K3_all.csv"
fi
Rscript "$script_dir/02_functional_enrichment.R" "$kegg_input" "$script_dir/GSEA"

echo "[3/3] figures"
Rscript "$script_dir/03_generate_figures.R"
