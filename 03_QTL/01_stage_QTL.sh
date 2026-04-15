#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

CONDA_ENV="${CONDA_ENV:-qtl}"
THREADS="${THREADS:-8}"
MAF_THR="${MAF_THR:-0.05}"
SEED="${SEED:-2026}"
STAGES="${STAGES:-stage1 stage2 stage3}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

usage() {
  cat <<'EOF'
Usage: bash 01_stage_QTL.sh [--stages "stage1 stage2"] [--threads 8] [--maf 0.05] [--seed 2026]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stages) STAGES="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --maf) MAF_THR="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p nominal permutation susie logs

if command -v micromamba >/dev/null 2>&1; then
  set +u
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate "$CONDA_ENV"
  set -u
fi

run_stage() {
  local stage="$1"
  local stage_label="${stage/stage/Stage}"
  local genotype="../01_pre/genotype.${stage}"
  local phenotype="../01_pre/phenotype.${stage_label}.bed.gz"
  local covariates="../01_pre/covariates_for_qtl.${stage_label}.txt"
  local nominal_prefix="nominal/${stage}"
  local perm_prefix="permutation/${stage}"
  local susie_prefix="susie/${stage}"

  for path in "${genotype}.pgen" "$phenotype" "$covariates"; do
    [[ -f "$path" ]] || { echo "[ERROR] Missing required input: $path" >&2; exit 1; }
  done

  echo "========================================================"
  echo "Running stage-specific TensorQTL workflow for ${stage}"
  echo "========================================================"

  if [[ ! -f "${nominal_prefix}.cis_qtl_pairs.txt.gz" ]]; then
    PYTHONWARNINGS="ignore" "$PYTHON_BIN" -m tensorqtl \
      --mode cis_nominal \
      --seed "$SEED" \
      --maf_threshold "$MAF_THR" \
      --covariates "$covariates" \
      "$genotype" "$phenotype" "$nominal_prefix" \
      >"logs/${stage}.tensorqtl.cis_nominal.log" 2>&1

    "$PYTHON_BIN" scripts/merge_tensorqtl_nominal.py \
      --prefix "$nominal_prefix" \
      --output "${nominal_prefix}.cis_qtl_pairs.txt.gz"
  else
    echo "[INFO] Reusing existing nominal output for ${stage}"
  fi

  if [[ ! -f "${perm_prefix}.cis_qtl.txt.gz" ]]; then
    PYTHONWARNINGS="ignore" "$PYTHON_BIN" -m tensorqtl \
      --mode cis \
      --seed "$SEED" \
      --maf_threshold "$MAF_THR" \
      --covariates "$covariates" \
      "$genotype" "$phenotype" "$perm_prefix" \
      >"logs/${stage}.tensorqtl.cis.log" 2>&1
  else
    echo "[INFO] Reusing existing permutation output for ${stage}"
  fi

  if [[ ! -f "${perm_prefix}.sig_egenes.fdr05.txt" ]]; then
    "$PYTHON_BIN" scripts/summarize_permutation.py \
      --input "${perm_prefix}.cis_qtl.txt.gz" \
      --sig-output "${perm_prefix}.sig_egenes.fdr05.txt"
  fi

  if [[ ! -f "${stage}_sig_QTL.txt.gz" ]]; then
    "$PYTHON_BIN" scripts/filter_significant_pairs.py \
      --permutation "${perm_prefix}.cis_qtl.txt.gz" \
      --nominal "${nominal_prefix}.cis_qtl_pairs.txt.gz" \
      --output "${stage}_sig_QTL.txt.gz"
  fi

  if [[ ! -f "${susie_prefix}.SuSiE_summary.txt.gz" ]]; then
    PYTHONWARNINGS="ignore" "$PYTHON_BIN" -m tensorqtl \
      --mode cis_susie \
      --seed "$SEED" \
      --maf_threshold "$MAF_THR" \
      --covariates "$covariates" \
      --cis_output "${perm_prefix}.cis_qtl.txt.gz" \
      "$genotype" "$phenotype" "$susie_prefix" \
      >"logs/${stage}.tensorqtl.cis_susie.log" 2>&1

    if [[ -f "${susie_prefix}.SuSiE.pickle" ]]; then
      "$PYTHON_BIN" scripts/export_susie_summary.py \
        --pickle "${susie_prefix}.SuSiE.pickle" \
        --output "${susie_prefix}.SuSiE_summary.txt.gz"
    fi
  else
    echo "[INFO] Reusing existing SuSiE summary for ${stage}"
  fi
}

for stage in $STAGES; do
  run_stage "$stage"
done
