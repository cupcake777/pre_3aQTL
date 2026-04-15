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
R_BIN="${R_BIN:-Rscript}"
VCF_FILE="${VCF_FILE:-/home/lyc/share_group_folder/raw_data/CHB_all_anno.vcf.gz}"

usage() {
  cat <<'EOF'
Usage: bash 02_downsample.sh [--stages "stage1 stage2"] [--threads 8] [--vcf /path/to.vcf.gz]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stages) STAGES="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --maf) MAF_THR="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --vcf) VCF_FILE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p downsample/pre_input downsample/nominal downsample/permutation downsample/susie downsample/logs

if command -v micromamba >/dev/null 2>&1; then
  set +u
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate "$CONDA_ENV"
  set -u
fi

"$R_BIN" downsample/01_prepare_downsample_inputs.R

for stage in $STAGES; do
  stage_label="${stage/stage/Stage}"
  sample_file="downsample/pre_input/sample.${stage_label}.downsample"
  genotype_prefix="downsample/pre_input/genotype.${stage}.downsample"
  phenotype_bed="downsample/pre_input/phenotype.${stage_label}.bed"
  phenotype_gz="${phenotype_bed}.gz"
  covariates="downsample/pre_input/covariates_for_qtl.${stage_label}.txt"
  nominal_prefix="downsample/nominal/${stage}"
  perm_prefix="downsample/permutation/${stage}"
  susie_prefix="downsample/susie/${stage}"

  [[ -f "$sample_file" ]] || { echo "[ERROR] Missing downsample sample file: $sample_file" >&2; exit 1; }

  if [[ ! -f "${genotype_prefix}.pgen" ]]; then
    plink2 --vcf "$VCF_FILE" \
      --keep "$sample_file" \
      --make-pgen \
      --out "$genotype_prefix" \
      --output-chr chrM \
      --threads "$THREADS" \
      >"downsample/logs/${stage}.plink2.log" 2>&1
  fi

  if [[ ! -f "$phenotype_gz" ]]; then
    bgzip -f "$phenotype_bed"
    tabix -p bed -f "$phenotype_gz"
  fi

  if [[ ! -f "${nominal_prefix}.cis_qtl_pairs.txt.gz" ]]; then
    PYTHONWARNINGS="ignore" "$PYTHON_BIN" -m tensorqtl \
      --mode cis_nominal \
      --seed "$SEED" \
      --maf_threshold "$MAF_THR" \
      --covariates "$covariates" \
      "$genotype_prefix" "$phenotype_gz" "$nominal_prefix" \
      >"downsample/logs/${stage}.tensorqtl.cis_nominal.log" 2>&1

    "$PYTHON_BIN" scripts/merge_tensorqtl_nominal.py \
      --prefix "$nominal_prefix" \
      --output "${nominal_prefix}.cis_qtl_pairs.txt.gz"
  else
    echo "[INFO] Reusing existing downsample nominal output for ${stage}"
  fi

  if [[ ! -f "${perm_prefix}.cis_qtl.txt.gz" ]]; then
    PYTHONWARNINGS="ignore" "$PYTHON_BIN" -m tensorqtl \
      --mode cis \
      --seed "$SEED" \
      --maf_threshold "$MAF_THR" \
      --covariates "$covariates" \
      "$genotype_prefix" "$phenotype_gz" "$perm_prefix" \
      >"downsample/logs/${stage}.tensorqtl.cis.log" 2>&1
  else
    echo "[INFO] Reusing existing downsample permutation output for ${stage}"
  fi

  if [[ ! -f "${perm_prefix}.sig_egenes.fdr05.txt" ]]; then
    "$PYTHON_BIN" scripts/summarize_permutation.py \
      --input "${perm_prefix}.cis_qtl.txt.gz" \
      --sig-output "${perm_prefix}.sig_egenes.fdr05.txt"
  fi

  if [[ ! -f "downsample/${stage}_sig_QTL.txt.gz" ]]; then
    "$PYTHON_BIN" scripts/filter_significant_pairs.py \
      --permutation "${perm_prefix}.cis_qtl.txt.gz" \
      --nominal "${nominal_prefix}.cis_qtl_pairs.txt.gz" \
      --output "downsample/${stage}_sig_QTL.txt.gz"
  fi

  if [[ ! -f "${susie_prefix}.SuSiE_summary.txt.gz" ]]; then
    PYTHONWARNINGS="ignore" "$PYTHON_BIN" -m tensorqtl \
      --mode cis_susie \
      --seed "$SEED" \
      --maf_threshold "$MAF_THR" \
      --covariates "$covariates" \
      --cis_output "${perm_prefix}.cis_qtl.txt.gz" \
      "$genotype_prefix" "$phenotype_gz" "$susie_prefix" \
      >"downsample/logs/${stage}.tensorqtl.cis_susie.log" 2>&1

    if [[ -f "${susie_prefix}.SuSiE.pickle" ]]; then
      "$PYTHON_BIN" scripts/export_susie_summary.py \
        --pickle "${susie_prefix}.SuSiE.pickle" \
        --output "${susie_prefix}.SuSiE_summary.txt.gz"
    fi
  else
    echo "[INFO] Reusing existing downsample SuSiE summary for ${stage}"
  fi
done
