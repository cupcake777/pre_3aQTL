#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

CONDA_ENV="${CONDA_ENV:-qtl}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
R_BIN="${R_BIN:-Rscript}"
THREADS="${THREADS:-8}"
MAF_THR="${MAF_THR:-0.05}"
SEED="${SEED:-2026}"
# Interaction encoding mode:
#   linear   — LifeStage encoded as 0/1/2 (default). Assumes the genetic effect changes
#              *linearly* across stages. High power for monotonic (Up/Down) patterns.
#              Limited power for non-linear (Peak/Valley) patterns identified in 02_pattern.
#   contrast — Two dummy variables (Stage2_vs_Stage1, Stage3_vs_Stage1) capture
#              non-linear effects independently. Requires running TensorQTL twice
#              (once per contrast column) and combining results.
INTERACTION_MODE="${INTERACTION_MODE:-linear}"

usage() {
  cat <<'EOF'
Usage: bash 04_interaction.sh [--interaction-mode linear|contrast]

Options:
  --interaction-mode  Encoding for the interaction term (default: linear).
                      linear   = LifeStage 0/1/2; detects monotonic effects.
                      contrast = dummy variables Stage2_vs_S1 and Stage3_vs_S1;
                                 also detects Peak/Valley non-linear effects.

Expected input files:
  ../01_pre/after.combat.txt
  ../01_pre/covariates_for_qtl.Stage{1,2,3}.txt
  ../../raw_data/genotype.{pgen,psam,pvar}
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --interaction-mode) INTERACTION_MODE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ "$INTERACTION_MODE" != "linear" && "$INTERACTION_MODE" != "contrast" ]]; then
  echo "[ERROR] --interaction-mode must be 'linear' or 'contrast'" >&2; exit 1
fi

mkdir -p interaction/pre_input interaction/nominal interaction/logs

if command -v micromamba >/dev/null 2>&1; then
  set +u
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate "$CONDA_ENV"
  set -u
fi

genotype="interaction/pre_input/genotype.interaction"
phenotype_bed="interaction/pre_input/phenotype.interaction.bed"
phenotype="${phenotype_bed}.gz"
covariates="interaction/pre_input/covariates_for_qtl.interaction.txt"
if [[ "$INTERACTION_MODE" == "contrast" ]]; then
  interaction_terms="interaction/pre_input/interaction_term_contrast.txt"
  echo "[INFO] Interaction mode: contrast (Stage2_vs_S1, Stage3_vs_S1 dummy variables)"
else
  interaction_terms="interaction/pre_input/interaction_term.txt"
  echo "[INFO] Interaction mode: linear (LifeStage 0/1/2)"
fi
nominal_prefix="interaction/nominal/interaction"
top_assoc="${nominal_prefix}.cis_qtl_top_assoc.txt.gz"
chunk_dir="interaction/nominal/chunks"
chunk_bed_dir="interaction/pre_input/chrom_beds"

"$PYTHON_BIN" scripts/prepare_interaction_inputs.py \
  --expression "../01_pre/after.combat.txt" \
  --cov-dir "../01_pre" \
  --genotype-prefix "../../raw_data/genotype" \
  --out-dir "interaction/pre_input" \
  >"interaction/logs/interaction.prepare_inputs.log" 2>&1

for path in "${genotype}.pgen" "$phenotype_bed" "$covariates" "$interaction_terms"; do
  [[ -f "$path" ]] || { echo "[ERROR] Missing required interaction input: $path" >&2; exit 1; }
done

if [[ ! -f "$phenotype" ]]; then
  bgzip -f "$phenotype_bed"
  tabix -p bed -f "$phenotype"
fi

if [[ ! -f "$top_assoc" ]]; then
  mkdir -p "$chunk_dir" "$chunk_bed_dir"
  rm -f "${chunk_dir}/"*.cis_qtl_top_assoc.txt.gz
  mapfile -t chroms < <(tabix -l "$phenotype")
  [[ ${#chroms[@]} -gt 0 ]] || { echo "[ERROR] No contigs found in $phenotype" >&2; exit 1; }
  running_jobs=0
  for chrom in "${chroms[@]}"; do
    chunk_bed="${chunk_bed_dir}/phenotype.interaction.${chrom}.bed.gz"
    chunk_tag="${chrom//[^[:alnum:]_.-]/_}"
    chunk_prefix="${chunk_dir}/interaction.${chunk_tag}"
    chunk_top_assoc="${chunk_prefix}.cis_qtl_top_assoc.txt.gz"
    chunk_log="interaction/logs/interaction.tensorqtl.${chunk_tag}.log"

    if [[ ! -f "$chunk_bed" ]]; then
      tabix -h "$phenotype" "$chrom" | bgzip -c >"$chunk_bed"
      tabix -p bed -f "$chunk_bed"
    fi

    (
      PYTHONWARNINGS="ignore" "$PYTHON_BIN" -m tensorqtl \
        --mode cis_nominal \
        --best_only \
        --seed "$SEED" \
        --maf_threshold 0 \
        --maf_threshold_interaction "$MAF_THR" \
        --covariates "$covariates" \
        --interaction "$interaction_terms" \
        "$genotype" "$chunk_bed" "$chunk_prefix" \
        >"$chunk_log" 2>&1
    ) &

    running_jobs=$((running_jobs + 1))
    if (( running_jobs >= THREADS )); then
      wait -n
      running_jobs=$((running_jobs - 1))
    fi
  done
  wait || { echo "[ERROR] One or more background tensorqtl chunk jobs failed" >&2; exit 1; }

  "$PYTHON_BIN" scripts/merge_interaction_top_assoc.py \
    --inputs "${chunk_dir}"/interaction.*.cis_qtl_top_assoc.txt.gz \
    --output "$top_assoc"
fi

"$PYTHON_BIN" scripts/summarize_interaction_top_assoc.py \
  --input "$top_assoc" \
  --output "interaction/interaction_sig_QTL.txt.gz"

# ---------------------------------------------------------------------------
# Generate full cis-nominal pairs (no --best_only) for colocalization (05_coloc).
# coloc.abf requires complete summary statistics across all variants in a locus,
# not just the top association per gene.  We only re-run chromosomes that contain
# at least one significant interaction gene to limit computational cost.
# ---------------------------------------------------------------------------
all_pairs_out="interaction/interaction_all_pairs.txt.gz"
all_pairs_dir="interaction/nominal_all"

if [[ ! -f "$all_pairs_out" ]]; then
  mkdir -p "$all_pairs_dir"

  # Collect chromosomes that have significant interaction genes
  sig_chroms=$(zcat "interaction/interaction_sig_QTL.txt.gz" \
    | awk -F'\t' 'NR>1 {print $1}' \
    | sort -u)
  [[ -n "$sig_chroms" ]] || { echo "[WARN] No significant interaction chromosomes found; skipping all-pairs run." >&2; }

  running_jobs=0
  while IFS= read -r chrom; do
    chunk_bed="${chunk_bed_dir}/phenotype.interaction.${chrom}.bed.gz"
    chunk_tag="${chrom//[^[:alnum:]_.-]/_}"
    all_prefix="${all_pairs_dir}/interaction.${chunk_tag}"
    all_log="interaction/logs/interaction.tensorqtl.all.${chunk_tag}.log"

    if [[ ! -f "$chunk_bed" ]]; then
      tabix -h "$phenotype" "$chrom" | bgzip -c >"$chunk_bed"
      tabix -p bed -f "$chunk_bed"
    fi

    (
      PYTHONWARNINGS="ignore" "$PYTHON_BIN" -m tensorqtl \
        --mode cis_nominal \
        --seed "$SEED" \
        --maf_threshold 0 \
        --maf_threshold_interaction "$MAF_THR" \
        --covariates "$covariates" \
        --interaction "$interaction_terms" \
        "$genotype" "$chunk_bed" "$all_prefix" \
        >"$all_log" 2>&1
    ) &

    running_jobs=$((running_jobs + 1))
    if (( running_jobs >= THREADS )); then
      wait -n
      running_jobs=$((running_jobs - 1))
    fi
  done <<< "$sig_chroms"
  wait || { echo "[ERROR] One or more all-pairs tensorqtl jobs failed" >&2; exit 1; }

  "$PYTHON_BIN" scripts/merge_tensorqtl_nominal.py \
    --prefix "${all_pairs_dir}/interaction" \
    --output "$all_pairs_out"
else
  echo "[INFO] Reusing existing interaction all-pairs output: $all_pairs_out"
fi

"$R_BIN" interaction/01_interaction_visualize.R
