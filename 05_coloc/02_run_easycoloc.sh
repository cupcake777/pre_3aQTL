#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

resolve_easycoloc_root() {
  if [[ -n "${EASYCOLOC_ROOT:-}" ]]; then
    if [[ -d "${EASYCOLOC_ROOT}" ]]; then
      printf '%s\n' "${EASYCOLOC_ROOT}"
      return 0
    fi
    echo "[ERROR] EASYCOLOC_ROOT does not exist: ${EASYCOLOC_ROOT}" >&2
    return 1
  fi

  for candidate in \
    "${ROOT_DIR}/../EasyColoc" \
    "/mnt/share_group_folder/work/EasyColoc" \
    "/home/lyc/share_group_folder/work/EasyColoc"; do
    if [[ -d "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done

  echo "[ERROR] Could not locate EasyColoc. Clone https://github.com/cupcake777/EasyColoc to ../EasyColoc or set EASYCOLOC_ROOT." >&2
  return 1
}

repo_root="$(resolve_easycoloc_root)"
results_dir="${ROOT_DIR}/05_coloc/results"
run_log="${results_dir}/run_easycoloc.log"
mkdir -p "${results_dir}"
cd "${repo_root}"
export EASYCOLOC_GLOBAL_CONFIG="${ROOT_DIR}/05_coloc/config/global.yaml"
export EASYCOLOC_GWAS_CONFIG="${ROOT_DIR}/05_coloc/config/gwas.yaml"
export EASYCOLOC_QTL_CONFIG="${ROOT_DIR}/05_coloc/config/qtl.yaml"
export EASYCOLOC_OUTPUT_DIR="${results_dir}"
export EASYCOLOC_LOG_FILE="${run_log}"
export EASYCOLOC_RUN_LABEL="pre3aqtl_05_coloc"
exec Rscript run_coloc.r > "${run_log}" 2>&1
