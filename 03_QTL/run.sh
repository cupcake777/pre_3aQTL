#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

STATUS_DIR="${STATUS_DIR:-pipeline_status}"
LOG_DIR="${LOG_DIR:-logs}"
RUN_TS="${RUN_TS:-$(date +%Y%m%d_%H%M%S)}"
RUN_LOG="${LOG_DIR}/full_pipeline_${RUN_TS}.log"
STATUS_FILE="${STATUS_DIR}/full_pipeline_${RUN_TS}.status"
LATEST_LINK="${STATUS_DIR}/latest.status"

mkdir -p "$STATUS_DIR" "$LOG_DIR"

log() {
  local msg="$1"
  local ts
  ts="$(date '+%Y-%m-%d %H:%M:%S')"
  printf '[%s] %s\n' "$ts" "$msg" | tee -a "$RUN_LOG"
}

write_status() {
  local state="$1"
  local step="${2:-}"
  local extra="${3:-}"
  {
    printf 'run_ts=%s\n' "$RUN_TS"
    printf 'state=%s\n' "$state"
    printf 'step=%s\n' "$step"
    printf 'pid=%s\n' "$$"
    printf 'updated_at=%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    printf 'log=%s\n' "$RUN_LOG"
    if [[ -n "$extra" ]]; then
      printf 'detail=%s\n' "$extra"
    fi
  } >"$STATUS_FILE"
  ln -sfn "$(basename "$STATUS_FILE")" "$LATEST_LINK"
}

run_step() {
  local step_name="$1"
  local cmd="$2"
  write_status "running" "$step_name"
  log "START ${step_name}: ${cmd}"
  if bash -lc "$cmd" >>"$RUN_LOG" 2>&1; then
    log "DONE ${step_name}"
  else
    local exit_code=$?
    log "FAIL ${step_name} exit_code=${exit_code}"
    write_status "failed" "$step_name" "exit_code=${exit_code}"
    exit "$exit_code"
  fi
}

write_status "starting" "bootstrap"
log "Pipeline runner started"

run_step "01_stage_QTL" "bash 01_stage_QTL.sh"
run_step "02_downsample" "bash 02_downsample.sh"
run_step "03_compare" "bash 03_compare.sh"
run_step "04_interaction" "bash 04_interaction.sh"
run_step "05_integrate" "bash 05_integrate.sh"

write_status "completed" "done"
log "Pipeline runner completed successfully"
