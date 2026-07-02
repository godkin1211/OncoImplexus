#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/generate_patient_report.R"

if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "ERROR: Cannot find ${R_SCRIPT}" >&2
  exit 1
fi

exec Rscript "${R_SCRIPT}" "$@"
