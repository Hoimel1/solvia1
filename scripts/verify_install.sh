#!/usr/bin/env bash
set -euo pipefail
LOG="reports/logs/proof_of_install.txt"
mkdir -p reports/logs
{
  echo "=== Verify Install $(date) ==="
  which gmx && gmx --version
  which mkdssp && mkdssp -V || true
  python3 -c "import numpy, pandas, MDAnalysis, xgboost, sklearn, shap; print('Python-Pakete OK')"
  echo "INSANE Pfad:"; bash scripts/tools/ensure_insane.sh
} | tee "$LOG"
