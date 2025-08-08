#!/usr/bin/env bash
set -euo pipefail
source ./config/local.env
SAMPLE_ID=${1:-}
REPL=${2:-1}
if [[ -z "$SAMPLE_ID" ]]; then echo "Usage: $0 <sample_id> [replicate]" >&2; exit 1; fi

SYSTEM_DIR="data/systems/${SAMPLE_ID}"
MD_DIR_SAMPLE="data/md/${SAMPLE_ID}/rep${REPL}"
mkdir -p "$MD_DIR_SAMPLE"

pushd "$SYSTEM_DIR" >/dev/null
  gmx grompp -f ../../mdp/prod_demo.mdp -c equil/npt.gro -p system.top -o "$MD_DIR_SAMPLE/prod.tpr" -maxwarn 2
popd >/dev/null

pushd "$MD_DIR_SAMPLE" >/dev/null
  gmx mdrun -deffnm prod
popd >/dev/null

echo "[OK] MD-Produktion abgeschlossen: $MD_DIR_SAMPLE"
