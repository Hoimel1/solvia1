#!/usr/bin/env bash
set -euo pipefail
source ./config/local.env
SAMPLE_ID=${1:-}
if [[ -z "$SAMPLE_ID" ]]; then echo "Usage: $0 <sample_id>" >&2; exit 1; fi

SYSTEM_DIR="data/systems/${SAMPLE_ID}"
MD_OUT_DIR="$SYSTEM_DIR/equil"
mkdir -p "$MD_OUT_DIR"

pushd "$SYSTEM_DIR" >/dev/null
  gmx grompp -f ../../mdp/minim.mdp -c system.gro -p system.top -o "$MD_OUT_DIR/minim.tpr" -maxwarn 2
  gmx mdrun -deffnm "$MD_OUT_DIR/minim"

  gmx grompp -f ../../mdp/nvt.mdp -c "$MD_OUT_DIR/minim.gro" -p system.top -o "$MD_OUT_DIR/nvt.tpr" -maxwarn 2
  gmx mdrun -deffnm "$MD_OUT_DIR/nvt"

  gmx grompp -f ../../mdp/npt.mdp -c "$MD_OUT_DIR/nvt.gro" -p system.top -o "$MD_OUT_DIR/npt.tpr" -maxwarn 2
  gmx mdrun -deffnm "$MD_OUT_DIR/npt"
popd >/dev/null

echo "[OK] Equilibration abgeschlossen: $MD_OUT_DIR"
