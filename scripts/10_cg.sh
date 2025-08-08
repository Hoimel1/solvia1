#!/usr/bin/env bash
set -euo pipefail
source ./config/local.env

SAMPLE_ID=${1:-}
if [[ -z "$SAMPLE_ID" ]]; then
  echo "Usage: $0 <sample_id>" >&2
  exit 1
fi

IN_PDB="data/interim/${SAMPLE_ID}/${SAMPLE_ID}_helix.pdb"
OUT_DIR="data/cg/${SAMPLE_ID}"
mkdir -p "$OUT_DIR"

martinize2 -f "$IN_PDB" -x "$OUT_DIR/${SAMPLE_ID}_cg.pdb" -o "$OUT_DIR/topol.top" \
  -ff "$MARTINI_FF" -water "$MARTINI_WATER" -dssp mkdssp --seed "$SEED"

# Optional: Position restraints und includes in top/
mkdir -p top

echo "[OK] B2 abgeschlossen: $OUT_DIR"
