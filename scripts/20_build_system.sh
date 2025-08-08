#!/usr/bin/env bash
set -euo pipefail
source ./config/local.env

SAMPLE_ID=${1:-}
if [[ -z "$SAMPLE_ID" ]]; then
  echo "Usage: $0 <sample_id>" >&2
  exit 1
fi

CG_DIR_SAMPLE="data/cg/${SAMPLE_ID}"
SYSTEM_DIR="data/systems/${SAMPLE_ID}"
mkdir -p "$SYSTEM_DIR"

INSANE=$(bash scripts/tools/ensure_insane.sh)

# Komposition lesen und Strings für -u und -l bauen (z.B. POPC:22,PSM:18,CHOL:60)
UPP=$(awk -F, 'NR>1 && $1=="upper" {printf (NR>2?",":""); printf "%s:%s", $2,$3}' "$LIPID_COMPOSITION_CSV")
LOW=$(awk -F, 'NR>1 && $1=="lower" {printf (NR>2?",":""); printf "%s:%s", $2,$3}' "$LIPID_COMPOSITION_CSV")

python3 "$INSANE" -u "$UPP" -l "$LOW" -x "$BOX_X_NM" -y "$BOX_Y_NM" -z "$BOX_Z_NM" \
  -o "$SYSTEM_DIR/system.gro" -p "$SYSTEM_DIR/system.top" -salt "$SALT_M" -sol W

echo "#include \"${CG_DIR_SAMPLE}/${SAMPLE_ID}_cg.itp\"" >> "$SYSTEM_DIR/system.top" || true

cp "$CG_DIR_SAMPLE/${SAMPLE_ID}_cg.pdb" "$SYSTEM_DIR/" || true

echo "[OK] B3 abgeschlossen: $SYSTEM_DIR"
