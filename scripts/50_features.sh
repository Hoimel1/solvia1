#!/usr/bin/env bash
set -euo pipefail
source ./config/local.env
SAMPLE_ID=${1:-}
REPL=${2:-1}
if [[ -z "$SAMPLE_ID" ]]; then echo "Usage: $0 <sample_id> [replicate]" >&2; exit 1; fi

MD_DIR_SAMPLE="data/md/${SAMPLE_ID}/rep${REPL}"
OUT_DIR="data/features/replicates/${SAMPLE_ID}/rep${REPL}"
mkdir -p "$OUT_DIR"

# Trajektorie zentrieren und PBC behandeln
pushd "$MD_DIR_SAMPLE" >/dev/null
  echo 0 0 | gmx trjconv -s prod.tpr -f prod.xtc -o prod_center.xtc -pbc mol -center -ur compact
popd >/dev/null

python3 scripts/extract_features.py \
  --tpr "$MD_DIR_SAMPLE/prod.tpr" \
  --xtc "$MD_DIR_SAMPLE/prod_center.xtc" \
  --lipid_csv "$LIPID_COMPOSITION_CSV" \
  --out_csv "$OUT_DIR/features.csv" \
  --last_ns "$FEATURE_WINDOW_NS"

echo "[OK] Feature-Extraktion abgeschlossen: $OUT_DIR/features.csv"
