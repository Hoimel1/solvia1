#!/usr/bin/env bash
set -euo pipefail

source ./config/local.env

SAMPLE_ID=${1:-}
if [[ -z "${SAMPLE_ID}" ]]; then
  echo "Usage: $0 <sample_id>" >&2
  exit 1
fi

# Sequenz aus CSV holen
SEQ=$(awk -F, -v id="$SAMPLE_ID" 'NR>1 && $1==id {print $2}' data/raw/peptides_demo.csv)
if [[ -z "${SEQ}" ]]; then
  echo "[ERROR] Kein Eintrag für ${SAMPLE_ID} in data/raw/peptides_demo.csv" >&2
  exit 2
fi

OUT_DIR="data/interim/${SAMPLE_ID}"
mkdir -p "$OUT_DIR"

# Helix-Fallback bauen
python3 scripts/build_helix.py --sequence "$SEQ" --out "$OUT_DIR/${SAMPLE_ID}_helix.pdb"

# Kurz-Minimierung in Vakuum (optional, als QC)
pushd "$OUT_DIR" >/dev/null
  echo "Protein" > index.ndx || true
  gmx pdb2gmx -f ${SAMPLE_ID}_helix.pdb -o ${SAMPLE_ID}_processed.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh <<EOF
1
1
EOF
  gmx grompp -f ../../mdp/minim.mdp -c ${SAMPLE_ID}_processed.gro -p topol.top -o em.tpr -maxwarn 2
  gmx mdrun -deffnm em -nt 1
popd >/dev/null

echo "[OK] B1 abgeschlossen: ${OUT_DIR}/${SAMPLE_ID}_helix.pdb"
