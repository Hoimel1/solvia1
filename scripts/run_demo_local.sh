#!/usr/bin/env bash
set -euo pipefail

SAMPLE_ID=${1:-}
if [[ -z "$SAMPLE_ID" ]]; then echo "Usage: $0 <sample_id>" >&2; exit 1; fi

./scripts/00_predict_structure.sh "$SAMPLE_ID"
./scripts/10_cg.sh "$SAMPLE_ID"
./scripts/20_build_system.sh "$SAMPLE_ID"
./scripts/30_equilibrate.sh "$SAMPLE_ID"
./scripts/40_md.sh "$SAMPLE_ID" 1
./scripts/50_features.sh "$SAMPLE_ID" 1
./scripts/60_ml.sh
