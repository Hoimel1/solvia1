#!/usr/bin/env bash
set -euo pipefail
TOOLS_DIR="scripts/tools"
INSANE_DIR="$TOOLS_DIR/insane"
if [[ ! -d "$INSANE_DIR" ]]; then
  mkdir -p "$INSANE_DIR"
  git clone --depth 1 https://github.com/marrink-lab/insane "$INSANE_DIR"
fi
if [[ ! -f "$INSANE_DIR/insane.py" ]]; then
  echo "[ERROR] insane.py nicht gefunden im geklonten Repo." >&2
  exit 1
fi
chmod +x "$INSANE_DIR/insane.py"
echo "$INSANE_DIR/insane.py"
