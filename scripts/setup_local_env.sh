#!/usr/bin/env bash
# ============================================================================
# setup_local_env.sh -- Lokale Installation aller Abhängigkeiten für Phase A
# Projekt: Solvia AMP-Tox Workflow (lokal-zuerst)
# ---------------------------------------------------------------------------
# Dieses Skript prüft Systemvoraussetzungen und installiert:
#   * Mambaforge (User-Space) + Conda-Umgebung "solvia_env"
#   * Python-Pakete (numpy, pandas, mdanalysis, xgboost, scikit-learn, shap, ...)
#   * Systempakete & Bibliotheken für Kompilierung
#   * GROMACS 2023.1 (GPU, CUDA 11.8, SM 89 für NVIDIA L4)
#   * DSSP (mkdssp)
#   * martinize2 + vermouth via pip (in der Conda-Umgebung)
#   * INSANE (Python2-Skript + Abhängigkeiten via "conda run -n py2_env")
# Weitere Tools (AlphaFold 3) werden später ggf. ergänzt.
# ----------------------------------------------------------------------------
# Aufruf:  bash scripts/setup_local_env.sh | tee reports/logs/setup_env_$(date +"%Y%m%d_%H%M%S").log
# ============================================================================
set -euo pipefail

## 0. Vorbereitungen ------------------------------------------------------------------
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR="$(dirname "${SCRIPT_DIR}")"
cd "${ROOT_DIR}"

printf "\n=== Solvia Lokale Umgebung – Setup startet ($(date)) ===\n"

## 1. System-Info ----------------------------------------------------------------------
uname -a
nproc --all || true
if command -v nvidia-smi &>/dev/null; then
  nvidia-smi -L || true
  nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv
else
  echo "[WARN] Keine NVIDIA-GPU gefunden – GROMACS wird ohne GPU gebaut."
fi

## 2. Abhängigkeiten: apt Packages -----------------------------------------------------
# sudo wird benötigt; Abbruch, falls nicht vorhanden
if ! command -v sudo &>/dev/null; then
  echo "[ERROR] sudo nicht verfügbar. Bitte mit Root-Rechten ausführen." >&2
  exit 1
fi

sudo apt-get update -y
sudo apt-get upgrade -y
sudo apt-get install -y build-essential cmake git wget curl libfftw3-dev libgsl-dev \
                        libhwloc-dev libxml2-dev libopenmpi-dev openmpi-bin ninja-build

## 3. Mambaforge installieren (falls nicht vorhanden) ---------------------------------
if ! command -v mamba &>/dev/null; then
  echo "[INFO] Installiere Mambaforge (user-space) …"
  MAMBA_INSTALLER="Mambaforge-Linux-x86_64.sh"
  wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/${MAMBA_INSTALLER}" -O /tmp/${MAMBA_INSTALLER}
  bash /tmp/${MAMBA_INSTALLER} -b -p "$HOME/mambaforge"
  eval "$($HOME/mambaforge/bin/conda shell.bash hook)"
  conda activate
fi

source "$HOME/mambaforge/etc/profile.d/conda.sh"

## 4. Conda-Umgebung erstellen --------------------------------------------------------
if ! conda info --envs | grep -q "solvia_env"; then
  echo "[INFO] Erstelle Conda-Umgebung solvia_env …"
  mamba create -y -n solvia_env -c conda-forge \
        python=3.11 numpy pandas mdanalysis xgboost scikit-learn shap \
        biopython seaborn jupyter matplotlib tqdm networkx dssp peptidebuilder
fi
conda activate solvia_env

## 5. martinize2, vermouth & weitere Python-Pakete -------------------------------------
pip install --upgrade pip
pip install martinize2==0.14.* vermouth==0.8.* PeptideBuilder

## 6. GROMACS 2023.1 kompilieren -------------------------------------------------------
GMX_VERSION="2023.1"
GMX_PREFIX="/usr/local/gromacs-${GMX_VERSION}"
if [ ! -d "${GMX_PREFIX}" ]; then
  echo "[INFO] Baue GROMACS ${GMX_VERSION} (GPU, CUDA 11.8) …"
  mkdir -p "$HOME/src" && cd "$HOME/src"
  git clone --branch v${GMX_VERSION} --depth 1 https://gitlab.com/gromacs/gromacs.git gromacs-${GMX_VERSION}
  mkdir -p gromacs-${GMX_VERSION}/build && cd gromacs-${GMX_VERSION}/build
  cmake .. \
    -DGMX_GPU=CUDA \
    -DGMX_CUDA_TARGET_SM=89 \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_DOUBLE=OFF \
    -DREGRESSIONTEST_DOWNLOAD=ON \
    -DCMAKE_INSTALL_PREFIX=${GMX_PREFIX}
  make -j$(nproc)
  sudo make install
fi
source ${GMX_PREFIX}/bin/GMXRC
cd "${ROOT_DIR}"

gmx --version
mkdssp -V

## 7. Abschluss & Verifikation ---------------------------------------------------------
echo "\n[INFO] Prüfe Binaries …"
which gmx && gmx --version
which mkdssp && mkdssp -V || echo "[WARN] mkdssp nicht im globalen PATH – in Conda-Umgebung verfügbar."

## 8. Skripte ausführbar machen --------------------------------------------------------
chmod +x "$ROOT_DIR"/scripts/*.sh || true
chmod +x "$ROOT_DIR"/scripts/tools/*.sh || true
chmod +x "$ROOT_DIR"/scripts/*.py || true

## 9. Fertig ---------------------------------------------------------------------------
echo "\n[OK] Installation abgeschlossen. Bitte neu einloggen oder 'source ${GMX_PREFIX}/bin/GMXRC' für gmx."

