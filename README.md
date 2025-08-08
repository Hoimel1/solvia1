# Solvia – Antimikrobielle Peptid-Toxizität (AMP-Tox) Workflow

Dieses Repository implementiert einen _lokal-zuerst_ Ansatz, um einen vollständigen End-to-End-Workflow zur Vorhersage der Membran-Toxizität antimikrobieller Peptide aufzubauen. Die Schritte orientieren sich an **docs/todos.md**.

## Schnelleinstieg

```bash
# 1. Repository klonen (falls nicht in diesem Verzeichnis)
# 2. Abhängigkeiten installieren
bash scripts/setup_local_env.sh

# 3. Demo-Lauf (nach erfolgreicher Installation)
./scripts/run_demo_local.sh buforin_II   # (wird in Phase C erzeugt)
```

## Verzeichnisstruktur (lokal)

```
data/
  raw/            # Peptidliste, Metadaten
  interim/        # Zwischenstufen (AF3-Modelle, RMSD etc.)
  cg/             # Coarse-Grain-Strukturen & Topologien
  systems/        # Membransysteme nach INSANE
  md/             # Trajektorien & Energies
  features/       # Extracted Features (CSV/Parquet)
  models/         # ML-Modelle
  reports/        # Visualisierungen, Logs
scripts/          # Step-wise Shell- und Python-Skripte
mdp/              # GROMACS-Parameterdateien
config/           # env- oder YAML-Konfigs
```

## Wissenschaftlicher Hintergrund

1. **Strukturvorhersage:** AlphaFold 3 oder Helix-Heuristik liefert Atommodelle kleiner Peptide.
2. **Coarse-Graining:** Mit `martinize2` werden Martini 3‐Topologien erzeugt, um MD-Zeiträume >10 µs abzudecken.
3. **Membranbau:** `INSANE` erstellt cholesterinreiche, asymmetrische Bilayer gemäß RBC-Nachbildung.
4. **Molekulardynamik:** GROMACS 2023/2024 (GPU) simuliert Einbettung & Dynamik (20–50 ns für Demo).
5. **Feature-Extraktion:** Geometrische u. energetische Deskriptoren (Insertionstiefe, H-Bonding etc.) via MDAnalysis.
6. **ML-Modellierung:** Gradient-Boosted Trees (XGBoost) zur Klassifikation/Regression von HC50.

Alle Parameter (Seeds, Force-Field, Cut-Offs) sind reproduzierbar dokumentiert. Siehe **docs/todos.md** für QC-Kriterien.

## Lizenz & Zitation

Dieses Projekt steht unter einer CC-BY-SA 4.0 Lizenz. Bitte zitiere entsprechende Primärliteratur für GROMACS, Martini 3, AlphaFold 3 und die verwendeten Lipidparametrisierungen.
