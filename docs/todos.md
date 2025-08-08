## To‑Do‑Liste (Lokal zuerst): End‑to‑End auf einer Maschine lauffähig, danach Docker & Automatisierung

Leitprinzip: Zuerst beweisen, dass der komplette Workflow lokal auf dieser Maschine (ohne Docker, ohne Snakemake) läuft. Danach schrittweise Dockerisierung und Automatisierung. Erst wenn der lokale E2E‑Lauf stabil ist, werden Container, Snakemake, CI/CD und Cloud angegangen.

- **Zielbild lokal (Definition):** Für mindestens 1–2 Beispielpeptide werden alle Schritte (Struktur → CG → Membranbau/Equilibration → MD → Feature‑Extraktion → Aggregation/ML) manuell über lokale Tools ausgeführt, Artefakte abgelegt und QC‑Kriterien erfüllt.

⸻

## Phase A: Lokale Voraussetzungen (ohne Container)

- **Systemabhängigkeiten**
  - [ ] Python ≥ 3.10 inkl. `numpy`, `pandas`, `mdanalysis`, `xgboost`, `scikit-learn`, `shap`
  - [ ] GROMACS 2023/2024 lokal installiert (GPU falls vorhanden); `gmx` im PATH
  - [ ] DSSP (`mkdssp`) im PATH für martinize2
  - [ ] martinize2 (v0.14.x) lokal installierbar (z.B. via pip/conda); Vermouth passend
  - [ ] INSANE (`insane.py`, Python2‑Runtime falls benötigt) lokal verfügbar
  - [ ] Optional: AlphaFold 3 lokal (GPU) ODER Fallback (helical builder/RoseTTAFold) für Demo

- **Projektstruktur (lokal minimal)**
  - [ ] Ordner: `data/{raw,interim,cg,systems,md,features,models,reports}`, `scripts/`, `mdp/`, `top/`, `lipids/`
  - [ ] Beispiel‑Eingabedateien: `data/raw/peptides_demo.csv` (2–3 Peptide mit HC50)
  - [ ] MDP‑Vorlagen anlegen: `mdp/minim.mdp`, `mdp/nvt.mdp`, `mdp/npt.mdp`, `mdp/prod_demo.mdp` (20–50 ns)

Akzeptanzkriterium: Alle Binaries sind aufrufbar (`gmx -version`, `mkdssp -V`), Python‑Pakete importierbar, Ordner existieren.

⸻

## Phase B: Manueller End‑to‑End‑Lauf (ohne Docker/Snakemake)

- **B1 Strukturvorhersage**
  - [ ] Für jedes Demo‑Peptid Struktur generieren (AF3 falls verfügbar; sonst einfacher Helix‑Builder)
  - [ ] QC: pLDDT/Heuristik protokollieren; Minimierung in Vakuum (kurz) durchführen
  - [ ] Artefakte: `data/interim/<id>/ranked_1.pdb` (+ Metadaten)

- **B2 Coarse‑Graining (martinize2, Martini 3)**
  - [ ] martinize2 Aufruf lokal mit fixen Seeds (42), `-ff martini3001`, `-water martini3`, `-dssp mkdssp`
  - [ ] Optional: adaptive ElNeDyn, IDP‑Modus je nach Peptid
  - [ ] QC: kurze Minimierung (Fmax < 1000 kJ/mol/nm), Backmapping‑RMSD ≤ 0.45 nm (wenn Tool vorhanden)
  - [ ] Artefakte: `data/cg/<id>/*_cg.pdb`, `*.top`, `*.itp`

- **B3 Membranbau (INSANE) + Equilibration (GROMACS)**
  - [ ] RBC‑ähnliche, cholesterinreiche, asymmetrische Bilayer gemäß `docs/concept.md` lokal bauen
  - [ ] Minimierung → 100 ps NVT (310 K, posres Peptid) → 5–10 ns NPT (semi‑isotrop)
  - [ ] QC: Dicke ~4.5–5.0 nm; APL plausibel; PS innen; keine Overlaps
  - [ ] Artefakte: `data/systems/<id>/system.gro`, `system.top`, Equilibration‑Outputs

- **B4 Produktions‑MD (Demo)**
  - [ ] 1–2 Replikate je Demo‑Peptid, je 20–50 ns (lokal realistisch), dt=0.02 ps, Cutoff 1.1 nm, RF εr=15
  - [ ] Energy‑Groups `Protein`/`Membrane` aktivieren; Checkpoints optional
  - [ ] Artefakte: `data/md/<id>/prod.xtc`, `.edr`, `.tpr`

- **B5 Trajektorien‑Vorbereitung & Feature‑Extraktion**
  - [ ] `gmx trjconv`: `-pbc mol -center -ur compact`; Membran auf z ausrichten falls nötig
  - [ ] `scripts/extract_features.py` lokal ausführen (letzte 5–10 ns für Demo) → CSV/Parquet
  - [ ] Features prüfen (Bereich/Einheiten); Artefakte nach `data/features/replicates/`

- **B6 Aggregation & Minimal‑ML**
  - [ ] Replikate mitteln; z‑Score Normalisierung; kleiner XGBoost‑Classifier/Regressor testweise
  - [ ] Ziel: End‑to‑End Skriptausführung ohne Fehler; Beispielmodelle unter `models/`

Akzeptanzkriterium: Für ≥1 Peptid liegen durchgängige Artefakte aller Stufen vor und Feature‑CSV wird erzeugt. Mindestens ein kurzer MD‑Lauf completes lokal ohne Instabilität.

⸻

## Phase C: Lokale Automatisierung (ohne Workflow‑Manager)

- **Einfaches Orchestrieren**
  - [ ] Bash‑Skripte in `scripts/` erstellen: `00_predict_structure.sh`, `10_cg.sh`, `20_build_system.sh`, `30_equilibrate.sh`, `40_md.sh`, `50_features.sh`, `60_ml.sh`
  - [ ] Gemeinsame Konfigurationsdatei `config/local.env` (Seeds, Pfade, Parameter)
  - [ ] Logging in `reports/logs/` mit Exit‑Code‑Checks

Akzeptanzkriterium: Ein Einzeiler startet den kompletten Demo‑Lauf lokal: `./scripts/run_demo_local.sh <sample_id>`.

⸻

## Phase D: Stabilisierung des lokalen Vollbetriebs

- [ ] Replikate (2) lokal testen; Parameter fixieren (MDP, Komposition, Seeds)
- [ ] QC‑Checks skriptbasiert (Dicke, APL, Insertions‑Tiefe, Wasser im Core)
- [ ] Ergebnisstruktur standardisieren (Ordner/Benennung)

Akzeptanzkriterium: Wiederholbare lokale Läufe liefern robuste Features für 2–3 Peptide.

⸻

## Phase E: Dockerisierung (erst nach erfolgreichem lokalem E2E)

- **Container‑Plan**
  - [ ] Dockerfile(s) erstellen: `containers/gromacs-martini/`, `containers/martinize2/`, `containers/ml/`, optional `containers/af3/`
  - [ ] Versionierung `<tool>:<version>-<git-sha>`; Seeds/Commits in Artefakten speichern
  - [ ] GPU‑Support testen (NVIDIA Toolkit) für GROMACS/AF3

- **Abgleich mit lokalen Skripten**
  - [ ] Skripte so anpassen, dass sie wahlweise lokal ODER im Container laufen (Flag)
  - [ ] Volumes/Pfade vereinheitlichen; gleiche Outputs wie lokal

Akzeptanzkriterium: Gleicher Demo‑Lauf in Containern reproduziert lokale Artefakte (Hash/Prüfsummen wo sinnvoll).

⸻

## Phase F: Snakemake‑Automatisierung (auf Basis validierter Container)

- **Workflow**
  - [ ] `Snakefile` mit Regeln: `predict_structure`, `coarse_grain`, `build_membrane`, `equilibrate`, `md_run`, `preprocess_traj`, `extract_features`, `aggregate_features`, `train_model`, `explain_model`
  - [ ] `config/config.yaml` (Seeds, Pfade, Replikate, MDP‑Profile)
  - [ ] Lokales Snakemake‑Profil; Logging/Checkpointing (GROMACS `.cpt`)

Akzeptanzkriterium: `snakemake -n` korrekter DAG; Demo‑Run per Snakemake entspricht lokalem Ergebnis.

⸻

## Phase G: CI/CD, Skalierung und Cloud (optional, nach F)

- **CI**
  - [ ] Mini‑Demo (CPU‑only, 1 ns) im CI; Regressionstests (martinize2 Topologie‑Hash, Feature‑Ranges)
  - [ ] Pre‑commit Hooks; Linting der Snakefiles

- **Skalierung**
  - [ ] SLURM‑Profil (HPC) ODER GCP‑Profil (Life Sciences/K8s)
  - [ ] Speicherpolitik (Trajektorien‑Teilspeicherung), `results.sqlite` und Reports

Akzeptanzkriterium: CI grün; kleiner Batch auf HPC/Cloud reproduziert lokale Ergebnisse.

⸻

## Artefakte & Dateien (Stand lokal‑zuerst)

- **Konfiguration**: `config/local.env` (zunächst), später `config/config.yaml`
- **MDP**: `mdp/minim.mdp`, `mdp/nvt.mdp`, `mdp/npt.mdp`, `mdp/prod_demo.mdp`
- **Skripte**: `scripts/*` (Steps B1–B6), QC‑Skripte (`qc_membrane.py`, `extract_energies.py`)
- **Topologien/Lipide**: `top/` (Includes/Posres), `lipids/` (Kompositionslisten)

⸻

## Risiken & offene Punkte (lokal)

- [ ] AF3 lokal evtl. nicht verfügbar → Fallback/Ensemble definieren
- [ ] INSANE (Py2) lokal lauffähig sicherstellen; Alternativen prüfen
- [ ] Stabilität hoher Cholesterin‑Anteile in kleinen Patches; ggf. Parameter anpassen
- [ ] Rechenzeit lokal: Demo‑MD ggf. auf 20–50 ns begrenzen

⸻

## Roadmap (umgestellt auf lokal‑zuerst)

- **Woche 1**: Phase A + B (1 Peptid, kompletter manueller E2E, 20–50 ns)
- **Woche 2**: Phase C + D (Skripte, Replikate, stabile Outputs)
- **Woche 3**: Phase E (Dockerisierung, Reproduzierbarkeit ggü. lokal)
- **Woche 4**: Phase F (Snakemake‑Workflow), Basis‑CI; optional Phase G (HPC/Cloud)

Definition of Done (lokal‑zuerst): Vollständiger lokaler E2E‑Lauf für ≥1 Peptid, reproduzierbar per Skript, QC erfüllt; danach identisches Ergebnis in Containern und im Snakemake‑Run.
