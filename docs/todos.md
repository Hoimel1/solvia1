## To‑Do-Liste: End‑to‑End‑Pipeline für MD‑basierte Toxizitätsvorhersage von AMPs (RBC, proteinfreie Membran)

Hinweis: Alle Aufgaben sind so formuliert, dass sie in Snakemake integriert, containerisiert und reproduzierbar sind. Verwende konsistente Seeds (42), speichere Container‑Digests/Commits mit, und halte Ordner‑/Dateinamen strikt am Schema `data/<stage>/<sample_id>/...`.

- **Globale Grundlagen**
  - [ ] Repository‑Struktur finalisieren: `data/{raw,interim,cg,systems,md,features,models,reports}`, `scripts/`, `mdp/`, `top/`, `lipids/`, `profiles/{local,slurm,gcp}`, `containers/`, `config/`
  - [ ] Namenskonventionen definieren: `sample_id`, Replikate (`_rep1`, `_rep2`), AF‑Seed (`_s{n}`), Stages
  - [ ] Projektkonfiguration anlegen: `config/config.yaml` (alle Flags/Parameter, s.u.)
  - [ ] Vorlagen für README/Runbook/DAG: `docs/usage.md`, `docs/runbook.md`, `snakemake --dag`
  - [ ] Pre‑commit + Code‑Style (Black, Ruff o.ä.), Lizenz, CITATION.cff

## Container & Umgebungen

- **Container-Registry & Basis**
  - [ ] Registry/Namespace festlegen (DockerHub/GHCR) und Zugänge konfigurieren
  - [ ] Gemeinsame Labels/Versionierung: `<tool>:<version>-<git-sha>`; Export der Image‑Digests in Artefakte

- **AlphaFold 3**
  - [ ] Multi‑Stage Dockerfile erstellen: OS, CUDA, PyTorch+JAX, AF3 Code (gepinnt), Entrypoint
  - [ ] NVIDIA Container Toolkit testen (GPU‑Support); Apptainer SIF via `docker2singularity.sh`
  - [ ] AF‑Datenbanken mountbar machen (Full + „lite“ Bundle); Auswahl per Flag `af_db: full|lite`
  - [ ] CI‑Smoke‑Test (CPU‑only Mini‑Peptid)

- **martinize2 v0.14.0 + Vermouth 0.14.x**
  - [ ] Leichtes Image bauen (Python, DSSP/mkdssp); Versionen pinnen
  - [ ] IDP‑Parameter und ElNeDyn‑Optionen enthalten; CLI‑Wrapper bereitstellen

- **GROMACS 2024 (Martini 3)**
  - [ ] GPU‑fähiges Image mit `gmx_mpi`, Python, MDAnalysis Tools
  - [ ] Martini 3.0.0 Dateien, `insane.py` (Py2) beilegen/zugänglich machen

- **ML/Analyse**
  - [ ] Python 3.10, `xgboost`, `scikit-learn`, `optuna`, `shap`, `mdanalysis`, `numpy/pandas`
  - [ ] Exportformate (XGBoost JSON, optional ONNX) testen

Akzeptanzkriterium: Alle Images baubar, versioniert, kurzer CI‑Lauf grün.

## Snakemake‑Workflow

- **Grundgerüst**
  - [ ] `Snakefile` mit Regeln: `ingest_data`, `predict_structure`, `coarse_grain`, `cg_qc`, `build_membrane`, `equilibrate`, `md_run`, `preprocess_traj`, `extract_features`, `aggregate_features`, `train_model`, `explain_model`, `reports`
  - [ ] `config/config.yaml`: Seeds, Pfade, `elastic_mode: adaptive|none`, `idp_mode: auto|off`, `membrane_type: planar|vesicle`, Cholesterin‑Anteil, Replikate, MDP‑Profile
  - [ ] Ressourcen/Clusterprofile: `profiles/local`, `profiles/slurm` (SBATCH‑Vorlage), `profiles/gcp` (Life Sciences/K8s)
  - [ ] Logging/Checkpointing: GROMACS `.cpt` alle 50 ns, Snakemake `checkpoint`s für Verlängerungen

Akzeptanzkriterium: `snakemake -n` zeigt korrekten DAG; Mini‑Demo läuft end‑to‑end (CPU‑only).

## Datenaufnahme & Metadaten

- **Eingangsdaten**
  - [ ] Peptidsequenzen + HC50 (CSV/TSV) nach `data/raw/peptides.csv`
  - [ ] Schema validieren (ID, Sequenz, Länge, HC50, Quelle); Pydantic/Schema‑Check
  - [ ] Sample‑Sheet mit Replikat‑Expansionslogik

- **Metadaten/Provenienz**
  - [ ] `cg_metadata.db`, `cg_qc.db`, zentrale `results.sqlite` (Features, Labels, Splits)

Akzeptanzkriterium: Konsistenzprüfungen grün; IDs eindeutig; Replikate korrekt generiert.

## Schritt 1: Strukturvorhersage (AF3)

- **Batching & Reproduzierbarkeit**
  - [ ] Batching kleiner Peptide (bis 32/16GB GPU) implementieren
  - [ ] 5 Seeds/Peptid; pLDDT/RMSD‑Clustering; Zentroid auswählen

- **Qualitätsfilter**
  - [ ] pLDDT‑Grenze (mean ≥ 80), MolProbity Clash‑Score ≤ 20 (falls verfügbar)
  - [ ] Vakuum‑Minimierung (AmberTools) zur Clash‑Auflösung
  - [ ] Ensemble‑Fallback (RoseTTAFold/helical builder) bei Grenzfällen; optional kurze Lösung‑MD

- **Outputs**
  - [ ] Speichere `ranked_*.pdb`, pLDDT‑Per‑Residue, Seed, Container‑Digest unter `data/interim/<id>/`

Akzeptanzkriterium: ≥95% Peptide bestehen QC; Metadaten vollständig; CI‑Beispiel ok.

## Schritt 2: Coarse‑Graining (martinize2, Martini 3, IDP)

- **Deterministische Ausführung**
  - [ ] martinize2 Aufruf fixieren (Seed 42, `-ff martini3001`, `-water martini3`, `-dssp mkdssp`)
  - [ ] Adaptive ElNeDyn: `elastic_mode: adaptive|none`; Skript `scripts/prune_elastic.py` (pLDDT < 70 aussparen)
  - [ ] IDP‑Umschaltung: Martini3‑IDP für disordered Peptide per `idp_mode: auto`

- **Artefakte**
  - [ ] `*_cg.pdb`, `*.top`, `*.itp`, ggf. `elastic.itp`, Checksummen → `data/cg/<id>/`

- **QC**
  - [ ] Kurze Vakuum‑Minimierung (100 Schritte, Fmax < 1000 kJ/mol/nm)
  - [ ] Backmapping‑RMSD (Backward 2025.2): Cα RMSD ≤ 0.45 nm
  - [ ] Statistik ElNeDyn (Anzahl/Längen) → `elastic_qc.tsv`; Persistenz in `cg_qc.db`
  - [ ] CI‑Regression: Hash‑Vergleich `*.top` für 2 Referenzpeptide

Akzeptanzkriterium: Alle CG‑Modelle bestehen QC; Hash‑CI stabil.

## Schritt 3: Membranbau (INSANE, RBC‑ähnlich, asymmetrisch)

- **Komposition & Asymmetrie**
  - [ ] Outer: überwiegend PC/SM‑Äquivalent + Chol (~45% gesamt)
  - [ ] Inner: PE‑reich (~50%), PS (~10–20%), Chol (~30–45%)
  - [ ] Leafflet‑Asymmetrie aktivieren (`-asym`), moderates Lipidzahl‑Ungleichgewicht (innen > außen, ~10–20%)

- **Systemaufbau**
  - [ ] `insane.py` Aufruf (Py2): Box ~10×10×Z nm, Salz 150 mM, `-chargeneutral`
  - [ ] Peptidplatzierung ~1 nm über äußerem Leaflet; Seeds für ≥2 Replikate/Peptid
  - [ ] Artefakte: `system.gro`, `system.top`, Index‑Dateien, Lipidzählungen → `data/systems/<id>/`

- **Equilibration**
  - [ ] Minimierung; 100 ps NVT (310 K, posres Peptid); 5–10 ns NPT (semi‑isotrop), posres lockern

- **QC**
  - [ ] Dicke (P‑P Abstand) ~4.5–5.0 nm; APL je Leaflet im plausiblen Bereich
  - [ ] Asymmetrie erhalten (PS nur innen), Chol‑Redistribution beobachten
  - [ ] Keine verbotenen Overlaps (gmx mindist)

Akzeptanzkriterium: Systeme physikalisch plausibel, keine Flip‑Flops außer Chol, stabile APL.

## Schritt 4: Produktions‑MD (0.5–1 µs, GROMACS 2024)

- **MDP‑Sätze in `mdp/` anlegen**
  - [ ] `minim.mdp`, `nvt.mdp`, `npt.mdp`, `prod.mdp` (dt=0.02 ps, cutoff 1.1 nm, RF ε_r=15)
  - [ ] Thermostat V‑Rescale (τ=1 ps), PR‑Barostat semi‑isotrop (τ=12 ps), LINCS‑Einstellungen
  - [ ] Energy‑Groups: `Protein` vs `Membrane` (Coul/LJ Interaktionen)

- **Läufe & Replikate**
  - [ ] 2 (optional 3) Replikate (unabhängige Velocities/INSANE‑Seeds) je Peptid, 500 ns (ggf. auf 1 µs verlängern)
  - [ ] Checkpoints alle 50 ns; Ressourcen: 1 GPU + 4–8 CPU Threads/Job

- **Monitoring**
  - [ ] Energie/Temp‑Stabilität, APL, Dicke, Insertions‑Tiefe, Wasser im Core, evtl. Pore‑Indikatoren

Akzeptanzkriterium: Läufe ohne Instabilitäten; Verlängerung nur bei Drift‑Indizien.

## Schritt 5: Trajektorien‑Vorverarbeitung

- **Bildgebung/Zentrierung**
  - [ ] `gmx trjconv`: `-pbc mol -center -ur compact`, Membran‑COM als Zentrum
  - [ ] Ausrichtung der Membran‑Normalen auf z (falls erforderlich)
  - [ ] Artefakte: `*_centered.xtc`, referenzierte `*.gro`

Akzeptanzkriterium: Peptid und Membran in derselben PBC‑Zelle, konsistente Z‑Ausrichtung.

## Schritt 5b: Feature‑Extraktion (letzte 100 ns)

- **Implementierung (`scripts/extract_features.py`)**
  - [ ] d_min (Peptid ↔ Phosphatbead), COM‑Distanz/Insertions‑Tiefe
  - [ ] Δ‑Dicke (nah Peptid vs fern), APL‑Änderung ggü. Kontrolle
  - [ ] Peptid RMSD, Radius der Gyration
  - [ ] Kontaktzahlen total/je Lipidtyp (PS/PC), Wasser‑Kontakte
  - [ ] Orientierungswinkel (Tilt) zur z‑Achse
  - [ ] Interaktionsenergien (Coul‑SR/LJ‑SR Protein‑Membrane aus `.edr`)
  - [ ] Wasser im Membrancore (z‑Fenster um Midplane)

- **Technik**
  - [ ] MDAnalysis/NumPy, optional Grid‑basierte Dicke (10×10)
  - [ ] Trajektorien‑Slicing (letzte 100 ns), Performance‑Optimierung (Chunking)
  - [ ] CSV/Parquet pro Replikat → `data/features/replicates/`, Validierung auf plausible Bereiche

Akzeptanzkriterium: Features für alle Replikate vorhanden, Einheiten/Definitionen dokumentiert.

## Schritt 6: Aggregation & ML

- **Aggregation**
  - [ ] Replikate mitteln; std optional als Zusatzfeature; Join mit HC50/Label
  - [ ] Z‑Score‑Normalisierung (train‑set‑fit, test‑set‑apply), log10(HC50) für Regression

- **Modelle (XGBoost)**
  - [ ] Classifier (toxisch ≤ 100 µM), Regressor (log10(HC50))
  - [ ] Class‑Imbalance: `scale_pos_weight`
  - [ ] Hyperparameter‑Tuning (Grid + Optuna, 5‑fold CV)
  - [ ] Metriken: ROC‑AUC, PR‑AUC, Recall toxisch, Spezifität, R²/RMSE
  - [ ] Persistenz: `models/xgb_classifier.json`, `models/xgb_regressor.json`; optional ONNX‑Export

Akzeptanzkriterium: Zielmetriken ~AUC ≥ 0.9 (Klass.), R² ~0.65–0.7 (Regr.) auf Test‑Split.

## Schritt 7: Interpretierbarkeit (SHAP)

- **Global & Lokal**
  - [ ] Mean |SHAP| Ranking, Summary Beeswarm
  - [ ] Dependence‑Plots mit Interaktionsfärbung (z.B. Insertions‑Tiefe × Coulomb‑Energie)
  - [ ] Force‑Plots für ausgewählte Peptide (TP/FP/FN)

- **Berichte**
  - [ ] Plot‑Sammlung → `reports/shap/`; Kurzinterpretation in `reports/model_report.md`

Akzeptanzkriterium: Kohärente Mechanismen (Tiefe/Δ‑Dicke/Energie) widerspiegeln sich in SHAP.

## Skalierung & Deployment

- **HPC (SLURM)**
  - [ ] Cluster‑Profil (`--cluster "sbatch ..."`), Ressourcenzuweisung je Regel (`gpus`, `threads`, `mem`, `time`)
  - [ ] Testlauf auf 32+ GPUs; Durchsatz messen

- **Cloud (GCP)**
  - [ ] Profil für Google Life Sciences/Kubernetes ODER feste GPU‑VM‑Pools (Terraform optional)
  - [ ] Artefakt‑Upload in Bucket, Preemptible‑Strategie, Wiederaufnahme via Checkpoints

- **Ergebnismanagement**
  - [ ] Speicherpolitik (Trajektorien‑Teilspeicherung: letzte 100 ns; Aufräum‑Regeln)
  - [ ] Zentrales `results.sqlite` + CSV‑Exports; Versionierte Reports

Akzeptanzkriterium: 1500 Peptide × 2 Replikate in ≤ 1 Woche auf 50–100 GPUs realistisch.

## CI/CD & Tests

- **Mini‑E2E**
  - [ ] 2 Demo‑Peptide, 1 ns Läufe, CPU‑only, vollständiger DAG in GitHub Actions

- **Regressions**
  - [ ] martinize2 Hash‑Vergleich; Feature‑Sanity‑Checks; Snakefile `--lint`

- **Qualität**
  - [ ] Pre‑commit Hooks; Static Checks; Repro‑Seeds in Artefakt‑Metadaten

Akzeptanzkriterium: CI grün, deterministische Outputs (innerhalb Toleranz).

## Dokumentation

- [ ] Nutzerleitfaden (`docs/usage.md`): Voraussetzungen, Profile, Starten/Stoppen
- [ ] `docs/concept.md` Querverweise, Parametertabellen, Standardwerte
- [ ] `docs/runbook.md`: Fehlersuche (z.B. Instabilitäten, Chol‑Verteilung, Pore‑Events)
- [ ] DAG‑Graph, Metrik‑Dashboard (Plotly/HTML‑Report), `snakemake --report`

## Artefakte & Dateien anzulegen/zu pflegen

- **Konfiguration**
  - [ ] `config/config.yaml`, `profiles/{local,slurm,gcp}/config.yaml`
- **MDP**
  - [ ] `mdp/minim.mdp`, `mdp/nvt.mdp`, `mdp/npt.mdp`, `mdp/prod.mdp`
- **Skripte**
  - [ ] `scripts/prune_elastic.py`, `scripts/extract_features.py`, `scripts/extract_energies.py`, `scripts/qc_membrane.py`
- **Container**
  - [ ] `containers/af3/Dockerfile`, `containers/martinize2/Dockerfile`, `containers/gromacs-martini/Dockerfile`, `containers/ml/Dockerfile`, `scripts/docker2singularity.sh`
- **Topologien/Lipide**
  - [ ] `top/` (Peptid‑Includes/Position‑Restraints), `lipids/` (Kompositions‑Listen, INSANE‑Vorlagen)

## Risiken & offene Punkte

- [ ] AF3 Lizenz/Redistribution der Gewichte und DBs (nur Mount?)
- [ ] INSANE (Py2) in Containerlaufzeit sicherstellen; Alternativen testen
- [ ] Martini3‑IDP Verfügbarkeit/Version exakt pinnen
- [ ] Hoher Cholesterin‑Anteil: Stabilität der Patch‑Asymmetrie (ggf. Vesikel‑Modus testen)
- [ ] Speicherbedarf/IO bei 3000 Läufen (Chunking, Kompression, Parquet)

## Roadmap (grobe Meilensteine)

- **M1 (Woche 1–2)**: Container + Snakemake‑Gerüst + Mini‑E2E (CPU)
- **M2 (Woche 3–4)**: AF3 + martinize2 + Membranbau + Equilibration stabil
- **M3 (Woche 5–6)**: Produktions‑MD + Feature‑Extraktion + Aggregation vollständig
- **M4 (Woche 7)**: ML + SHAP + Berichte; CI Vollabdeckung
- **M5 (Woche 8)**: Skalierung HPC/Cloud, Kosten/Throughput bestätigt

Definition of Done: Alle Akzeptanzkriterien pro Abschnitt erfüllt; `snakemake --report` liefert vollständigen, reproduzierbaren Lauf inkl. Modelle, SHAP‑Berichte und QC‑Dashboards.
