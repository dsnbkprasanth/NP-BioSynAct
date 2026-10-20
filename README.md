# NP-BioSynAct

**NP-BioSynAct** is a reproducible two-arm computational workflow for large-scale natural-product polypharmacology analysis. It combines biosynthesis-aware ligand-based prediction with structure-based prioritization using docking, molecular dynamics, and end-point free-energy estimation.

The release contains executable analysis code, frozen derived datasets, figures, configuration, provenance records, and reproducibility metadata.

## What the pipeline does

### Arm A — biosynthesis-aware polypharmacology

- Curates natural products from COCONUT, LOTUS, and IMPPAT.
- Carries NPClassifier biosynthetic pathway, superclass, and class annotations.
- Builds a multi-target activity matrix from ChEMBL.
- Trains a multitask molecular-fingerprint MLP.
- Evaluates generalization with scaffold-separated cross-validation.
- Reports a 1-nearest-neighbour similarity baseline, applicability-domain status, calibration, uncertainty, and label-permutation controls.
- Quantifies biosynthetic-class × target coupling with effect sizes and multiple-testing correction.
- Performs ADMET/drug-likeness filtering, diversity selection, matched-decoy construction, and biotransformation analysis.

### Arm B — structure-based prioritization

- Prepares selected target–natural-product pairs.
- Performs smina docking and GNINA rescoring.
- Validates docking protocols by crystal-ligand redocking.
- Supports GROMACS molecular dynamics.
- Supports Uni-GBSA and gmx_MMPBSA analyses.
- Produces a compact evidence matrix combining ligand-based and structure-based readouts.

**Important:** docking, molecular dynamics, and MM/GBSA are computational prioritization methods. They do not constitute experimental validation or proof of biological activity.

## Release snapshot

- Natural products screened: **766,045**
- Activity targets in the multitask panel: **29**
- Scaffold-split macro ROC-AUC: **~0.886** in the released training analysis
- Random-split macro ROC-AUC: **0.935**
- Scaffold-split comparison: **0.883** in the independent split-comparison analysis
- Final structure-based panel: **5 target–lead systems**
- Hardware used for the released computational workflow: Linux x86_64, 32 GB RAM, NVIDIA RTX 5070 12 GB

The exact values associated with each analysis are preserved in `data/frozen/` and described in `metadata/DATA_DICTIONARY.md`.

## Repository layout

```text
NP-BioSynAct/
├── README.md
├── LICENSE
├── CITATION.cff
├── CHANGELOG.md
├── environment.yml
├── requirements.txt
├── run_all.sh
├── npbs_common.py
├── configs/
│   └── release_config.yaml
├── data/
│   └── frozen/
├── figures/
├── scripts/
├── docs/
│   ├── WORKFLOW.md
│   ├── REPRODUCIBILITY.md
│   └── references.ris
└── metadata/
    ├── DATA_DICTIONARY.md
    ├── DATA_MANIFEST.csv
    ├── DATA_MANIFEST.json
    ├── REPRODUCIBILITY_MANIFEST.json
    ├── SOURCE_PROVENANCE.md
    ├── PUBLIC_RELEASE_CHECKLIST.md
    └── CHECKSUMS.sha256
```

## Reproduction

### 1. Create the environment

```bash
conda env create -f environment.yml
conda activate npbs
```

External programs used by the structure arm are installed separately:

- GROMACS
- smina
- GNINA
- acpype
- Uni-GBSA
- gmx_MMPBSA

Version information is recorded in `metadata/REPRODUCIBILITY_MANIFEST.json`.

### 2. Obtain external datasets

Large source files are intentionally not stored in Git. See:

```text
metadata/SOURCE_PROVENANCE.md
docs/REPRODUCIBILITY.md
```

Expected source files include the COCONUT, LOTUS, and IMPPAT SDF snapshots and the ChEMBL activity export.

### 3. Run the ligand-based workflow

```bash
bash run_all.sh
```

The driver executes the reproducible ligand-based analyses and generates intermediate files under `scripts/`. Released frozen results are under `data/frozen/`.

### 4. Run the structure-based workflow

Prepare the target structures and selected ligands, then follow the commands in:

```text
docs/WORKFLOW.md
```

The structure arm requires external molecular-dynamics and docking software and is intentionally separated from the lightweight ligand-based workflow.

## Reproducibility and auditability

- Fixed random seeds are recorded in `configs/release_config.yaml`.
- Source versions and retrieval dates are recorded in `metadata/SOURCE_PROVENANCE.md`.
- Frozen derived data are accompanied by SHA-256 checksums.
- The data dictionary documents the released CSV/TSV fields.
- Numerical tolerances and tested software versions are recorded in the reproducibility manifest.
- Large source datasets and raw trajectories are excluded from Git through `.gitignore`.

Reproducibility is expected to be numerically equivalent rather than bitwise identical across different CPU/GPU and software stacks.

## Data and software licensing

The repository code is released under the MIT License. External databases retain their respective licenses and terms. See `metadata/SOURCE_PROVENANCE.md` before redistributing source datasets.

## Citation

Use `CITATION.cff` for software citation metadata. The repository is intended to provide the computational implementation and released derived data associated with the NP-BioSynAct study.
