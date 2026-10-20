# Reproducibility

## External data

Large source datasets are not redistributed in Git. Source names, versions, retrieval dates, provenance, and licensing information are recorded in `metadata/SOURCE_PROVENANCE.md`.

The release was built from:

- COCONUT 3D SDF, July 2026 snapshot
- LOTUS v9
- IMPPAT 2.0
- ChEMBL Release 34 activity records retrieved through the REST API

The source files are expected to be placed where the corresponding scripts can access them. Their expected filenames are documented in `scripts/11_parse_np_sdf.py` and `metadata/SOURCE_PROVENANCE.md`.

## Environment

The primary Python environment is defined in `environment.yml`. A minimal pip dependency list is provided in `requirements.txt`.

External structure-based tools are installed separately because their binary builds and GPU configurations vary:

- GROMACS
- smina
- GNINA
- acpype
- Uni-GBSA
- gmx_MMPBSA

Tested software versions are recorded in `metadata/REPRODUCIBILITY_MANIFEST.json`.

## Determinism

Fixed seeds are stored in `configs/release_config.yaml`.

Bitwise identity is not expected across different GPU architectures, BLAS implementations, or software versions. The release manifest records practical numerical tolerances for key outputs.

## Frozen results

`data/frozen/` contains the released derived tables used to describe the computational results. SHA-256 checksums are provided in `metadata/CHECKSUMS.sha256`.

`metadata/DATA_DICTIONARY.md` describes the fields in the frozen tables.

## Re-running

From the repository root:

```bash
conda env create -f environment.yml
conda activate npbs
bash run_all.sh
```

The structure arm is run separately after target structures and ligand parameters have been prepared; see `docs/WORKFLOW.md`.

## Scope

This repository is a computational research resource. It does not contain experimental validation data, patient data, clinical data, or proprietary source datasets.
