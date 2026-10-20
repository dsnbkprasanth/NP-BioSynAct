# NP-BioSynAct computational workflow

NP-BioSynAct is organized as two complementary computational arms:

1. **Arm A:** natural-product curation, activity modelling, biosynthesis × target coupling, applicability-domain filtering, lead triage, and uncertainty/calibration analysis.
2. **Arm B:** docking, redocking validation, molecular dynamics, and end-point free-energy analysis for a focused set of target–lead pairs.

The two arms provide complementary evidence for prioritization. Agreement between computational methods should be interpreted as **orthogonal computational support**, not experimental validation.

## Arm A — ligand-based analysis

| Step | Script | Main output | Purpose |
|---|---|---|---|
| 11 | `11_parse_np_sdf.py` | `np_compounds.csv` | Merge and deduplicate natural-product structures and biosynthetic annotations. |
| 12 | `12_get_chembl_api.py` | `chembl_activities.csv` | Retrieve curated human single-protein activity records from ChEMBL. |
| 13 | `13_build_dataset.py` | `training_matrix.csv` | Construct active/inactive/unknown target labels. |
| 20 | `20_train_multitask.py` | `np_predictions.csv` | Train the multitask molecular-fingerprint model and score the NP library. |
| 50 | `50_cv.py` | `cv_results.csv` | Repeated scaffold-separated cross-validation. |
| 60 | `60_splits.py` | `split_comparison.csv` | Compare random and scaffold-separated evaluation. |
| 62 | `62_permutation.py` | permutation outputs | Estimate label and coupling null distributions. |
| 63 | `63_calibration.py` | `calibration.csv`, figure | Quantify calibration and MC-dropout uncertainty. |
| 30/33 | `30_coupling.py`, `33_atlas.py` | coupling tables | Construct the biosynthesis × target atlas. |
| 32 | `32_coupling_stats.py` | `coupling_stats.csv` | Effect sizes, multiple-testing correction, and scaffold controls. |
| 51 | `51_admet.py` | `admet_leads.csv` | Drug-likeness and developability filters. |
| 34 | `34_final_leads.py` | `final_leads.csv` | Produce the focused lead set. |
| 64 | `64_diversity_decoys.py` | diverse leads and decoys | Diversity selection and matched negative controls. |
| 31 | `31_validate.py` | recovery results | Retrospective recovery of known natural-product activities. |
| 52/53 | `52_biotransform.py`, `53_biotransform_control.py` | biotransformation tables | Compare parent and predicted-metabolite activity profiles. |
| 65 | `65_evidence_matrix.py` | evidence matrix | Integrate computational evidence across the selected systems. |
| 40 | `40_figures.py` | PNG figures | Generate released analysis figures. |

## Arm B — structure-based analysis

### Docking

Use `54_dock.sh` for smina docking and GNINA rescoring. `66_dock_prep.py` prepares coded ligands and target–ligand pairs.

Crystal-ligand redocking is used to check whether the selected docking setup reproduces a known pose. The released RMSD threshold is 2 Å.

### Molecular dynamics

Use `56_md.sh` with a compatible GROMACS installation. The released configuration specifies:

- AMBER99SB-ILDN for protein
- GAFF2 for ligands
- TIP3P water
- 0.15 M NaCl
- 100 ps NVT
- 100 ps NPT
- 100 ns production
- 2 fs integration step
- PME electrostatics
- 310 K and 1 bar

### End-point free energy

- `61_unigbsa.sh` provides single-pose Uni-GBSA estimates.
- `57_mmpbsa_analyze.sh` supports ensemble gmx_MMPBSA analysis.

MM/GBSA values are approximate end-point estimates and should not be interpreted as experimental binding free energies.

## Outputs

Released frozen results are in `data/frozen/`. Intermediate calculations generated during reproduction remain separate from the frozen release.

See `metadata/DATA_DICTIONARY.md` for field definitions and `docs/REPRODUCIBILITY.md` for external inputs and numerical reproducibility guidance.
