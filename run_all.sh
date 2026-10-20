#!/usr/bin/env bash
# NP-BioSynAct -- full in-silico pipeline, one shot. Idempotent: heavy steps that
# already produced output are skipped. Re-run any single step by deleting its output.
set -euo pipefail
cd "$(dirname "$0")"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True   # avoid 5070 VRAM fragmentation
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate npbs

echo "== data =="
[ -f scripts/np_compounds.csv  ] || python scripts/11_parse_np_sdf.py
[ -f scripts/chembl_activities.csv  ] || python scripts/12_get_chembl_api.py
[ -f scripts/training_matrix.csv  ] || python scripts/13_build_dataset.py

echo "== model + predictions =="
[ -f scripts/np_predictions.csv  ] || python scripts/20_train_multitask.py

echo "== rigor: repeated scaffold CV =="
python scripts/50_cv.py

echo "== coupling atlas (honest) =="
python scripts/30_coupling.py
python scripts/33_atlas.py
python scripts/32_coupling_stats.py

echo "== leads: ADMET -> drug-like filter =="
python scripts/51_admet.py
python scripts/34_final_leads.py

echo "== validation: retrospective recovery =="
python scripts/31_validate.py

echo "== gut biotransformation =="
pip install --no-build-isolation sygma >/dev/null 2>&1 || echo "  (sygma install failed; skipping 52)"
python scripts/52_biotransform.py || echo "  (52 skipped)"

echo "== robustness checks and controls =="
python scripts/60_splits.py
python scripts/62_permutation.py
python scripts/63_calibration.py
python scripts/64_diversity_decoys.py
python scripts/65_evidence_matrix.py

echo "== figures =="
python scripts/40_figures.py

echo "== structure arm is manual: fill targets.tsv then 'bash 54_dock.sh', then MD/MM-PBSA (WORKFLOW.md) =="

echo "== ALL DONE =="
