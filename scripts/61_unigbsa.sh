#!/usr/bin/env bash
# Uni-GBSA prescreen on the 5 MD-shortlisted complexes.
# Env: gbsa. Runs from INSIDE each MD dir (proven to work manually).
# Uses whatever receptor.pdb is in the dir (user-repaired ones stay untouched).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

DIRS=(
  "MD/DPP4_NP06_6-6-9-Trimethyl-3-pentyl"
  "MD/MAPK14_NP08"
  "MD/NTRK1_NP07"
  "MD/OPRM1_NP08"
  "MD/PTGS2_NP08"
)
NT=${NT:-1}
mkdir -p unigbsa

for d in "${DIRS[@]}"; do
  T=$(basename "$d" | cut -d_ -f1-2)
  OUT="unigbsa/${T}.csv"
  [ -f "$OUT" ] && { echo "-- $T done, skip"; continue; }
  [ -f "$d/receptor.pdb" ] && [ -f "$d/ligand.mol" ] || { echo "!! $d missing files, skip"; continue; }
  # ligand: pH 7.4 (fixes odd-electron sqm crash)
  [ -f "$d/ligand.sdf" ] || obabel "$d/ligand.mol" -O "$d/ligand.sdf" -p 7.4 >/dev/null 2>&1
  echo "== $T =="
  ( cd "$d" && unigbsa-pipeline -i receptor.pdb -l ligand.sdf -o "${T}_out.csv" -nt "$NT" --verbose ) \
      || { echo "  !! $T failed"; continue; }
  # copy result up to unigbsa/
  [ -f "$d/${T}_out.csv" ] && cp "$d/${T}_out.csv" "$OUT"
done

# --- collate ---
python - <<'PY'
import glob, os, pandas as pd
rows=[]
for f in sorted(glob.glob("unigbsa/*.csv")):
    T=os.path.basename(f)[:-4]
    try: df=pd.read_csv(f); df.insert(0,"target_code",T); rows.append(df)
    except Exception as e: print(f"skip {f}: {e}")
if rows:
    out=pd.concat(rows, ignore_index=True); out.to_csv("unigbsa_scores.csv", index=False)
    print("wrote unigbsa_scores.csv"); print(out[["target_code","mode","TOTAL","status"]].to_string(index=False))
else:
    print("no results yet")
PY
