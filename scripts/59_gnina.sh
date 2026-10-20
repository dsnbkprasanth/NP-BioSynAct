#!/usr/bin/env bash
# gnina CNN rescoring on the 5 redock-valid targets, then push gnina's best pose
# into MD/<T>_<code>/ligand.mol (replacing the smina pose). Reuses rec/ from 54_dock.sh.
# Needs: gnina, rdkit. Runs on GPU automatically.
set -euo pipefail
cd "$(dirname "$0")"
TARGETS="DPP4 MAPK14 NTRK1 OPRM1 PTGS2"

for T in $TARGETS; do
  echo "== gnina $T =="
  gnina -r rec/${T}_rec.pdbqt -l ligands.sdf --autobox_ligand rec/${T}_ref.sdf \
        --num_modes 9 --seed 0 -o dock/${T}_gnina.sdf     # default CNN scoring
done

python - <<'PY'
import csv, os, re
from rdkit import Chem
def slug(s): return re.sub(r"[^A-Za-z0-9]+","-",str(s))[:24].strip("-")
short={r["target"]:(r["best_lead"],r["best_lead_name"]) for r in csv.DictReader(open("md_shortlist.csv"))}
rows=[]
for T in "DPP4 MAPK14 NTRK1 OPRM1 PTGS2".split():
    sdf=f"dock/{T}_gnina.sdf"
    if not os.path.exists(sdf): print("!! missing",sdf); continue
    best={}   # code -> first (best) mol
    for m in Chem.SDMolSupplier(sdf,sanitize=False):
        if m is None: continue
        nm=m.GetProp("_Name") if m.HasProp("_Name") else "?"
        if nm not in best: best[nm]=m
    for nm,m in best.items():
        g=lambda p: m.GetProp(p) if m.HasProp(p) else ""
        rows.append({"target":T,"code":nm,"CNNaffinity":g("CNNaffinity"),
                     "CNNscore":g("CNNscore"),"vina":g("minimizedAffinity")})
    # update MD ligand.mol with gnina best pose of the shortlisted NP
    if T in short:
        code,name=short[T]; d=f"MD/{T}_{code}"
        if code in best and os.path.isdir(d):
            Chem.MolToMolFile(best[code], f"{d}/ligand.mol")
            print(f"  {T}: MD ligand.mol updated to gnina pose ({code})")
csv.DictWriter(open("gnina_scores.csv","w"),fieldnames=["target","code","CNNaffinity","CNNscore","vina"]).writerows([{"target":"target","code":"code","CNNaffinity":"CNNaffinity","CNNscore":"CNNscore","vina":"vina"}]) if False else None
import pandas as pd; pd.DataFrame(rows).to_csv("gnina_scores.csv",index=False)
print("\ngnina_scores.csv written")
print(pd.DataFrame(rows)[pd.DataFrame(rows).code.isin([c for c,_ in short.values()])].to_string(index=False))
PY
echo "done -> gnina_scores.csv + MD/*/ligand.mol (gnina best pose). Now run MD in StreaMD."
