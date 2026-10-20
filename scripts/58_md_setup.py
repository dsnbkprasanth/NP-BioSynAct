#!/usr/bin/env python
"""Build MD/ tree for manual StreaMD runs. Per valid complex (md_shortlist.csv):
  MD/<TARGET>_<CODE>_<name>/ligand.mol   (best docked pose, +H, coords kept)
                            receptor.pdb  (cleaned receptor from docking)
Run: python 58_md_setup.py
"""
import os, re, shutil
import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

sl = pd.read_csv(os.path.join(HERE, "md_shortlist.csv"))
def slug(s): return re.sub(r"[^A-Za-z0-9]+", "-", str(s))[:24].strip("-")

n = 0
for _, r in sl.iterrows():
    T, code = r.target, r.best_lead
    name = r.get("best_lead_name", code)
    d = os.path.join(HERE, "MD", f"{T}_{code}_{slug(name)}"); os.makedirs(d, exist_ok=True)
    # receptor
    rec = os.path.join(HERE, "rec", f"{T}_rec.pdb")
    if os.path.exists(rec): shutil.copy(rec, os.path.join(d, "receptor.pdb"))
    else: print(f"  !! {rec} missing")
    # best pose of this NP from the smina docking sdf
    sdf = os.path.join(HERE, "dock", f"{T}_ligands_smina.sdf")
    lig = None
    for m in Chem.SDMolSupplier(sdf):
        if m and (m.GetProp("_Name") if m.HasProp("_Name") else "") == code:
            lig = m; break
    if lig is None:
        print(f"  !! pose {code} not found in {sdf}"); continue
    ligH = Chem.AddHs(lig, addCoords=True)                 # H with 3D coords
    Chem.MolToMolFile(ligH, os.path.join(d, "ligand.mol"))
    n += 1
    print(f"  {T}_{code}_{slug(name)}: ligand.mol + receptor.pdb")

print(f"\n{n} complexes ready under MD/  (run each in StreaMD)")
assert n > 0, "nothing built"
print("check ok")
