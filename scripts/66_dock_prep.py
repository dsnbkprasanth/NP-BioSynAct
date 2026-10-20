#!/usr/bin/env python
"""Dock prep + human-readable naming. Assigns clean codes NP01.. to the diverse
leads (the public-facing compound name), builds a master key table, and writes
ligands.sdf (3D, _Name=code) + pairs.csv using codes. Channels + NFKB1 excluded.

Run: python 66_dock_prep.py
Out: compound_key.csv (code<->id<->name<->class<->SMILES), ligands.sdf, pairs.csv
"""
import os, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
try: import requests
except Exception: requests = None
HERE = os.path.dirname(os.path.abspath(__file__))
DOCKABLE = ["OPRM1","CNR2","PPARG","DPP4","PTGS2","GSK3B","MAPK14","PTPN1","AKR1B1",
            "TNF","ADRB2","NTRK1","SIRT1","NOS2","ADORA1","HTR1A"]
TOPK = 3

lead = pd.read_csv(os.path.join(HERE, "leads_diverse.csv"))
meta = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                   usecols=["inchikey", "id", "source", "bio_class"], low_memory=False
                   ).drop_duplicates("inchikey").set_index("inchikey")
lead = lead.merge(meta, left_on="inchikey", right_index=True, how="left", suffixes=("", "_m"))
lead = lead.sort_values("mean_prob", ascending=False).reset_index(drop=True)
lead["code"] = [f"NP{i+1:02d}" for i in range(len(lead))]

def pubchem_name(ik):
    if requests is None: return None
    try:
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{ik}/property/Title/JSON", timeout=15)
        if r.ok: return r.json()["PropertyTable"]["Properties"][0].get("Title")
    except Exception: return None
    return None
lead["pubchem_name"] = [pubchem_name(ik) for ik in lead.inchikey]
lead["name"] = lead.pubchem_name.fillna(lead.id).fillna(lead.code)   # best readable label

key = lead[["code", "name", "pubchem_name", "id", "source", "bio_class", "inchikey", "smiles"]]
key.to_csv(os.path.join(HERE, "compound_key.csv"), index=False)

# pairs: code x top-K dockable targets
pred = pd.read_csv(os.path.join(HERE, "np_predictions.csv")).set_index("inchikey")
have = [t for t in DOCKABLE if t in pred.columns]
rows = []
for _, r in lead.iterrows():
    if r.inchikey not in pred.index: continue
    top = pred.loc[r.inchikey, have].sort_values(ascending=False).head(TOPK)
    for t, p in top.items():
        if p > 0.5:
            rows.append({"code": r.code, "name": r["name"], "target": t,
                         "prob": round(float(p), 3), "inchikey": r.inchikey, "smiles": r.smiles})
pairs = pd.DataFrame(rows).sort_values(["target", "prob"], ascending=[True, False])
pairs.to_csv(os.path.join(HERE, "pairs.csv"), index=False)

# ligands.sdf named by code
w = Chem.SDWriter(os.path.join(HERE, "ligands.sdf")); n = 0
for _, r in lead.iterrows():
    m = Chem.MolFromSmiles(r.smiles)
    if not m: continue
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) != 0: continue
    try: AllChem.MMFFOptimizeMolecule(m)
    except Exception: pass
    m.SetProp("_Name", r.code)                       # readable code, not InChIKey
    for c in ("name", "bio_class"): m.SetProp(c, str(r[c]))
    w.write(m); n += 1
w.close()
print(f"compound_key.csv: {len(key)} leads coded NP01..NP{len(key):02d}")
print(key[["code", "name", "bio_class"]].to_string(index=False))
print(f"\nligands.sdf: {n} 3D leads (named by code) | pairs.csv: {len(pairs)} pairs")
print(pairs[["code", "name", "target", "prob"]].to_string(index=False))
assert n > 0 and len(pairs) > 0, "no ligands/pairs"
print("check ok -> compound_key.csv, ligands.sdf, pairs.csv")
