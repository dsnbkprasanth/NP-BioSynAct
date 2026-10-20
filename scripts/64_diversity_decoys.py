#!/usr/bin/env python
"""Chemical-diversity lead selection and matched decoys.
Butina-cluster the drug-like AD leads and take one representative per cluster (so the
structural arm tests DIVERSE predictions, not analogues). Build a property-matched
low-confidence decoy set so docking/MD can show the ML selection ENRICHES for binders.

Run: python 64_diversity_decoys.py    Out: leads_diverse.csv, decoys.csv
"""
import numpy as np, pandas as pd, os
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Descriptors, MolFromSmiles
from npbs_common import fpmol, HERE
N_SEL = 8

smi = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                  usecols=["inchikey", "smiles"], low_memory=False).set_index("inchikey").smiles
lead = pd.read_csv(os.path.join(HERE, "final_leads.csv")).head(100)
lead["smiles"] = lead.inchikey.map(smi)
lead = lead.dropna(subset=["smiles"]).reset_index(drop=True)

fps = [fpmol(s) for s in lead.smiles]
ok = [i for i, f in enumerate(fps) if f]; lead = lead.iloc[ok].reset_index(drop=True); fps = [fps[i] for i in ok]
# Butina clustering on Tanimoto distance
dists = []
for i in range(1, len(fps)):
    s = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i]); dists += [1 - x for x in s]
clusters = Butina.ClusterData(dists, len(fps), 0.4, isDistData=True)
reps = [c[0] for c in clusters][:N_SEL]                     # first = highest-ranked in cluster
div = lead.iloc[reps].copy()
div.to_csv(os.path.join(HERE, "leads_diverse.csv"), index=False)
print(f"{len(clusters)} chemical clusters -> {len(div)} diverse leads")
print(div[["inchikey", "bio_class", "n_hits", "mean_prob", "QED"]].to_string(index=False))

# matched decoys: AD-inside, low predicted activity, MW-matched to the leads
ad = pd.read_csv(os.path.join(HERE, "atlas_leads_AD.csv"))
ad["smiles"] = ad.inchikey.map(smi)
ad = ad.dropna(subset=["smiles"])
low = ad[(ad.n_hits == 0) & (ad.mean_prob < 0.2)].copy()
low["MW"] = [Descriptors.MolWt(MolFromSmiles(s)) if MolFromSmiles(s) else np.nan for s in low.smiles]
low = low.dropna(subset=["MW"])
div_mw = [Descriptors.MolWt(MolFromSmiles(s)) for s in div.smiles]
picks = []
for mw in div_mw:                                           # nearest-MW decoy per lead
    cand = low[~low.inchikey.isin(picks)]
    if cand.empty: break
    picks.append(cand.iloc[(cand.MW - mw).abs().argmin()].inchikey)
dec = low[low.inchikey.isin(picks)]
dec.to_csv(os.path.join(HERE, "decoys.csv"), index=False)
print(f"matched decoys: {len(dec)} (low prob, AD-inside, MW-matched)")
assert len(div) > 0 and len(dec) > 0, "diversity/decoy selection empty"
print("check ok -> leads_diverse.csv, decoys.csv")
