#!/usr/bin/env python
"""Biotransform negative control. The 96% "gain" may be size/decoration bias
(bigger metabolite -> higher max-Tanimoto to ANY set). Re-score the same
parent/metabolite pairs against RANDOM NP reference pools (not real actives).
 - real gain% (from biotransform.csv) vs random-reference gain%.
 - if similar -> inflation, reframe. If real >> random -> genuine signal.

Run: python 53_biotransform_control.py    Out: biotransform_control.csv
"""
import numpy as np, pandas as pd, os
from rdkit import DataStructs
from npbs_common import fpmol, HERE, RELIABLE

bt = pd.read_csv(os.path.join(HERE, "biotransform.csv"))
smi = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                  usecols=["inchikey", "smiles"], low_memory=False).set_index("inchikey").smiles
bt["parent_smiles"] = bt.inchikey.map(smi)

# random reference: 20 pools of 200 random NPs (mirrors 20-target scorer shape)
rng = np.random.default_rng(0)
pool_smi = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                       usecols=["smiles"], low_memory=False).smiles.dropna()
pools = [[f for f in (fpmol(s) for s in pool_smi.sample(200, random_state=i)) if f]
         for i in range(len(RELIABLE))]

def polypharm_rand(s):
    q = fpmol(s)
    if not q: return None
    return sum(max(DataStructs.BulkTanimotoSimilarity(q, p)) for p in pools)

rows = []
for _, r in bt.iterrows():
    if not isinstance(r.parent_smiles, str) or not isinstance(r.best_metabolite, str) or not r.best_metabolite:
        continue
    pp, mp = polypharm_rand(r.parent_smiles), polypharm_rand(r.best_metabolite)
    if pp is None or mp is None: continue
    rows.append({"inchikey": r.inchikey, "gain_real": r.gain, "gain_random": mp - pp})
out = pd.DataFrame(rows); out.to_csv(os.path.join(HERE, "biotransform_control.csv"), index=False)

real_pct = (bt.gain > 0).mean() * 100
rand_pct = (out.gain_random > 0).mean() * 100
print(f"metabolites gaining vs REAL actives : {real_pct:.0f}%  (mean gain {bt.gain.mean():.3f})")
print(f"metabolites gaining vs RANDOM pools : {rand_pct:.0f}%  (mean gain {out.gain_random.mean():.3f})")
delta = real_pct - rand_pct
print(f"delta = {delta:.0f} pp")
if delta < 15:
    print("VERDICT: mostly structural inflation -> reframe as 'metabolites stay in active space', not 'metabolism activates'")
else:
    print("VERDICT: real activity-specific gain beyond size bias -> claim holds")
assert len(out) > 0, "no pairs scored"
print("check ok -> biotransform_control.csv")
