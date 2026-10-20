#!/usr/bin/env python
"""
NP-BioSynAct -- drug-like lead filter. The raw AD-inside ranking rewards big
promiscuous NPs (median QED 0.07). Keep only drug-like, synthesizable, non-PAINS
multi-target leads. Reuses admet_leads.csv (already has all descriptors).

Run: conda activate npbs && python 34_final_leads.py
Out: final_leads.csv
"""
import os, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))

df = pd.read_csv(os.path.join(HERE, "admet_leads.csv"))
keep = df[(df.QED >= 0.5) & (df.Lipinski_viol <= 1) & (~df.PAINS) & (df.SAscore <= 6)]
keep = keep.sort_values(["n_hits", "mean_prob"], ascending=False)
keep.to_csv(os.path.join(HERE, "final_leads.csv"), index=False)
print(f"drug-like multi-target leads: {len(keep)} / {len(df)}")
print(keep.head(20)[["inchikey","bio_class","n_hits","mean_prob","QED","SAscore","BBB"]].to_string(index=False))
assert len(keep) > 0, "no drug-like leads survived -- loosen QED/SAscore cutoffs"
print("check ok -> final_leads.csv")
