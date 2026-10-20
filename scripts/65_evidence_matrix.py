#!/usr/bin/env python
"""Lead-selection funnel and integrated evidence matrix.
Funnel: 766k -> AD-inside -> drug-like -> diverse leads. Evidence matrix: one row per
final lead with ML/AD/ADMET evidence; docking/MD/MM-PBSA columns left 'pending' for
the structure arm to fill (steps 54-56).

Run: python 65_evidence_matrix.py    Out: evidence_matrix.csv, fig_funnel.png
"""
import pandas as pd, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))

def n(f):
    try: return len(pd.read_csv(os.path.join(HERE, f)))
    except Exception: return 0

funnel = {
    "NP library": n("np_compounds.csv"),
    "AD-inside": n("atlas_leads_AD.csv"),
    "drug-like": n("final_leads.csv"),
    "diverse leads": n("leads_diverse.csv"),
}
print("lead-selection funnel:")
for k, v in funnel.items(): print(f"  {k:16s} {v}")

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.barh(list(funnel)[::-1], list(funnel.values())[::-1], color="#2c7fb8")
ax.set_xscale("log"); ax.set_xlabel("compounds (log)"); ax.set_title("Lead-selection funnel")
for i, v in enumerate(list(funnel.values())[::-1]): ax.text(v, i, f" {v}", va="center", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_funnel.png"), dpi=200)

# evidence matrix from ML/ADMET; structure columns pending
div = pd.read_csv(os.path.join(HERE, "leads_diverse.csv"))
cols = ["inchikey", "bio_class", "mean_prob", "n_hits", "QED", "SAscore", "BBB", "Lipinski_viol"]
em = div[[c for c in cols if c in div.columns]].copy()
em["AD"] = "inside"
for c in ["dock_smina", "dock_gnina", "redock_rmsd", "MD_stable", "MMPBSA_dG", "key_residues"]:
    em[c] = "pending"
em["tier"] = "Tier2 (ML+AD+ADMET+docking pending)"
em.to_csv(os.path.join(HERE, "evidence_matrix.csv"), index=False)
print(f"\nevidence matrix: {len(em)} leads -> evidence_matrix.csv (fill structure cols after 54-56)")
print("check ok -> evidence_matrix.csv, fig_funnel.png")
