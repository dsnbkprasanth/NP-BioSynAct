#!/usr/bin/env python
"""Analyse dock_scores.csv: redock pass/fail, best lead vs best decoy per target,
lead-vs-decoy plot, and pick MD complexes (best lead at redock-VALID targets only).

Run: python 55_dock_analyze.py    Out: dock_summary.csv, md_shortlist.csv, fig_dock.png
"""
import pandas as pd, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(HERE, "dock_scores.csv"))
key = pd.read_csv(os.path.join(HERE, "compound_key.csv")).set_index("code").name.to_dict()

d["score"] = pd.to_numeric(d.score, errors="coerce")
redock = d[d.ligand_set == "REDOCK"].set_index("target").score      # RMSD here
lig = d[(d.ligand_set == "ligands") & (d.engine == "smina")]
dec = d[(d.ligand_set == "decoys") & (d.engine == "smina")]

rows = []
for t in sorted(lig.target.unique()):
    rmsd = redock.get(t, float("nan"))
    L = lig[lig.target == t].score; D = dec[dec.target == t].score
    best = lig[lig.target == t].sort_values("score").iloc[0]
    rows.append({"target": t, "redock_rmsd": round(rmsd, 2),
                 "valid": rmsd <= 2.0,
                 "best_lead": best["name"], "best_lead_name": key.get(best["name"], best["name"]),
                 "best_lead_score": round(L.min(), 2), "mean_lead": round(L.mean(), 2),
                 "best_decoy_score": round(D.min(), 2), "mean_decoy": round(D.mean(), 2),
                 "lead_beats_decoy_mean": L.mean() < D.mean()})
s = pd.DataFrame(rows)
s.to_csv(os.path.join(HERE, "dock_summary.csv"), index=False)
print(s.to_string(index=False))

# MD shortlist: best lead at VALID targets only
md = s[s.valid][["target", "best_lead", "best_lead_name", "best_lead_score"]]
md.to_csv(os.path.join(HERE, "md_shortlist.csv"), index=False)
print(f"\nMD complexes (redock-valid targets): {len(md)}")
print(md.to_string(index=False))

# lead vs decoy plot (valid targets) -- guarded (numpy2/matplotlib clash in some envs)
try:
    valid = s[s.valid].target.tolist()
    fig, ax = plt.subplots(figsize=(8, 4))
    for i, t in enumerate(valid):
        Ls = [float(x) for x in lig[lig.target == t].score]
        Ds = [float(x) for x in dec[dec.target == t].score]
        ax.scatter([i-0.15]*len(Ls), Ls, c="#2c7fb8", s=18, label="lead" if i == 0 else "")
        ax.scatter([i+0.15]*len(Ds), Ds, c="#d95f02", s=18, label="decoy" if i == 0 else "")
    ax.set_xticks(range(len(valid))); ax.set_xticklabels(valid, rotation=45)
    ax.set_ylabel("smina affinity (kcal/mol)"); ax.legend()
    ax.set_title("Leads vs matched decoys (redock-valid targets)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_dock.png"), dpi=200)
    print("fig_dock.png written")
except Exception as e:
    print(f"(plot skipped: {e} -- run in npbs env for the figure; CSVs are done)")
print("check ok -> dock_summary.csv, md_shortlist.csv")
