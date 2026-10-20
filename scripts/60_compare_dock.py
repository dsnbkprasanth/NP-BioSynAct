#!/usr/bin/env python
"""smina vs gnina comparison + MD selection.
Reads dock_summary.csv (smina) + gnina_scores.csv (gnina CNN) + md_shortlist.csv.
Writes smina_vs_gnina.csv (per target x compound) and md_selection.csv (final MD set).
MD pick = the shortlisted best lead per target; flagged CONCORDANT when it is also
gnina's top-CNN compound at that target (both methods agree -> strongest evidence).
Run: python 60_compare_dock.py
"""
import os, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
def p(f): return os.path.join(HERE, f)

for f in ("dock_summary.csv", "gnina_scores.csv", "md_shortlist.csv"):
    if not os.path.exists(p(f)): sys.exit(f"!! missing {f} -- run docking/gnina first")

sm = pd.read_csv(p("dock_summary.csv"))
gn = pd.read_csv(p("gnina_scores.csv"))
sl = pd.read_csv(p("md_shortlist.csv"))

# --- normalise column names (dock_summary layout can vary) ---
def col(df, *cands):
    for c in cands:
        if c in df.columns: return c
    return None
sm_t = col(sm, "target"); sm_c = col(sm, "code", "lead", "compound", "best_lead")
sm_s = col(sm, "smina", "smina_score", "score", "best_lead_score", "affinity", "minimizedAffinity")
if not (sm_t and sm_c and sm_s): sys.exit(f"dock_summary cols unclear: {list(sm.columns)}")
sm = sm[[sm_t, sm_c, sm_s]].rename(columns={sm_t:"target", sm_c:"code", sm_s:"smina_kcal"})

gn = gn.rename(columns={c: c for c in gn.columns})
# gnina: lower vina = better; higher CNNaffinity = better
merged = sm.merge(gn, on=["target","code"], how="outer")

# per-target rank: smina (ascending kcal) vs gnina CNNaffinity (descending)
merged["smina_rank"] = merged.groupby("target")["smina_kcal"].rank(method="min")
if "CNNaffinity" in merged:
    merged["gnina_rank"] = merged.groupby("target")["CNNaffinity"].rank(method="min", ascending=False)
merged = merged.sort_values(["target","smina_rank"])
merged.to_csv(p("smina_vs_gnina.csv"), index=False)

# --- MD selection: shortlisted lead per target + concordance flag ---
gtop = {}  # target -> gnina best code by CNNaffinity
if "CNNaffinity" in gn:
    for t, g in gn.groupby("target"):
        gtop[t] = g.sort_values("CNNaffinity", ascending=False).iloc[0]["code"]

rows = []
for _, r in sl.iterrows():
    t, code = r["target"], r["best_lead"]
    m = merged[(merged.target==t) & (merged.code==code)]
    g = m.iloc[0] if len(m) else {}
    rows.append({
        "target": t, "code": code, "name": r.get("best_lead_name", code),
        "smina_kcal": r.get("best_lead_score", g.get("smina_kcal","")),
        "gnina_CNNaffinity": g.get("CNNaffinity",""),
        "gnina_CNNscore": g.get("CNNscore",""),
        "gnina_vina": g.get("vina",""),
        "gnina_top_at_target": gtop.get(t,""),
        "concordant": "YES" if gtop.get(t)==code else "no",
    })
md = pd.DataFrame(rows)
md.to_csv(p("md_selection.csv"), index=False)

print("=== smina vs gnina (MD set) ===")
print(md.to_string(index=False))
print(f"\nconcordant (both methods agree): {(md.concordant=='YES').sum()}/{len(md)}")
print("\nwrote: smina_vs_gnina.csv  md_selection.csv")
print("MD these 5 targets (best lead each). Concordant ones = strongest prior.")
assert len(md) == len(sl), "row count mismatch"
