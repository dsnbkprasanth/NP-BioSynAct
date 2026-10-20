#!/usr/bin/env python
"""
NP-BioSynAct -- HONEST reframe (path A). Produces the defensible atlas:
  * targets trimmed to model AUC in [0.80, 0.995] (beats 1-NN, not small-n perfect)
  * predictions restricted to applicability domain (maxTanimoto to ChEMBL >= 0.3)
  * couplings kept only if: q<0.05, effect>0.6, AND scaffold-dedup robust
    (median_dedup/median_full >= 0.7)  -> not a single analog series
  * biosynthetic classes significant for >50% of targets flagged PROMISCUOUS
    and excluded from the specificity claim (kills the oligopeptide catch-all)

Run: conda activate npbs && python 33_atlas.py
Out: atlas_couplings_specific.csv, atlas_couplings_all.csv, atlas_leads_AD.csv
"""
import os
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# model AUC (this run); keep beats-baseline band, drop 1.000 small-n + <0.80
AUC = {"P2RX3":.989,"TRPA1":.990,"HTR1A":.981,"TNF":.965,"ADRB2":.964,"NTRK1":.961,
       "OPRM1":.960,"PPARG":.958,"SCN9A":.954,"PTPN1":.952,"AKR1B1":.950,"CNR2":.944,
       "SIRT1":.935,"GSK3B":.933,"NOS2":.912,"ADORA1":.883,"SCN10A":.860,"MAPK14":.840,
       "DPP4":.839,"NFKB1":.819,"NLRP3":.791,"GRM5":.732,"PTGS2":.732,"TRPV1":.614,
       "IL6":1.,"RAGE":1.,"TLR4":1.,"IL1B":.349}
TARGETS = [t for t, a in AUC.items() if 0.80 <= a <= 0.995]
AD_CUT, EFF, QCUT, ROBUST = 0.30, 0.60, 0.05, 0.70


def bh(p):
    p = np.asarray(p); n = len(p); o = p.argsort(); q = np.empty(n)
    q[o] = p[o] * n / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def fp(s):
    m = Chem.MolFromSmiles(s); return FPGEN.GetFingerprint(m) if m else None


def run():
    print(f"targets kept: {len(TARGETS)}")
    pred = pd.read_csv(os.path.join(HERE, "np_predictions.csv")).set_index("inchikey")
    bio = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                      usecols=["inchikey", "smiles", "bio_class", "murcko_scaffold"],
                      low_memory=False).dropna(subset=["bio_class"]).set_index("inchikey")
    df = pred[TARGETS].join(bio, how="inner").dropna(subset=["smiles"])
    print(f"labelled NPs: {len(df)}")

    # ---- applicability domain: maxTanimoto to a ChEMBL reference sample -----
    ref = pd.read_csv(os.path.join(HERE, "chembl_activities.csv"),
                      usecols=["smiles"]).smiles.dropna().drop_duplicates().sample(
                      3000, random_state=0)
    ref_fp = [x for x in map(fp, ref) if x]
    ad = np.empty(len(df))
    for i, s in enumerate(df.smiles.values):
        q = fp(s)
        ad[i] = max(DataStructs.BulkTanimotoSimilarity(q, ref_fp)) if q else 0.0
    df["ad"] = ad
    ind = df[df.ad >= AD_CUT].copy()
    print(f"in-domain (Tanimoto>={AD_CUT}): {len(ind)} / {len(df)} "
          f"({100*len(ind)/len(df):.0f}%)")

    # ---- couplings on in-domain set ---------------------------------------
    rows = []
    big = ind.bio_class.value_counts()
    for cls in big[big >= 100].index:
        m = ind.bio_class == cls
        sub, dd = ind[m], ind[m].drop_duplicates("murcko_scaffold")
        for t in TARGETS:
            a, b = ind.loc[m, t].values, ind.loc[~m, t].values
            if len(np.unique(a)) < 2:
                continue
            U, p = mannwhitneyu(a, b, alternative="greater")
            mf, md = float(np.median(a)), float(dd[t].median())
            rows.append({"bio_class": cls, "target": t, "n": int(m.sum()),
                         "n_scaf": dd.murcko_scaffold.nunique(),
                         "effect": round(U/(len(a)*len(b)), 3),
                         "median_full": round(mf, 3), "median_dedup": round(md, 3),
                         "robust": round(md/(mf+1e-9), 2), "p": p})
    res = pd.DataFrame(rows); res["q"] = bh(res.p.values)
    res.to_csv(os.path.join(HERE, "atlas_couplings_all.csv"), index=False)

    keep = res[(res.q < QCUT) & (res.effect > EFF) & (res.robust >= ROBUST)]
    # promiscuity: class significant-robust for >50% of targets = catch-all
    frac = keep.bio_class.value_counts() / len(TARGETS)
    promisc = set(frac[frac > 0.5].index)
    specific = keep[~keep.bio_class.isin(promisc)].sort_values("effect", ascending=False)
    specific.to_csv(os.path.join(HERE, "atlas_couplings_specific.csv"), index=False)
    print(f"\nrobust couplings: {len(keep)}  | promiscuous classes excluded: {sorted(promisc)}")
    print(f"SPECIFIC couplings (the honest atlas): {len(specific)}")
    print(specific.head(20)[["bio_class","target","n","effect","median_full","median_dedup","robust"]].to_string(index=False))

    # ---- AD-inside polypharmacology leads ---------------------------------
    ind["n_hits"] = (ind[TARGETS] > 0.5).sum(1)
    ind["mean_prob"] = ind[TARGETS].mean(1).round(3)
    leads = ind.sort_values(["n_hits", "mean_prob"], ascending=False)[
        ["bio_class", "ad", "n_hits", "mean_prob"]]
    leads.to_csv(os.path.join(HERE, "atlas_leads_AD.csv"))
    print(f"\nAD-inside multi-target leads (>=3 hits): {int((ind.n_hits>=3).sum())}")

    assert not specific.empty, "no specific couplings survived -- claim collapses"
    print("check ok -> atlas_couplings_specific.csv, atlas_leads_AD.csv")


if __name__ == "__main__":
    run()
