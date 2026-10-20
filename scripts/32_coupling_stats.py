#!/usr/bin/env python
"""
NP-BioSynAct -- statistical hardening of the biosynthesis<->activity coupling.

Statistical safeguards:
  1. significance: Mann-Whitney U per (biosynthetic class x target) + BH-FDR
  2. scaffold confound: recompute enrichment after Murcko-scaffold dedup within
     each class -> shows coupling is not one analog series
  3. applicability domain: max Tanimoto of top-ranked NPs to ChEMBL space

Run: conda activate npbs && python 32_coupling_stats.py
Out: coupling_stats.csv
"""
import os
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

AUC = {"P2RX3":.993,"TRPA1":.990,"HTR1A":.986,"TNF":.981,"NTRK1":.969,"CNR2":.958,
       "OPRM1":.957,"SCN9A":.954,"PPARG":.952,"ADRB2":.950,"SIRT1":.949,"GSK3B":.944,
       "PTPN1":.940,"AKR1B1":.933,"NOS2":.916,"ADORA1":.887,"SCN10A":.886,"DPP4":.859,
       "NLRP3":.806,"MAPK14":.785,"GRM5":.781,"NFKB1":.742,"PTGS2":.740,"TRPV1":.624,
       "IL6":1.,"RAGE":1.,"TLR4":1.,"IL1B":.333}
RELIABLE = [t for t, a in AUC.items() if 0.70 <= a <= 0.99]


def bh_fdr(p):
    p = np.asarray(p); n = len(p); order = p.argsort()
    q = np.empty(n); q[order] = (p[order] * n / (np.arange(n) + 1))
    # enforce monotonicity
    q[order] = np.minimum.accumulate(q[order][::-1])[::-1]
    return np.clip(q, 0, 1)


def fp(smi):
    m = Chem.MolFromSmiles(smi)
    return FPGEN.GetFingerprint(m) if m else None


def run():
    pred = pd.read_csv(os.path.join(HERE, "np_predictions.csv")).set_index("inchikey")
    bio = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                      usecols=["inchikey", "smiles", "bio_class", "murcko_scaffold"]
                      ).dropna(subset=["bio_class"]).set_index("inchikey")
    df = pred[RELIABLE].join(bio, how="inner")
    print(f"labelled NPs: {len(df)}  classes: {df.bio_class.nunique()}")

    # ---- 1 & 2: significance + scaffold-dedup enrichment -------------------
    rows = []
    big = df.bio_class.value_counts()
    for cls in big[big >= 100].index:
        inm = df.bio_class == cls
        sub = df[inm]
        dedup = sub.drop_duplicates("murcko_scaffold")     # one per scaffold
        for t in RELIABLE:
            a, b = df.loc[inm, t].values, df.loc[~inm, t].values
            if len(np.unique(a)) < 2:
                continue
            U, p = mannwhitneyu(a, b, alternative="greater")
            cl = U / (len(a) * len(b))                      # P(in-class > out) effect
            rows.append({"bio_class": cls, "target": t, "n": int(inm.sum()),
                         "n_scaffolds": dedup.murcko_scaffold.nunique(),
                         "effect_cl": round(cl, 3),
                         "median_full": round(np.median(a), 3),
                         "median_dedup": round(float(dedup[t].median()), 3),
                         "p": p})
    res = pd.DataFrame(rows)
    res["q"] = bh_fdr(res.p.values)
    sig = res[(res.q < 0.05) & (res.effect_cl > 0.6)].sort_values("effect_cl", ascending=False)
    res.sort_values("q").to_csv(os.path.join(HERE, "coupling_stats.csv"), index=False)
    print(f"\nsignificant enriched couplings (q<0.05, effect>0.6): {len(sig)}")
    print(sig.head(15).to_string(index=False))
    # confound check: does dedup kill it? correlate full vs dedup medians on sig set
    if len(sig) > 3:
        r = np.corrcoef(sig.median_full, sig.median_dedup)[0, 1]
        print(f"\nfull-vs-scaffold-dedup median corr on sig couplings: {r:.2f} "
              f"({'robust' if r > 0.7 else 'CONFOUNDED -- driven by analog series'})")

    # ---- 3: applicability domain of the top polypharmacology shortlist -----
    poly = pd.read_csv(os.path.join(HERE, "np_polypharm_ranked.csv")).head(200)
    poly = poly.merge(df[["smiles"]], left_on="inchikey", right_index=True, how="left")
    ref = pd.read_csv(os.path.join(HERE, "chembl_activities.csv"),
                      usecols=["smiles"]).smiles.dropna().sample(3000, random_state=0)
    ref_fp = [x for x in (fp(s) for s in ref) if x]
    ad = []
    for s in poly.smiles.dropna():
        q = fp(s)
        ad.append(max(DataStructs.BulkTanimotoSimilarity(q, ref_fp)) if q else 0)
    ad = np.array(ad)
    print(f"\ntop-200 leads applicability domain: "
          f"{(ad >= 0.3).mean()*100:.0f}% inside (maxTanimoto>=0.3), median {np.median(ad):.2f}")

    assert not res.empty, "no couplings computed"
    print("check ok -> coupling_stats.csv")


if __name__ == "__main__":
    run()
