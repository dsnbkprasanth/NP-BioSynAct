#!/usr/bin/env python
"""
NP-BioSynAct -- the novel analysis: biosynthesis <-> polypharmacology coupling.

Joins NP-library predictions (20) with NPClassifier biosynthetic labels (11),
then asks: which biosynthetic pathways/classes are ENRICHED for activity against
pain / chronic-disease targets? Also ranks the most polypharmacological NPs.

Reliable targets only (0.70 <= scaffold-AUC <= 0.99): drops broken IL1B (0.33)
and small-n perfect IL6/RAGE/TLR4 (1.000).

Run: conda activate npbs && python 30_coupling.py
Out: coupling_pathway.csv, coupling_class.csv, np_polypharm_ranked.csv, fig_coupling.png
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# scaffold-split AUC from 20_train_multitask.py -- keep only trustworthy heads
AUC = {"P2RX3":.993,"TRPA1":.990,"HTR1A":.986,"TNF":.981,"NTRK1":.969,"CNR2":.958,
       "OPRM1":.957,"SCN9A":.954,"PPARG":.952,"ADRB2":.950,"SIRT1":.949,"GSK3B":.944,
       "PTPN1":.940,"AKR1B1":.933,"NOS2":.916,"ADORA1":.887,"SCN10A":.886,"DPP4":.859,
       "NLRP3":.806,"MAPK14":.785,"GRM5":.781,"NFKB1":.742,"PTGS2":.740,"TRPV1":.624,
       "IL6":1.0,"RAGE":1.0,"TLR4":1.0,"IL1B":.333}
RELIABLE = [t for t, a in AUC.items() if 0.70 <= a <= 0.99]


def run():
    pred = pd.read_csv(os.path.join(HERE, "np_predictions.csv")).set_index("inchikey")
    bio = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                      usecols=["inchikey", "bio_pathway", "bio_superclass", "bio_class"])
    bio = bio.dropna(subset=["bio_pathway"]).set_index("inchikey")
    df = pred[RELIABLE].join(bio, how="inner")
    print(f"labelled NPs with predictions: {len(df)}  |  targets: {len(RELIABLE)}")

    glob = df[RELIABLE].mean()                         # global mean prob per target

    def couple(level, out, min_n=50):
        grp = df.groupby(level)
        big = grp.size()[grp.size() >= min_n].index    # skip tiny classes
        mean = grp[RELIABLE].mean().loc[big]
        enr = np.log2(mean.div(glob) + 1e-9)           # log2 enrichment vs global
        mean.to_csv(os.path.join(HERE, f"{out}_mean.csv"))
        enr.round(3).to_csv(os.path.join(HERE, f"{out}_log2enrich.csv"))
        return mean, enr

    pw_mean, pw_enr = couple("bio_pathway", "coupling_pathway", min_n=50)
    cl_mean, cl_enr = couple("bio_class", "coupling_class", min_n=100)

    # top biosynthesis<->target couplings (class level)
    top = (cl_enr.stack().rename("log2enr").reset_index()
             .sort_values("log2enr", ascending=False).head(20))
    top.columns = ["bio_class", "target", "log2enr"]
    print("\ntop biosynthetic-class / target couplings (log2 enrichment):")
    print(top.to_string(index=False))

    # polypharmacology: breadth per NP over reliable targets
    poly = pd.DataFrame({
        "bio_pathway": df.bio_pathway, "bio_class": df.bio_class,
        "mean_prob": df[RELIABLE].mean(1).round(3),
        "n_hits_p50": (df[RELIABLE] > 0.5).sum(1),
    }).sort_values(["n_hits_p50", "mean_prob"], ascending=False)
    poly.to_csv(os.path.join(HERE, "np_polypharm_ranked.csv"))
    print(f"\ntop polypharmacological NPs (multi-target): "
          f"{int((poly.n_hits_p50 >= 3).sum())} hit >=3 targets")

    # headline figure: pathway x target mean-probability heatmap
    fig, ax = plt.subplots(figsize=(11, max(3, 0.5 * len(pw_mean))))
    im = ax.imshow(pw_mean.values, aspect="auto", cmap="magma", vmin=0, vmax=pw_mean.values.max())
    ax.set_xticks(range(len(RELIABLE))); ax.set_xticklabels(RELIABLE, rotation=90, fontsize=8)
    ax.set_yticks(range(len(pw_mean))); ax.set_yticklabels(pw_mean.index, fontsize=9)
    ax.set_title("Biosynthetic pathway x target: mean predicted activity")
    fig.colorbar(im, label="mean P(active)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_coupling.png"), dpi=200)

    assert not pw_mean.empty and not pw_enr.isna().all().all(), "empty coupling"
    print("check ok -> coupling_*.csv, np_polypharm_ranked.csv, fig_coupling.png")


if __name__ == "__main__":
    run()
