#!/usr/bin/env python
"""
NP-BioSynAct -- publication figures (300 dpi).
Fig1 biosynthesis x target enrichment heatmap | Fig2 per-target AUC | Fig3 recovery.

Run: conda activate npbs && python 40_figures.py
Out: fig1_coupling.png, fig2_auc.png, fig3_recovery.png
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300,
                     "savefig.bbox": "tight", "axes.titlesize": 10})

AUC = {"P2RX3":.993,"TRPA1":.990,"HTR1A":.986,"TNF":.981,"NTRK1":.969,"CNR2":.958,
       "OPRM1":.957,"SCN9A":.954,"PPARG":.952,"ADRB2":.950,"SIRT1":.949,"GSK3B":.944,
       "PTPN1":.940,"AKR1B1":.933,"NOS2":.916,"ADORA1":.887,"SCN10A":.886,"DPP4":.859,
       "NLRP3":.806,"MAPK14":.785,"GRM5":.781,"NFKB1":.742,"PTGS2":.740,"TRPV1":.624,
       "IL6":1.0,"RAGE":1.0,"TLR4":1.0,"IL1B":.333}
RECOVERY = [("capsaicin","TRPV1",96.0),("morphine","OPRM1",97.5),("curcumin","PTGS2",0.0),
            ("resveratrol","SIRT1",50.8),("berberine","DPP4",30.9),("quercetin","PTGS2",47.6),
            ("cannabidiol","CNR2",90.6),("caffeine","ADORA1",7.0),("serotonin","HTR1A",96.3),
            ("cinnamaldehyde","TRPA1",27.4),("THC","CNR2",90.6)]


def fig1():
    f = os.path.join(HERE, "coupling_pathway_log2enrich.csv")
    if not os.path.exists(f):
        print("  fig1 skipped: run 30_coupling.py first"); return
    m = pd.read_csv(f, index_col=0)
    v = np.nanmax(np.abs(m.values))
    fig, ax = plt.subplots(figsize=(11, 0.5 * len(m) + 1.5))
    im = ax.imshow(m.values, aspect="auto", cmap="RdBu_r", vmin=-v, vmax=v)
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=90)
    ax.set_yticks(range(len(m))); ax.set_yticklabels(m.index)
    ax.set_title("Biosynthetic pathway x target enrichment (log2 vs global mean)")
    fig.colorbar(im, label="log2 enrichment")
    fig.savefig(os.path.join(HERE, "fig1_coupling.png")); plt.close(fig)
    print("  fig1_coupling.png")


def fig2():
    s = pd.Series(AUC).sort_values()
    reliable = (s >= 0.70) & (s <= 0.99)
    colors = ["#b0b0b0" if not r else "#2c7fb8" for r in reliable]
    fig, ax = plt.subplots(figsize=(6, 7))
    ax.barh(range(len(s)), s.values, color=colors)
    ax.axvline(0.7, ls="--", c="k", lw=.8)
    ax.set_yticks(range(len(s))); ax.set_yticklabels(s.index)
    ax.set_xlabel("scaffold-split ROC-AUC"); ax.set_xlim(0, 1.02)
    ax.set_title("Per-target performance (grey = excluded: broken or small-n)")
    fig.savefig(os.path.join(HERE, "fig2_auc.png")); plt.close(fig)
    print("  fig2_auc.png")


def fig3():
    r = pd.DataFrame(RECOVERY, columns=["cpd", "tgt", "pct"]).sort_values("pct")
    lab = r.cpd + " -> " + r.tgt
    col = ["#d95f02" if p >= 90 else "#7570b3" for p in r.pct]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.hlines(range(len(r)), 0, r.pct, color=col, lw=2)
    ax.plot(r.pct, range(len(r)), "o", color="k", ms=4)
    ax.axvline(90, ls="--", c="k", lw=.8)
    ax.set_yticks(range(len(r))); ax.set_yticklabels(lab)
    ax.set_xlabel("percentile among all NPs (target)"); ax.set_xlim(0, 100)
    ax.set_title("Retrospective recovery of known NP drugs (orange = top-10%)")
    fig.savefig(os.path.join(HERE, "fig3_recovery.png")); plt.close(fig)
    print("  fig3_recovery.png")


if __name__ == "__main__":
    fig1(); fig2(); fig3()
    print("done")
