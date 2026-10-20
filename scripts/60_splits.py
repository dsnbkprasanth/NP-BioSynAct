#!/usr/bin/env python
"""Compare random and scaffold-separated splits to quantify analogue-leakage optimism.
Scaffold is the honest generalisation number. (Temporal split omitted: ChEMBL dates
not pulled; note as limitation.)

Run: python 60_splits.py    Out: split_comparison.csv
"""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from npbs_common import load_matrix, train, predict, scaffold_split, HERE
import os

df, genes, X, Y, M = load_matrix()

def macro(te):
    net = train(X, Y, M, ~te, genes, epochs=30); P = predict(net, X[te])
    Yte, Mte = Y[te], M[te]; a = []
    for j in range(len(genes)):
        msk = Mte[:, j]; yt = Yte[msk, j]
        if yt.size > 1 and yt.min() != yt.max():
            a.append(roc_auc_score(yt, P[msk, j]))
    return float(np.mean(a))

rng = np.random.default_rng(0); n = len(df)
teR = np.zeros(n, bool); teR[rng.choice(n, int(.2*n), replace=False)] = True
teS = scaffold_split(df, 0.2, 0)
r, s = macro(teR), macro(teS)
pd.DataFrame({"split": ["random", "scaffold"], "macro_auc": [r, s]}).to_csv(
    os.path.join(HERE, "split_comparison.csv"), index=False)
print(f"random-split   macro AUC: {r:.3f}")
print(f"scaffold-split macro AUC: {s:.3f}")
print(f"optimism gap (random-scaffold): {r-s:.3f}  -> report scaffold as honest performance")
assert s > 0.6, "scaffold AUC too low"
print("check ok -> split_comparison.csv")
