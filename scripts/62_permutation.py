#!/usr/bin/env python
"""Permutation null controls.
 (A) ML signal: shuffle target labels, retrain, 50x -> null AUC vs real.
 (B) Coupling: shuffle biosynthetic-class labels, 1000x -> null enrichment vs observed,
     proving the biosynthesis<->target association is not random class structure.

Run: python 62_permutation.py    Out: perm_model.csv, perm_coupling.csv
"""
import numpy as np, pandas as pd, os
from sklearn.metrics import roc_auc_score
from npbs_common import load_matrix, train, predict, scaffold_split, HERE, RELIABLE

# ---- (A) model label-permutation ------------------------------------------
df, genes, X, Y, M = load_matrix()
te = scaffold_split(df, 0.2, 0); tr = ~te
def macro(Ymat):
    net = train(X, Ymat, M, tr, genes, epochs=25); P = predict(net, X[te])
    a = []
    for j in range(len(genes)):
        msk = M[te][:, j]; yt = Ymat[te][msk, j]
        if yt.size > 1 and yt.min() != yt.max():
            a.append(roc_auc_score(yt, P[msk, j]))
    return float(np.mean(a))
real = macro(Y)
rng = np.random.default_rng(0); null = []
for k in range(50):
    Yp = Y.copy(); idx = rng.permutation(len(Yp)); Yp = Yp[idx]  # break compound-label link
    null.append(macro(Yp))
null = np.array(null); p = (np.sum(null >= real) + 1) / (len(null) + 1)
pd.DataFrame({"real_auc": [real], "null_mean": [null.mean()], "null_max": [null.max()],
             "p_value": [p]}).to_csv(os.path.join(HERE, "perm_model.csv"), index=False)
print(f"(A) model: real {real:.3f} vs null {null.mean():.3f}+/-{null.std():.3f} (max {null.max():.3f}), p={p:.3f}")

# ---- (B) coupling class-permutation ---------------------------------------
pred = pd.read_csv(os.path.join(HERE, "np_predictions.csv")).set_index("inchikey")
bio = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                  usecols=["inchikey", "bio_class"], low_memory=False).dropna().set_index("inchikey")
d = pred[RELIABLE].join(bio, how="inner")
big = d.bio_class.value_counts(); classes = big[big >= 100].index
# observed statistic: mean abs deviation of class-mean from global mean, summed
glob = d[RELIABLE].mean()
def stat(labels):
    g = d[RELIABLE].groupby(labels).mean().loc[classes]
    return float((g - glob).abs().values.sum())
obs = stat(d.bio_class.values)
rng = np.random.default_rng(1); nc = []
lab = d.bio_class.values.copy()
for k in range(1000):
    nc.append(stat(rng.permutation(lab)))
nc = np.array(nc); pC = (np.sum(nc >= obs) + 1) / (len(nc) + 1)
pd.DataFrame({"observed": [obs], "null_mean": [nc.mean()], "null_p995": [np.percentile(nc, 99.5)],
             "p_value": [pC]}).to_csv(os.path.join(HERE, "perm_coupling.csv"), index=False)
print(f"(B) coupling: observed {obs:.2f} vs null {nc.mean():.2f} (99.5th {np.percentile(nc,99.5):.2f}), p={pC:.3f}")
assert p < 0.05 and pC < 0.05, "a permutation test is not significant -- report honestly"
print("check ok -> perm_model.csv, perm_coupling.csv")
