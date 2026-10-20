#!/usr/bin/env python
"""Prediction calibration and uncertainty.
MC-dropout (30 stochastic passes) -> per-prediction mean & SD (epistemic uncertainty).
Calibration curve, Brier score, Expected Calibration Error on the scaffold-test set.

Run: python 63_calibration.py    Out: calibration.csv, fig_calibration.png
"""
import numpy as np, pandas as pd, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.metrics import brier_score_loss
from sklearn.calibration import calibration_curve
from npbs_common import load_matrix, train, predict, scaffold_split, HERE

df, genes, X, Y, M = load_matrix()
te = scaffold_split(df, 0.2, 0); tr = ~te
net = train(X, Y, M, tr, genes, epochs=30)

# MC-dropout: keep dropout ON, 30 passes
runs = np.stack([predict(net, X[te], mc=True) for _ in range(30)])   # (30, n_te, T)
mean, sd = runs.mean(0), runs.std(0)

Yte, Mte = Y[te], M[te]
yt, pt, ut = [], [], []
for j in range(len(genes)):
    m = Mte[:, j]
    yt.append(Yte[m, j]); pt.append(mean[m, j]); ut.append(sd[m, j])
yt = np.concatenate(yt); pt = np.concatenate(pt); ut = np.concatenate(ut)

brier = brier_score_loss(yt, pt)
frac_pos, mean_pred = calibration_curve(yt, pt, n_bins=10, strategy="quantile")
ece = np.mean(np.abs(frac_pos - mean_pred))
pd.DataFrame({"metric": ["Brier", "ECE", "mean_uncertainty"],
             "value": [brier, ece, ut.mean()]}).to_csv(os.path.join(HERE, "calibration.csv"), index=False)
print(f"Brier {brier:.3f}  ECE {ece:.3f}  mean MC-dropout SD {ut.mean():.3f}")

fig, ax = plt.subplots(1, 2, figsize=(9, 4))
ax[0].plot([0, 1], [0, 1], "k--", lw=.8); ax[0].plot(mean_pred, frac_pos, "o-")
ax[0].set_xlabel("mean predicted prob"); ax[0].set_ylabel("observed frequency")
ax[0].set_title(f"Calibration (ECE={ece:.3f})")
ax[1].hist(ut, bins=40); ax[1].set_xlabel("MC-dropout SD (uncertainty)"); ax[1].set_ylabel("count")
ax[1].set_title("Prediction uncertainty")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_calibration.png"), dpi=200)
print("check ok -> calibration.csv, fig_calibration.png")
