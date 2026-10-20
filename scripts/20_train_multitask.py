#!/usr/bin/env python
"""
NP-BioSynAct -- multitask polypharmacology model (fingerprint baseline / engine).

Morgan FP -> shared MLP -> 29 sigmoid heads (one per target).
Scaffold split (no leakage). Reports per-target ROC-AUC on held-out scaffolds,
then predicts the whole NP library.

Run: conda activate npbs && python 20_train_multitask.py
Out: np_predictions.csv  (inchikey, is_np, <29 target probabilities>)
"""
import os
import numpy as np, pandas as pd
import torch, torch.nn as nn
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DEV = "cuda" if torch.cuda.is_available() else "cpu"
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fp(smi):
    m = Chem.MolFromSmiles(smi)
    return np.array(FPGEN.GetFingerprint(m), dtype=np.float32) if m else None


def scaffold(smi):
    try:
        return MurckoScaffoldSmiles(smi, includeChirality=False)
    except Exception:
        return smi


class Net(nn.Module):
    def __init__(self, n_out):
        super().__init__()
        self.body = nn.Sequential(nn.Linear(2048, 512), nn.ReLU(), nn.Dropout(0.3),
                                  nn.Linear(512, 256), nn.ReLU(), nn.Dropout(0.3))
        self.head = nn.Linear(256, n_out)

    def forward(self, x):
        return self.head(self.body(x))


def main():
    print("device:", DEV)
    df = pd.read_csv(os.path.join(HERE, "training_matrix.csv"))
    genes = [c for c in df.columns if c not in ("inchikey", "smiles", "is_np")]

    # featurize
    df["fpv"] = df.smiles.map(fp)
    df = df.dropna(subset=["fpv"]).reset_index(drop=True)
    X = np.stack(df.fpv.values)
    Y = df[genes].values.astype(np.float32)          # 1 / 0 / NaN(unknown)
    M = ~np.isnan(Y)                                  # mask: measured entries only
    Y = np.nan_to_num(Y)                             # NaN->0 (masked out of loss)

    # scaffold split 80/20: a whole scaffold lands in one side (no leakage)
    sc = df.smiles.map(scaffold)
    uniq = sc.drop_duplicates().sample(frac=1, random_state=0)
    test_sc = set(uniq.iloc[:int(0.2 * len(uniq))])
    te = sc.isin(test_sc).values
    tr = ~te
    print(f"compounds {len(df)}  train {tr.sum()}  test {te.sum()}")

    Xtr = torch.tensor(X[tr], device=DEV); Ytr = torch.tensor(Y[tr], device=DEV)
    Mtr = torch.tensor(M[tr].astype(np.float32), device=DEV)
    Xte = torch.tensor(X[te], device=DEV)

    # per-target pos_weight from MEASURED train entries only
    pos = (Y[tr] * M[tr]).sum(0); neg = ((1 - Y[tr]) * M[tr]).sum(0)
    pw = torch.tensor(neg / np.clip(pos, 1, None), dtype=torch.float32, device=DEV)

    net = Net(len(genes)).to(DEV)
    opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-5)
    lossf = nn.BCEWithLogitsLoss(pos_weight=pw, reduction="none")
    net.train()
    for ep in range(40):
        for i in range(0, len(Xtr), 512):
            opt.zero_grad()
            m = Mtr[i:i+512]
            per = lossf(net(Xtr[i:i+512]), Ytr[i:i+512]) * m   # ignore unknowns
            loss = per.sum() / m.sum().clamp(min=1)
            loss.backward(); opt.step()
        if ep % 10 == 9:
            print(f"  epoch {ep+1}  loss {loss.item():.3f}")

    # eval: per-target AUC on held-out scaffolds
    net.eval()
    with torch.no_grad():
        P = torch.sigmoid(net(Xte)).cpu().numpy()
    Yte = Y[te]; Mte = M[te]
    aucs = {}
    for j, g in enumerate(genes):
        msk = Mte[:, j]                                  # measured test entries only
        yt = Yte[msk, j]
        if yt.min() != yt.max():                         # both classes present
            aucs[g] = roc_auc_score(yt, P[msk, j])
    macro = np.mean(list(aucs.values()))
    print(f"\nscaffold-split macro ROC-AUC: {macro:.3f}  ({len(aucs)} targets)")
    for g, a in sorted(aucs.items(), key=lambda x: -x[1]):
        print(f"  {g:10s} {a:.3f}")
    assert macro > 0.6, "macro AUC too low -- something is wrong"

    # ---- similarity baseline: 1-NN max-Tanimoto to train actives ------------
    # Compare the multitask model against a nearest-neighbour similarity baseline.
    with torch.no_grad():
        c = Xte @ Xtr.t()
        tan = (c / (Xte.sum(1, keepdim=True) + Xtr.sum(1)[None, :] - c + 1e-6)).cpu().numpy()
    Ytr_np, Mtr_np = Y[tr], M[tr]
    base = {}
    for j, g in enumerate(genes):
        ai = np.where((Ytr_np[:, j] == 1) & Mtr_np[:, j])[0]
        msk = Mte[:, j]; yt = Yte[msk, j]
        if len(ai) and yt.min() != yt.max():
            base[g] = roc_auc_score(yt, tan[:, ai].max(1)[msk])
    print(f"\nsimilarity-baseline macro AUC: {np.mean(list(base.values())):.3f}  "
          f"(model {macro:.3f})")

    # ---- bootstrap 95% CI on model macro AUC -------------------------------
    rng = np.random.default_rng(0); ms = []
    for _ in range(200):
        aa = []
        for j, g in enumerate(genes):
            if g not in aucs:
                continue
            idx = np.where(Mte[:, j])[0]
            bs = rng.choice(idx, len(idx), replace=True)
            yt = Yte[bs, j]
            if yt.min() != yt.max():
                aa.append(roc_auc_score(yt, P[bs, j]))
        ms.append(np.mean(aa))
    lo, hi = np.percentile(ms, [2.5, 97.5])
    print(f"model macro AUC 95% CI: [{lo:.3f}, {hi:.3f}]")

    # predict the whole NP library (chunked)
    np_df = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                        usecols=["inchikey", "smiles"])
    out = os.path.join(HERE, "np_predictions.csv")
    first = True
    for k in range(0, len(np_df), 20000):
        chunk = np_df.iloc[k:k+20000].copy()
        chunk["fpv"] = chunk.smiles.map(fp)
        chunk = chunk.dropna(subset=["fpv"])
        xb = torch.tensor(np.stack(chunk.fpv.values), device=DEV)
        with torch.no_grad():
            pb = torch.sigmoid(net(xb)).cpu().numpy()
        res = pd.DataFrame(pb, columns=genes, index=chunk.inchikey).round(3)
        res.to_csv(out, mode="w" if first else "a", header=first)
        first = False
        print(f"  predicted {k+len(chunk)}/{len(np_df)}")
    print("done ->", out)


if __name__ == "__main__":
    main()
