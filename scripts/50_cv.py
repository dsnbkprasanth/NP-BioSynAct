#!/usr/bin/env python
"""
NP-BioSynAct -- repeated scaffold CV (rigor). 5 independent scaffold splits ->
mean +/- sd per-target and macro AUC, vs 1-NN Tanimoto baseline. Kills the
single-split / instability critique.

Run: conda activate npbs && python 50_cv.py    (~5 min GPU)
Out: cv_results.csv
"""
import os
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np, pandas as pd, torch, torch.nn as nn
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
DEV = "cuda" if torch.cuda.is_available() else "cpu"
G = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
REPEATS = 5


def fp(s):
    m = Chem.MolFromSmiles(s); return np.array(G.GetFingerprint(m), np.float32) if m else None
def scaf(s):
    try: return MurckoScaffoldSmiles(s, includeChirality=False)
    except Exception: return s


class Net(nn.Module):
    def __init__(s, o):
        super().__init__()
        s.b = nn.Sequential(nn.Linear(2048,512), nn.ReLU(), nn.Dropout(.3),
                            nn.Linear(512,256), nn.ReLU(), nn.Dropout(.3)); s.h = nn.Linear(256,o)
    def forward(s, x): return s.h(s.b(x))


def main():
    df = pd.read_csv(os.path.join(HERE, "training_matrix.csv"))
    genes = [c for c in df.columns if c not in ("inchikey","smiles","is_np")]
    df["fpv"] = df.smiles.map(fp); df = df.dropna(subset=["fpv"]).reset_index(drop=True)
    X = np.stack(df.fpv); Y = df[genes].values.astype(np.float32); M = ~np.isnan(Y); Y = np.nan_to_num(Y)
    sc = df.smiles.map(scaf).values
    rec = {g: [] for g in genes}; macro, bmacro = [], []
    for r in range(REPEATS):
        uq = pd.Series(np.unique(sc)).sample(frac=1, random_state=r)
        te = np.isin(sc, uq.iloc[:int(.2*len(uq))].values); tr = ~te  # ndarray, not set
        Xtr = torch.tensor(X[tr], device=DEV); Ytr = torch.tensor(Y[tr], device=DEV)
        Mtr = torch.tensor(M[tr].astype(np.float32), device=DEV); Xte = torch.tensor(X[te], device=DEV)
        pos = (Y[tr]*M[tr]).sum(0); neg = ((1-Y[tr])*M[tr]).sum(0)
        pw = torch.tensor(neg/np.clip(pos,1,None), dtype=torch.float32, device=DEV)
        net = Net(len(genes)).to(DEV); opt = torch.optim.Adam(net.parameters(), 1e-3, weight_decay=1e-5)
        lf = nn.BCEWithLogitsLoss(pos_weight=pw, reduction="none"); net.train()
        for _ in range(25):
            for i in range(0, len(Xtr), 512):
                opt.zero_grad(); m = Mtr[i:i+512]
                (lf(net(Xtr[i:i+512]), Ytr[i:i+512])*m).sum().div(m.sum().clamp(min=1)).backward(); opt.step()
        net.eval()
        with torch.no_grad():
            P = torch.sigmoid(net(Xte)).cpu().numpy()
        Yte, Mte = Y[te], M[te]; aa, bb = [], []
        te_sum = Xte.sum(1, keepdim=True)
        for j, g in enumerate(genes):
            msk = Mte[:, j]; yt = Yte[msk, j]
            if yt.size < 2 or yt.min() == yt.max(): continue
            rec[g].append(roc_auc_score(yt, P[msk, j])); aa.append(rec[g][-1])
            # baseline: max Tanimoto to this target's train ACTIVES only (small, no OOM)
            ai = np.where((Y[tr][:, j] == 1) & M[tr][:, j])[0]
            if len(ai):
                with torch.no_grad():
                    Xa = Xtr[torch.tensor(ai, device=DEV)]; c = Xte @ Xa.t()
                    sim = (c / (te_sum + Xa.sum(1)[None, :] - c + 1e-6)).max(1).values.cpu().numpy()
                bb.append(roc_auc_score(yt, sim[msk]))
        macro.append(np.mean(aa)); bmacro.append(np.mean(bb))
        del Xtr, Ytr, Mtr, Xte
        if DEV == "cuda": torch.cuda.empty_cache()
        print(f"  repeat {r+1}: model {macro[-1]:.3f}  baseline {bmacro[-1]:.3f}")
    out = pd.DataFrame({"target": list(rec), "auc_mean": [np.mean(v) for v in rec.values()],
                        "auc_sd": [np.std(v) for v in rec.values()]}).sort_values("auc_mean", ascending=False)
    out.to_csv(os.path.join(HERE, "cv_results.csv"), index=False)
    print(out.to_string(index=False))
    print(f"\nMODEL macro AUC {np.mean(macro):.3f} +/- {np.std(macro):.3f}")
    print(f"BASELINE macro  {np.mean(bmacro):.3f} +/- {np.std(bmacro):.3f}")
    print(f"model beats baseline in {sum(m>b for m,b in zip(macro,bmacro))}/{REPEATS} repeats")
    assert np.mean(macro) > 0.7, "CV AUC too low"
    print("check ok -> cv_results.csv")


if __name__ == "__main__":
    main()
