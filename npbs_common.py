"""Shared helpers for NP-BioSynAct rigor scripts (60-65). Import, don't duplicate."""
import os
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np, pandas as pd, torch, torch.nn as nn
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DEV = "cuda" if torch.cuda.is_available() else "cpu"
_G = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
RELIABLE = ["P2RX3","TRPA1","HTR1A","TNF","ADRB2","NTRK1","OPRM1","PPARG","SCN9A",
            "PTPN1","AKR1B1","CNR2","SIRT1","GSK3B","NOS2","ADORA1","SCN10A",
            "MAPK14","DPP4","NFKB1"]


def fp(s):
    m = Chem.MolFromSmiles(s); return np.array(_G.GetFingerprint(m), np.float32) if m else None
def fpmol(s):
    m = Chem.MolFromSmiles(s); return _G.GetFingerprint(m) if m else None
def scaf(s):
    try: return MurckoScaffoldSmiles(s, includeChirality=False)
    except Exception: return s


class Net(nn.Module):
    def __init__(self, o, p=0.3):
        super().__init__()
        self.b = nn.Sequential(nn.Linear(2048,512), nn.ReLU(), nn.Dropout(p),
                               nn.Linear(512,256), nn.ReLU(), nn.Dropout(p))
        self.h = nn.Linear(256, o)
    def forward(self, x): return self.h(self.b(x))


def load_matrix():
    df = pd.read_csv(os.path.join(HERE, "training_matrix.csv"))
    genes = [c for c in df.columns if c not in ("inchikey", "smiles", "is_np")]
    df["fpv"] = df.smiles.map(fp); df = df.dropna(subset=["fpv"]).reset_index(drop=True)
    X = np.stack(df.fpv); Y = df[genes].values.astype(np.float32)
    M = ~np.isnan(Y); Y = np.nan_to_num(Y)
    return df, genes, X, Y, M


def train(X, Y, M, tr, genes, epochs=30, p=0.3):
    Xtr = torch.tensor(X[tr], device=DEV); Ytr = torch.tensor(Y[tr], device=DEV)
    Mtr = torch.tensor(M[tr].astype(np.float32), device=DEV)
    pos = (Y[tr]*M[tr]).sum(0); neg = ((1-Y[tr])*M[tr]).sum(0)
    pw = torch.tensor(neg/np.clip(pos, 1, None), dtype=torch.float32, device=DEV)
    net = Net(len(genes), p).to(DEV); opt = torch.optim.Adam(net.parameters(), 1e-3, weight_decay=1e-5)
    lf = nn.BCEWithLogitsLoss(pos_weight=pw, reduction="none"); net.train()
    for _ in range(epochs):
        for i in range(0, len(Xtr), 512):
            opt.zero_grad(); m = Mtr[i:i+512]
            (lf(net(Xtr[i:i+512]), Ytr[i:i+512])*m).sum().div(m.sum().clamp(min=1)).backward(); opt.step()
    return net


def predict(net, Xsub, mc=False):
    (net.train() if mc else net.eval())
    with torch.no_grad():
        return torch.sigmoid(net(torch.tensor(Xsub, device=DEV))).cpu().numpy()


def scaffold_split(df, frac=0.2, seed=0):
    sc = df.smiles.map(scaf).values
    uq = pd.Series(np.unique(sc)).sample(frac=1, random_state=seed)
    return np.isin(sc, uq.iloc[:int(frac*len(uq))].values)
