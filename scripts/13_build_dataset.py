#!/usr/bin/env python
"""
NP-BioSynAct -- assemble the ML training matrix with REAL negatives.

chembl_activities.csv (actives + inactives, all carry pchembl) -> dedup by
InChIKey -> compound x 29-target label matrix:
    active   pchembl >= 6.5  -> 1
    inactive pchembl <  5.0  -> 0
    ambiguous 5.0-6.5 / not measured -> blank (unknown, masked in training)

Out: training_matrix.csv  ->  inchikey, smiles, is_np, <one col per gene: 1/0/blank>

Run: conda activate npbs && python 13_build_dataset.py
"""
import os
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ACT = os.path.join(HERE, "chembl_activities.csv")
NP = os.path.join(HERE, "np_compounds.csv")
OUT = os.path.join(HERE, "training_matrix.csv")
ACT_HI, INACT_LO = 6.5, 5.0


def inchikey(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToInchiKey(m) if m else None


def run():
    act = pd.read_csv(ACT)
    act = act.dropna(subset=["pchembl"])
    print(f"activities: {len(act)} rows, {act.gene.nunique()} targets, "
          f"pchembl {act.pchembl.min():.1f}-{act.pchembl.max():.1f}")

    smi2ik = {s: inchikey(s) for s in act.smiles.unique()}
    act["inchikey"] = act.smiles.map(smi2ik)
    act = act.dropna(subset=["inchikey"])

    # max pchembl per (compound, target)
    piv = act.pivot_table(index="inchikey", columns="gene",
                          values="pchembl", aggfunc="max")
    # threshold -> 1 / 0 / NaN(unknown)
    lab = pd.DataFrame(np.where(piv >= ACT_HI, 1.0,
                        np.where(piv < INACT_LO, 0.0, np.nan)),
                       index=piv.index, columns=piv.columns)
    lab = lab.dropna(how="all")                       # drop all-unknown compounds

    smi = act.drop_duplicates("inchikey").set_index("inchikey").smiles
    lab.insert(0, "smiles", smi.reindex(lab.index))
    np_ik = set(pd.read_csv(NP, usecols=["inchikey"]).inchikey)
    lab.insert(1, "is_np", lab.index.isin(np_ik).astype("int8"))

    lab.to_csv(OUT)
    genes = [c for c in lab.columns if c not in ("smiles", "is_np")]
    pos = (lab[genes] == 1).sum()
    neg = (lab[genes] == 0).sum()
    print(f"compounds: {len(lab)}  |  NP overlap: {int(lab.is_np.sum())}")
    print("per target  active / inactive:")
    for g in genes:
        print(f"  {g:10s} {int(pos[g]):6d} / {int(neg[g]):6d}")

    # self-check: every target has BOTH classes (else AUC undefined / leakage risk)
    assert (pos > 0).all() and (neg > 0).all(), "a target lacks one class"
    print("check ok ->", OUT)


if __name__ == "__main__":
    run()
