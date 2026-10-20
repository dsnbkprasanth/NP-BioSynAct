#!/usr/bin/env python
"""
NP-BioSynAct -- retrospective recovery validation (no wet lab).
Do known NP drugs rank high on their real target in np_predictions?
Match on InChIKey connectivity block (first 14 chars) so stereo/salt differences
don't break the lookup.

Run: conda activate npbs && python 31_validate.py
"""
import os
import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

# (name, SMILES, expected target) -- textbook NP pharmacology
KNOWN = [
    ("capsaicin",    "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",                    "TRPV1"),
    ("morphine",     "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "OPRM1"),
    ("curcumin",     "O=C(/C=C/c1ccc(O)c(OC)c1)CC(=O)/C=C/c1ccc(O)c(OC)c1",   "PTGS2"),
    ("resveratrol",  "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",                         "SIRT1"),
    ("berberine",    "COc1ccc2cc3[n+](cc2c1OC)CCc1cc2c(cc1-3)OCO2",           "DPP4"),
    ("quercetin",    "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",            "PTGS2"),
    ("cannabidiol",  "CCCCCc1cc(O)c(C2CC(C)=CCC2C(=C)C)c(O)c1",               "CNR2"),
    ("caffeine",     "Cn1cnc2c1c(=O)n(C)c(=O)n2C",                             "ADORA1"),
    ("serotonin",    "NCCc1c[nH]c2ccc(O)cc12",                                 "HTR1A"),
    ("cinnamaldehyde","O=C/C=C/c1ccccc1",                                      "TRPA1"),
    ("THC",          "CCCCCc1cc(O)c2c(c1)OC(C)(C)[C@@H]1CCC(C)=C[C@H]21",      "CNR2"),
]


def block(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToInchiKey(m)[:14] if m else None


def run():
    pred = pd.read_csv(os.path.join(HERE, "np_predictions.csv"))
    pred["blk"] = pred.inchikey.str[:14]
    # training-set membership: an anchor seen in training is memorised, not generalised
    train_blk = set(pd.read_csv(os.path.join(HERE, "chembl_activities.csv"),
                    usecols=["smiles"]).smiles.map(block).dropna())
    hits = gen = 0
    print(f"{'compound':13s} {'target':7s} {'prob':>6s} {'pctile':>7s}  in_train")
    for name, smi, tgt in KNOWN:
        b = block(smi)
        rows = pred[pred.blk == b]
        if rows.empty or tgt not in pred.columns:
            print(f"{name:13s} {tgt:7s}   not in library"); continue
        p = rows[tgt].max()
        pct = (pred[tgt] < p).mean() * 100          # percentile among all NPs
        top = pct >= 90
        seen = b in train_blk
        hits += top
        gen += top and not seen                      # genuine = held-out recovery
        print(f"{name:13s} {tgt:7s} {p:6.3f} {pct:6.1f}%{'  <-- top10%' if top else ''}"
              f"  {'yes' if seen else 'no'}")
    print(f"\n{hits}/{len(KNOWN)} recovered top-10%; {gen} GENUINE (held-out, not in training)")
    print("note: in_train=yes recoveries are memorisation, report only held-out ones")


if __name__ == "__main__":
    run()
