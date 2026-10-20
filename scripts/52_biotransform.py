#!/usr/bin/env python
"""
NP-BioSynAct -- gut/hepatic biotransformation (SI's gut-brain-pain theme).
Predict phase-1/2 metabolites of top AD-inside leads (SyGMa), then re-score
polypharmacology of each metabolite vs its parent using a 1-NN-to-ChEMBL-actives
similarity scorer (model~=baseline, established). Flags metabolites that GAIN
multi-target activity after metabolism -> novel, on-theme finding.

Needs: pip install sygma
Run: conda activate npbs && python 52_biotransform.py
Out: biotransform.csv
"""
import os, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")
import sygma
HERE = os.path.dirname(os.path.abspath(__file__))
G = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
TARGETS = ["P2RX3","TRPA1","HTR1A","TNF","ADRB2","NTRK1","OPRM1","PPARG","SCN9A",
           "PTPN1","AKR1B1","CNR2","SIRT1","GSK3B","NOS2","ADORA1","SCN10A","MAPK14",
           "DPP4","NFKB1"]
N_LEADS = 100
SCN = sygma.Scenario([[sygma.ruleset['phase1'], 1], [sygma.ruleset['phase2'], 1]])


def fpmol(m): return G.GetFingerprint(m) if m else None
def fps(s):
    m = Chem.MolFromSmiles(s); return fpmol(m)


def build_actives():
    act = pd.read_csv(os.path.join(HERE, "chembl_activities.csv"))
    act = act[act.pchembl >= 6.5]
    d = {}
    for t in TARGETS:
        smis = act.loc[act.gene == t, "smiles"].dropna().unique()
        d[t] = [x for x in (fps(s) for s in smis) if x]
    return d


def polypharm(smi, actives):
    q = fps(smi)
    if not q:
        return None, {}
    per = {t: (max(DataStructs.BulkTanimotoSimilarity(q, actives[t])) if actives[t] else 0)
           for t in TARGETS}
    return sum(per.values()), per


def run():
    actives = build_actives()
    leads = pd.read_csv(os.path.join(HERE, "atlas_leads_AD.csv")).head(N_LEADS)
    smi = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                      usecols=["inchikey", "smiles"], low_memory=False).set_index("inchikey").smiles
    rows = []
    for ik in leads.inchikey:
        ps = smi.get(ik)
        m = Chem.MolFromSmiles(ps) if isinstance(ps, str) else None
        if not m:
            continue
        p_sum, _ = polypharm(ps, actives)
        tree = SCN.run(m); tree.calc_scores()
        mets = tree.to_smiles()                       # [(smiles, prob), ...] incl. parent
        best_gain = 0.0; best = None
        for msmi, prob in mets:
            if msmi == ps:
                continue
            s, _ = polypharm(msmi, actives)
            if s is not None and s - p_sum > best_gain:
                best_gain, best = s - p_sum, (msmi, prob, s)
        rows.append({"inchikey": ik, "parent_polypharm": round(p_sum, 3),
                     "n_metabolites": len(mets) - 1,
                     "best_metab_polypharm": round(best[2], 3) if best else round(p_sum, 3),
                     "gain": round(best_gain, 3),
                     "best_metabolite": best[0] if best else ""})
    out = pd.DataFrame(rows).sort_values("gain", ascending=False)
    out.to_csv(os.path.join(HERE, "biotransform.csv"), index=False)
    ng = int((out.gain > 0).sum())
    print(f"leads processed: {len(out)}")
    print(f"metabolites that GAIN polypharmacology vs parent: {ng} ({100*ng/len(out):.0f}%)")
    print(out.head(15)[["inchikey","parent_polypharm","best_metab_polypharm","gain"]].to_string(index=False))
    assert len(out) > 0, "no leads metabolised"
    print("check ok -> biotransform.csv")


if __name__ == "__main__":
    run()
