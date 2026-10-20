#!/usr/bin/env python
"""
NP-BioSynAct -- ADMET / drug-likeness profiling of AD-inside multi-target leads.
RDKit only. MW, logP, HBD/HBA, TPSA, RotB, QED, SAscore, Lipinski/Veber,
PAINS flag, BBB (BOILED-Egg approx). Answers "are the leads drug-like & CNS-ready?"

Run: conda activate npbs && python 51_admet.py
Out: admet_leads.csv
"""
import os, sys, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, QED, Crippen, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score")); import sascorer
RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

_p = FilterCatalogParams()
for c in ("PAINS_A", "PAINS_B", "PAINS_C"):
    _p.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, c))
PAINS = FilterCatalog(_p)


def profile(smi):
    m = Chem.MolFromSmiles(smi)
    if not m:
        return None
    mw, logp, tpsa = Descriptors.MolWt(m), Crippen.MolLogP(m), rdMolDescriptors.CalcTPSA(m)
    hbd, hba = rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(m)
    lip = sum([mw > 500, logp > 5, hbd > 5, hba > 10])          # Lipinski violations
    veber = rotb <= 10 and tpsa <= 140
    bbb = (tpsa <= 79 and 0.4 <= logp <= 6.0)                   # BOILED-Egg yolk approx
    return dict(MW=round(mw,1), logP=round(logp,2), HBD=hbd, HBA=hba, TPSA=round(tpsa,1),
                RotB=rotb, QED=round(QED.qed(m),3), SAscore=round(sascorer.calculateScore(m),2),
                Lipinski_viol=lip, Veber=veber, BBB=bbb, PAINS=PAINS.HasMatch(m))


def run():
    leads = pd.read_csv(os.path.join(HERE, "atlas_leads_AD.csv")).head(500)
    smi = pd.read_csv(os.path.join(HERE, "np_compounds.csv"),
                      usecols=["inchikey", "smiles"], low_memory=False).set_index("inchikey").smiles
    leads["smiles"] = leads.inchikey.map(smi)
    rows = []
    for _, r in leads.dropna(subset=["smiles"]).iterrows():
        p = profile(r.smiles)
        if p:
            p.update(inchikey=r.inchikey, bio_class=r.get("bio_class"),
                     n_hits=r.get("n_hits"), mean_prob=r.get("mean_prob"))
            rows.append(p)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "admet_leads.csv"), index=False)
    n = len(out)
    print(f"profiled {n} AD-inside leads")
    print(f"  drug-like (Lipinski<=1, Veber, not PAINS): {int(((out.Lipinski_viol<=1)&out.Veber&~out.PAINS).sum())} ({100*((out.Lipinski_viol<=1)&out.Veber&~out.PAINS).mean():.0f}%)")
    print(f"  BBB-permeant (CNS-ready): {int(out.BBB.sum())} ({100*out.BBB.mean():.0f}%)")
    print(f"  PAINS-flagged: {int(out.PAINS.sum())} ({100*out.PAINS.mean():.0f}%)")
    print(f"  median QED {out.QED.median():.2f}  median SAscore {out.SAscore.median():.2f}")
    assert n > 0, "no leads profiled"
    print("check ok -> admet_leads.csv")


if __name__ == "__main__":
    run()
