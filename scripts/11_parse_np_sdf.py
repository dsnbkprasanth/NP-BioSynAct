#!/usr/bin/env python
"""
NP-BioSynAct -- parse COCONUT + LOTUS + IMPPAT SDFs into one training table.
Streams (never loads a 2GB file into RAM), dedups by InChIKey.

Out: np_compounds.csv  ->  id, source, smiles, inchikey, murcko_scaffold,
                            bio_pathway, bio_superclass, bio_class

Run: conda activate npbs && python 11_parse_np_sdf.py
"""
import os, re, csv, sys
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
SDFS = {  # source: filename
    "COCONUT": "coconut_sdf_3d-07-2026.sdf",
    "LOTUS":   "acd_lotusv9.sdf",
    "IMPPAT":  "imppat_phytocompounds.sdf",
}
OUT = os.path.join(HERE, "np_compounds.csv")

# NPClassifier tag names differ per source; match loosely.
PATH_RE  = re.compile(r"(np.?classifier.*pathway|^pathway$)", re.I)
SUPER_RE = re.compile(r"(np.?classifier.*superclass|^superclass$)", re.I)
CLASS_RE = re.compile(r"(np.?classifier.*class$|^class$)", re.I)


def pick(props, rx):
    for k in props:
        if rx.search(k):
            return props[k]
    return ""


def run():
    seen = set()
    n_in = n_out = 0
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "source", "smiles", "inchikey", "murcko_scaffold",
                    "bio_pathway", "bio_superclass", "bio_class"])
        for source, fname in SDFS.items():
            path = os.path.join(HERE, fname)
            if not os.path.exists(path):
                print(f"  skip {source}: {fname} not found"); continue
            print(f"  parsing {source} ...")
            supp = Chem.ForwardSDMolSupplier(path, sanitize=True)
            for mol in supp:
                n_in += 1
                if mol is None:
                    continue
                try:
                    ik = Chem.MolToInchiKey(mol)
                except Exception:
                    continue
                if not ik or ik in seen:
                    continue
                seen.add(ik)
                props = {k: mol.GetProp(k) for k in mol.GetPropNames()}
                smi = Chem.MolToSmiles(mol)
                try:
                    scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
                except Exception:
                    scaf = ""
                cid = (mol.GetProp("_Name") or "").strip() or ik
                w.writerow([cid, source, smi, ik, scaf,
                            pick(props, PATH_RE), pick(props, SUPER_RE),
                            pick(props, CLASS_RE)])
                n_out += 1
                if n_out % 50000 == 0:
                    print(f"    {n_out} unique written")
    print(f"done: {n_in} read, {n_out} unique -> {OUT}")
    # self-check: file has a header + rows, every row has a SMILES
    with open(OUT) as fh:
        rows = list(csv.DictReader(fh))
    assert rows and all(r["smiles"] for r in rows), "empty or missing SMILES"
    print(f"check ok: {len(rows)} rows, all have SMILES")
    # how many carry biosynthesis labels (COCONUT should; LOTUS/IMPPAT may not)
    lab = sum(1 for r in rows if r["bio_pathway"])
    print(f"biosynthesis-labelled: {lab}/{len(rows)}"
          f" ({100*lab/len(rows):.0f}%)")


if __name__ == "__main__":
    run()
