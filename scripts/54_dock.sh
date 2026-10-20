#!/usr/bin/env bash
# Full docking: auto-make targets.tsv from pairs.csv (+suggested PDBs), fetch PDB,
# AUTO-pick co-crystal ligand (largest non-water HETATM), split receptor, redock
# validate, dock leads+decoys with smina AND gnina. No manual ligand resnames.
#
# Override a PDB/chain by editing targets.tsv (TARGET<TAB>PDB<TAB>CHAIN) then re-run.
# Needs: wget, obabel, obrms, smina, gnina.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p dock rec

# decoys.sdf if missing
[ -f decoys.sdf ] || python - <<'PY'
import pandas as pd; from rdkit import Chem; from rdkit.Chem import AllChem
d=pd.read_csv("decoys.csv"); w=Chem.SDWriter("decoys.sdf")
for _,r in d.iterrows():
    m=Chem.MolFromSmiles(r.smiles)
    if not m: continue
    m=Chem.AddHs(m)
    if AllChem.EmbedMolecule(m,AllChem.ETKDGv3())!=0: continue
    AllChem.MMFFOptimizeMolecule(m); m.SetProp("_Name","decoy_"+str(r.inchikey)); w.write(m)
w.close(); print("decoys.sdf written")
PY

# auto-generate targets.tsv from pairs.csv unique targets + suggested holo PDBs
[ -f targets.tsv ] || python - <<'PY'
import pandas as pd
PDB={"OPRM1":"5C1M","CNR2":"6PT0","PPARG":"2PRG","DPP4":"1X70","PTGS2":"5KIR",
     "GSK3B":"1Q5K","MAPK14":"1A9U","PTPN1":"2QBS","AKR1B1":"1US0","TNF":"2AZ5",
     "ADRB2":"2RH1","NTRK1":"4YNE","SIRT1":"4IG9","NOS2":"3E7G","ADORA1":"5UEN","HTR1A":"7E2X"}
t=sorted(pd.read_csv("pairs.csv").target.unique())
with open("targets.tsv","w") as o:
    for x in t:
        if x in PDB: o.write(f"{x}\t{PDB[x]}\tA\n")
print("targets.tsv written (verify PDB/chain; edit if needed):")
print(open("targets.tsv").read())
PY

echo "target,ligand_set,engine,name,score" > dock_scores.csv
while IFS=$'\t' read -r T PDB CH; do
  [ -z "${T:-}" ] && continue; case "$T" in \#*) continue;; esac
  CH=${CH:-A}
  echo "== $T ($PDB chain $CH) =="
  wget -q "https://files.rcsb.org/download/${PDB}.pdb" -O rec/${T}.pdb || { echo "  fetch failed"; continue; }
  # split receptor (ATOM, chain) + auto-pick largest non-water HETATM ligand
  python - "$T" "$CH" <<'PY'
import sys,collections
t,ch=sys.argv[1],sys.argv[2]
SKIP={"HOH","WAT","NA","CL","K","MG","ZN","CA","MN","SO4","PO4","GOL","EDO","ACT",
      "DMS","PEG","BME","IOD","FMT","TRS","CIT","NAG","BMA","MAN"}
rec=[]; het=collections.defaultdict(list)
for L in open(f"rec/{t}.pdb"):
    if L.startswith("ATOM") and L[21]==ch: rec.append(L)
    elif L.startswith("HETATM") and L[21]==ch:
        rn=L[17:20].strip()
        if rn not in SKIP: het[(rn,L[22:26])].append(L)
open(f"rec/{t}_rec.pdb","w").writelines(rec)
if het:
    big=max(het,key=lambda k:len(het[k]))
    open(f"rec/{t}_lig.pdb","w").writelines(het[big])
    print(f"  ligand auto-picked: {big[0]} ({len(het[big])} atoms)")
else:
    open(f"rec/{t}_lig.pdb","w").write(""); print("  !! no HETATM ligand found -- set chain or PDB")
PY
  [ -s rec/${T}_lig.pdb ] || { echo "  skip $T (no ligand)"; continue; }
  obabel rec/${T}_rec.pdb -O rec/${T}_rec.pdbqt -xr -p 7.4 >/dev/null 2>&1
  obabel rec/${T}_lig.pdb -O rec/${T}_ref.sdf >/dev/null 2>&1
  smina -r rec/${T}_rec.pdbqt -l rec/${T}_ref.sdf --autobox_ligand rec/${T}_ref.sdf \
        --exhaustiveness 16 -o dock/${T}_redock.sdf >/dev/null 2>&1
  RMSD=$(obrms rec/${T}_ref.sdf dock/${T}_redock.sdf 2>/dev/null | awk '{print $NF}' | head -1)
  echo "  redock RMSD = ${RMSD:-NA} A (want <=2.0)"; echo "$T,REDOCK,smina,ref,${RMSD:-NA}" >> dock_scores.csv
  for SET in ligands decoys; do
    [ -f ${SET}.sdf ] || continue
    smina -r rec/${T}_rec.pdbqt -l ${SET}.sdf --autobox_ligand rec/${T}_ref.sdf \
          --exhaustiveness 16 --num_modes 5 -o dock/${T}_${SET}_smina.sdf >/dev/null 2>&1 || true
    gnina -r rec/${T}_rec.pdbqt -l ${SET}.sdf --autobox_ligand rec/${T}_ref.sdf \
          --num_modes 5 -o dock/${T}_${SET}_gnina.sdf >/dev/null 2>&1 || true
    python - "$T" "$SET" <<'PY'
import sys,os
from rdkit import Chem
t,s=sys.argv[1],sys.argv[2]
def dump(sdf,eng):
    if not os.path.exists(sdf): return
    seen={}
    for m in Chem.SDMolSupplier(sdf,sanitize=False):
        if m is None: continue
        nm=m.GetProp("_Name") if m.HasProp("_Name") else "?"
        for p in ("minimizedAffinity","CNNaffinity","CNNscore","affinity"):
            if m.HasProp(p):
                if nm not in seen: seen[nm]=m.GetProp(p)
                break
    open("dock_scores.csv","a").write("".join(f"{t},{s},{eng},{nm},{v}\n" for nm,v in seen.items()))
dump(f"dock/{t}_{s}_smina.sdf","smina"); dump(f"dock/{t}_{s}_gnina.sdf","gnina")
print(f"  {t} {s}: scored")
PY
  done
done < targets.tsv
echo "done -> dock_scores.csv (check redock<=2A; leads should beat decoys)"
