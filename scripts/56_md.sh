#!/usr/bin/env bash
# MD for the redock-valid complexes in md_shortlist.csv.
# Protein: amber99sb-ildn + TIP3P. Ligand: GAFF2 via acpype (bcc charges).
# Per complex x NREP replicas: build -> EM -> NVT -> NPT -> production.
# Ligand pose = best smina pose from dock/<T>_ligands_smina.sdf.
#
# vars: NS production length (ns), NREP replicas. Needs: gromacs (gmx), acpype, obabel, rdkit.
set -euo pipefail
cd "$(dirname "$0")"
NS=${NS:-100}; NREP=${NREP:-3}; FF=amber99sb-ildn
mkdir -p md mdp

# --- mdp files (CHARMM/AMBER-compatible, PME, V-rescale, Parrinello-Rahman) ---
cat > mdp/em.mdp <<'M'
integrator=steep
emtol=1000.0
nsteps=50000
nstlist=10
cutoff-scheme=Verlet
coulombtype=PME
rvdw=1.0
rcoulomb=1.0
M
cat > mdp/nvt.mdp <<'M'
integrator=md
nsteps=50000
dt=0.002
tcoupl=V-rescale
tc-grps=Protein Non-Protein
tau-t=0.1 0.1
ref-t=310 310
constraints=h-bonds
constraint-algorithm=lincs
cutoff-scheme=Verlet
coulombtype=PME
rvdw=1.0
rcoulomb=1.0
gen-vel=yes
gen-temp=310
define=-DPOSRES
M
cat > mdp/npt.mdp <<'M'
integrator=md
nsteps=50000
dt=0.002
tcoupl=V-rescale
tc-grps=Protein Non-Protein
tau-t=0.1 0.1
ref-t=310 310
pcoupl=Parrinello-Rahman
pcoupltype=isotropic
tau-p=2.0
ref-p=1.0
compressibility=4.5e-5
constraints=h-bonds
cutoff-scheme=Verlet
coulombtype=PME
rvdw=1.0
rcoulomb=1.0
define=-DPOSRES
M
NSTEPS=$(python -c "print(int($NS*1000/0.002))")
cat > mdp/md.mdp <<M
integrator=md
nsteps=$NSTEPS
dt=0.002
tcoupl=V-rescale
tc-grps=Protein Non-Protein
tau-t=0.1 0.1
ref-t=310 310
pcoupl=Parrinello-Rahman
pcoupltype=isotropic
tau-p=2.0
ref-p=1.0
compressibility=4.5e-5
constraints=h-bonds
cutoff-scheme=Verlet
coulombtype=PME
rvdw=1.0
rcoulomb=1.0
nstxout-compressed=5000
M

tail -n +2 md_shortlist.csv | while IFS=, read -r T CODE NAME SCORE; do
  echo "==== $T / $CODE ===="
  # 1. extract best NP pose -> mol2
  python - "$T" "$CODE" <<'PY'
import sys
from rdkit import Chem
t,code=sys.argv[1],sys.argv[2]
w=Chem.SDWriter(f"md/{t}_lig.sdf"); done=False
for m in Chem.SDMolSupplier(f"dock/{t}_ligands_smina.sdf",sanitize=False):
    if m is None: continue
    if (m.GetProp("_Name") if m.HasProp("_Name") else "")==code and not done:
        w.write(m); done=True
w.close(); print("pose extracted" if done else "!! pose not found")
PY
  obabel md/${T}_lig.sdf -O md/${T}_lig.mol2 -p 7.4 --partialcharge gasteiger >/dev/null 2>&1
  # 2. ligand params (GAFF2)
  ( cd md && acpype -i ${T}_lig.mol2 -b ${T}LIG -n 0 -a gaff2 -c bcc >/dev/null 2>&1 )
  LIGDIR=md/${T}LIG.acpype
  # 3. protein topology
  gmx pdb2gmx -f rec/${T}_rec.pdb -ff $FF -water tip3p -ignh \
      -o md/${T}_prot.gro -p md/${T}_topol.top -i md/${T}_posre.itp >/dev/null 2>&1
  # 4. merge protein + ligand coords + topology
  python - "$T" "$LIGDIR" <<'PY'
import sys
t,ld=sys.argv[1],sys.argv[2]
P=open(f"md/{t}_prot.gro").read().splitlines()
L=open(f"{ld}/{t}LIG_GMX.gro").read().splitlines()
pat=P[2:-1]; lat=L[2:-1]; box=P[-1]
open(f"md/{t}_complex.gro","w").write(
    "complex\n%d\n%s\n%s\n" % (len(pat)+len(lat), "\n".join(pat+lat), box))
# topology: include ligand itp after ff, add molecule
top=open(f"md/{t}_topol.top").read().splitlines(); out=[]
for l in top:
    out.append(l)
    if l.strip().endswith('forcefield.itp"'):
        out += [f'#include "{ld}/{t}LIG_GMX.itp"',
                '#ifdef POSRES', f'#include "{ld}/posre_{t}LIG.itp"', '#endif']
out.append(f"{t}LIG   1")
open(f"md/{t}_topol.top","w").write("\n".join(out)+"\n")
print("merged complex + topology")
PY
  cd md
  gmx editconf -f ${T}_complex.gro -o ${T}_box.gro -bt dodecahedron -d 1.0 >/dev/null 2>&1
  gmx solvate -cp ${T}_box.gro -cs spc216.gro -o ${T}_solv.gro -p ${T}_topol.top >/dev/null 2>&1
  gmx grompp -f ../mdp/em.mdp -c ${T}_solv.gro -p ${T}_topol.top -o ${T}_ion.tpr -maxwarn 2 >/dev/null 2>&1
  echo SOL | gmx genion -s ${T}_ion.tpr -o ${T}_ions.gro -p ${T}_topol.top -pname NA -nname CL -neutral -conc 0.15 >/dev/null 2>&1
  for R in $(seq 1 $NREP); do
    echo "  replica $R"
    gmx grompp -f ../mdp/em.mdp  -c ${T}_ions.gro -p ${T}_topol.top -o ${T}_em_$R.tpr -maxwarn 2 >/dev/null 2>&1
    gmx mdrun -deffnm ${T}_em_$R -nb gpu >/dev/null 2>&1
    gmx grompp -f ../mdp/nvt.mdp -c ${T}_em_$R.gro -r ${T}_em_$R.gro -p ${T}_topol.top -o ${T}_nvt_$R.tpr -maxwarn 2 >/dev/null 2>&1
    gmx mdrun -deffnm ${T}_nvt_$R -nb gpu >/dev/null 2>&1
    gmx grompp -f ../mdp/npt.mdp -c ${T}_nvt_$R.gro -r ${T}_nvt_$R.gro -t ${T}_nvt_$R.cpt -p ${T}_topol.top -o ${T}_npt_$R.tpr -maxwarn 2 >/dev/null 2>&1
    gmx mdrun -deffnm ${T}_npt_$R -nb gpu >/dev/null 2>&1
    gmx grompp -f ../mdp/md.mdp -c ${T}_npt_$R.gro -t ${T}_npt_$R.cpt -p ${T}_topol.top -o ${T}_md_$R.tpr -maxwarn 2 >/dev/null 2>&1
    gmx mdrun -deffnm ${T}_md_$R -nb gpu -bonded gpu
  done
  cd ..
done
echo "done -> md/<T>_md_<R>.xtc  (next: 57_mmpbsa.sh, then analysis)"
# ponytail: bcc charges + single protein FF; 3 replicas. tune protonation/box per system if needed.
