#!/usr/bin/env bash
# MM-GBSA (gmx_MMPBSA) + stability analysis, LAYOUT-AGNOSTIC.
# Reads md_inputs.csv (fill from your StreaMD outputs), one row per trajectory:
#   target,tpr,xtc,topol,ligres
#   DPP4,path/md.tpr,path/md.xtc,path/topol.top,LIG
#   MAPK14_r1,...,LIG
# ligres = the ligand residue name in the topology (check: grep -A2 '\[ molecules' topol.top).
# Needs: gmx, gmx_MMPBSA.
set -euo pipefail
cd "$(dirname "$0")"
[ -f md_inputs.csv ] || { cat > md_inputs.csv <<'C'
target,tpr,xtc,topol,ligres
DPP4,,,,
MAPK14,,,,
NTRK1,,,,
OPRM1,,,,
PTGS2,,,,
C
echo "md_inputs.csv template written -- fill tpr/xtc/topol/ligres from your StreaMD run, then re-run"; exit 0; }

mkdir -p mmpbsa
cat > mmpbsa/mmpbsa.in <<'M'
&general
  sys_name="NPBS", startframe=1, interval=5, verbose=2,
/
&gb
  igb=5, saltcon=0.150,
/
M

tail -n +2 md_inputs.csv | while IFS=, read -r T TPR XTC TOP LIG; do
  [ -z "${T:-}" ] && continue
  [ -f "${TPR:-}" ] && [ -f "${XTC:-}" ] || { echo "  $T: tpr/xtc missing, skip"; continue; }
  echo "== $T (lig $LIG) =="
  # index: Protein (1) + ligand group named LIG (20)
  printf "r %s\nname 20 LIG\nq\n" "$LIG" | gmx make_ndx -f "$TPR" -o mmpbsa/${T}.ndx >/dev/null 2>&1
  printf "Protein\nSystem\n" | gmx trjconv -s "$TPR" -f "$XTC" -o mmpbsa/${T}_c.xtc -center -pbc mol -ur compact >/dev/null 2>&1
  ( cd mmpbsa && gmx_MMPBSA -O -i mmpbsa.in -cs "../$TPR" -ci ${T}.ndx -cg 1 20 \
        -ct ${T}_c.xtc -cp "../$TOP" -o ${T}_MMPBSA.dat -eo ${T}_MMPBSA.csv >/dev/null 2>&1 ) \
        || echo "  mmpbsa failed ($T) -- verify -cg group numbers vs ${T}.ndx"
  # stability
  printf "Backbone\nBackbone\n" | gmx rms   -s "$TPR" -f mmpbsa/${T}_c.xtc -o mmpbsa/${T}_rmsd.xvg >/dev/null 2>&1
  printf "Protein\n"            | gmx rmsf  -s "$TPR" -f mmpbsa/${T}_c.xtc -o mmpbsa/${T}_rmsf.xvg >/dev/null 2>&1
  printf "Protein\n"            | gmx gyrate -s "$TPR" -f mmpbsa/${T}_c.xtc -o mmpbsa/${T}_rg.xvg >/dev/null 2>&1
  printf "Protein\nLIG\n"       | gmx hbond -s "$TPR" -f mmpbsa/${T}_c.xtc -num mmpbsa/${T}_hb.xvg >/dev/null 2>&1 || true
  echo "  done: mmpbsa/${T}_MMPBSA.dat + rmsd/rmsf/rg/hb"
done
echo "all done -> mmpbsa/  (next: 58 collate deltaG vs ML prob)"
# ponytail: if 'name 20 LIG' clashes with existing groups, check gmx make_ndx numbering and adjust -cg.
