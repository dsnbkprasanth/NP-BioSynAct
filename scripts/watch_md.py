#!/usr/bin/env python
"""Live MD progress across MD/*/ (StreaMD/GROMACS). ASCII bars, no deps.
Shows the newest running stage (em/nvt/npt/md) per complex + % done.
Run: python watch_md.py            (from the repository root)
Ctrl-C to quit.
"""
import glob, os, re, time, sys
ROOT = sys.argv[1] if len(sys.argv) > 1 else "MD"
STAGES = ["em", "nvt", "npt", "md_out", "md"]          # order = pipeline order

def parse(logf):
    """(current_step, total_nsteps) from a GROMACS .log, else (None,None)."""
    try: txt = open(logf, errors="ignore").read()
    except OSError: return None, None
    tot = re.search(r"nsteps\s*=\s*(\d+)", txt)
    tot = int(tot.group(1)) if tot else None
    steps = re.findall(r"^\s*(\d+)\s+[\d.]+\s*$", txt, re.M)   # "Step Time" rows
    cur = int(steps[-1]) if steps else (tot if "Finished mdrun" in txt else 0)
    return cur, tot

def bar(frac, n=30):
    frac = max(0.0, min(1.0, frac)); f = int(frac*n)
    return "[" + "#"*f + "-"*(n-f) + f"] {frac*100:5.1f}%"

def stage_of(d):
    logs = glob.glob(f"{d}/**/*.log", recursive=True)
    logs = [l for l in logs if os.path.basename(l)[:-4] in STAGES or "md_out" in l]
    if not logs: return None, None, None
    newest = max(logs, key=os.path.getmtime)
    name = "md" if "md_out" in newest or os.path.basename(newest)=="md.log" else os.path.basename(newest)[:-4]
    cur, tot = parse(newest)
    return name, cur, tot

def main():
    dirs = sorted(d for d in glob.glob(f"{ROOT}/*") if os.path.isdir(d))
    if not dirs: sys.exit(f"no complexes under {ROOT}/")
    while True:
        out = [f"MD progress  {time.strftime('%H:%M:%S')}   ({len(dirs)} complexes)\n"]
        for d in dirs:
            name, cur, tot = stage_of(d)
            tag = os.path.basename(d)
            if name is None:
                out.append(f"{tag:16} waiting (prep)")
            elif tot:
                out.append(f"{tag:16} {name:6} {bar(cur/tot)}  {cur}/{tot}")
            else:
                out.append(f"{tag:16} {name:6} step {cur}")
        os.system("clear")
        print("\n".join(out))
        print("\nCtrl-C to quit. equilibration = em/nvt/npt; production = md.")
        time.sleep(15)

if __name__ == "__main__":
    try: main()
    except KeyboardInterrupt: print()
