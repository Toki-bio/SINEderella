#!/usr/bin/env python3
"""composite_units.py RUN_ROOT [THREADS=4]
   composite_units.py --bank BANK.fa OUT.tsv [THREADS=4]

Which units each consensus of the searched bank is made of, for the report's Element hierarchy (owner 2026-10-18: show
the parts of composites that no unit explains, in grey; r4 sits inside P1 although P1 is named r1 + r3).

All consensuses are compared with all (ssearch36, E <= 1e-5, local hits >= MIN_LEN 50 bp at >= MIN_ID 75 %).
A consensus is a COMPOSITE when it holds a near-complete copy of a shorter consensus (>= 70 % of that one's length at
>= 90 % identity) and >= 120 bp beyond it; every other consensus is a UNIT. Each composite is then described by units only: at every position
the unit with the highest-identity hit covering it; stretches covered by no unit are "-" (drawn grey). Runs shorter than
10 bp join the run before them.

Out: <step2_output>/composite_units.tsv
  consensus  kind  start  end  unit  identity  alternatives
  (kind = unit: one row, the whole unit; kind = composite: one row per stretch; alternatives = the next best units over
   that stretch with their identity, e.g. "r4_32seqs 97.8")
"""
import glob
import os
import subprocess
import sys
import tempfile

MIN_ID, MIN_LEN, MIN_RUN = 75.0, 50, 10
FULL_ID, FULL_FRAC, MIN_REST = 90.0, 0.70, 120

if len(sys.argv) > 3 and sys.argv[1] == "--bank":
    # composite_units.py --bank BANK.fa OUT.tsv [THREADS]: the same table from a bank alone (step 1 uses it to search the
    # composites' stretches as queries of their own, tools/composite_pieces.py)
    bank, outp = sys.argv[2], sys.argv[3]
    threads = sys.argv[4] if len(sys.argv) > 4 else "4"
else:
    run = sys.argv[1].rstrip("/")
    threads = sys.argv[2] if len(sys.argv) > 2 else "4"
    s2s = sorted(glob.glob(os.path.join(run, "step2", "step2_output*")), key=os.path.getmtime)
    if not s2s:
        sys.exit("composite_units: no step2/step2_output* in %s" % run)
    bank = os.path.join(run, "consensuses.clean.fa")
    outp = os.path.join(s2s[-1], "composite_units.tsv")
if not os.path.isfile(bank):
    sys.exit("composite_units: missing %s" % bank)

seq, order, name = {}, [], None
for line in open(bank):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        seq[name] = []
        order.append(name)
    elif name:
        seq[name].append(line)
L = {k: len("".join(v)) for k, v in seq.items()}

with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
    with open(os.path.join(tmp, "b.fa"), "w") as fh:
        for k in order:
            fh.write(">%s\n%s\n" % (k, "".join(seq[k])))
    out = subprocess.run(["ssearch36", "-3", "-m", "8", "-E", "1e-5", "-T", threads, "b.fa", "b.fa"], cwd=tmp,
                         capture_output=True, text=True, check=True).stdout

hits = {}                                      # target -> [(start, end, query, identity)]
for line in out.splitlines():
    f = line.split("\t")
    q, t, pid, ln = f[0], f[1], float(f[2]), int(f[3])
    if q == t or pid < MIN_ID or ln < MIN_LEN or L[q] >= L[t]:
        continue
    a, b = sorted((int(f[8]), int(f[9])))
    hits.setdefault(t, []).append((a, b, q, pid))


def is_composite(t):
    """t holds a near-complete copy of a shorter consensus (>= FULL_FRAC of that consensus's length at >= FULL_ID %) and
    >= MIN_REST bp beyond it. Related monomers match each other in pieces at 75-95 %, so a looser rule (two partial hits)
    called rsi r1 and r3 composites; a length version (r3 = r4 + a 63 bp tail) is not a composite either."""
    for a, b, q, pid in hits.get(t, []):
        if pid >= FULL_ID and (b - a + 1) >= FULL_FRAC * L[q] and L[t] - (b - a + 1) >= MIN_REST:
            return True
    return False


comp = {t for t in order if is_composite(t)}
units = [t for t in order if t not in comp]
rows = []
for t in order:
    if t not in comp:
        rows.append((t, "unit", 1, L[t], t, "100.0", "-"))
        continue
    uh = [x for x in hits.get(t, []) if x[2] in units]
    best = []
    for p in range(1, L[t] + 1):
        cov = [(x[3], x[2]) for x in uh if x[0] <= p <= x[1]]
        best.append(max(cov) if cov else (0.0, "-"))
    runs = []                                  # [unit, start, end]
    for p, (pid, u) in enumerate(best, 1):
        if runs and runs[-1][0] == u:
            runs[-1][2] = p
        else:
            runs.append([u, p, p])
    merged = []
    for r in runs:
        if merged and r[2] - r[1] + 1 < MIN_RUN:
            merged[-1][2] = r[2]
        else:
            merged.append(r)
    for u, a, b in merged:
        mid = (a + b) // 2
        ids = [x[3] for x in uh if x[2] == u and x[0] <= mid <= x[1]]
        alts = sorted({(x[3], x[2]) for x in uh if x[0] <= mid <= x[1] and x[2] != u}, reverse=True)[:2]
        rows.append((t, "composite", a, b, u, "%.1f" % max(ids) if ids else "-",
                     ", ".join("%s %.1f" % (n, p) for p, n in alts) or "-"))

with open(outp, "w") as fo:
    fo.write("consensus\tkind\tstart\tend\tunit\tidentity\talternatives\n")
    for r in rows:
        fo.write("\t".join(str(x) for x in r) + "\n")
print("COMPOSITE_UNITS units=%d composites=%d -> %s" % (len(units), len(comp), outp))
