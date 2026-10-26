#!/usr/bin/env python3
"""fs9_bodycheck.py GENOME.fa COPIES.bed FS9_OUT_DIR [--random 20000] [--seed 1] [--alpha 0.05] [--min-copies 20]

Element (body) test of the twin pairs fs9_twins.sh found from the flanks; run by fs9_twins.sh at the end (BODYCHECK=1).

Copies multiplied together with their neighbourhood (segmental duplication, array unit) carry SINE bodies as similar as
their flanks. Two independent insertions next to similar sequence (a hot spot: two copies of a young repeat, or of another
SINE) carry bodies only as similar as two random copies of the family. A flank pair stays a twin (or array) pair only if
its two bodies are more similar than the 1 - alpha/k quantile of --random random pairs of the family's copies that are in
no flank pair (seeded; fewer than --min-background such copies: no check, the flank result stands), k = the
larger number of flank partners of its two copies: a copy at a hot spot has many partners and so many chances to pass by
luck (toy: 9 of 40 hot-spot copies passed a plain 95 % cut; with the partner correction 0). The other pairs get class
"flank_only" and do not make a copy a twin. Families with fewer than --min-copies copies are not checked.

Measured 2026-10-14 (flankscan/tests/real_twins/, share of pairs with the body above the random 95th percentile):
  Python bivittatus Squam3C, 6,713 copies: strict-pass pairs 79 %, pairs only found with masked bases kept (SOFTMASK=0) 4 %
  (870 pairs).  Anolis carolinensis Squam1, 79,225 copies: strict pairs 87 %, SOFTMASK=0-only pairs 18 % (151,502 pairs).

In place: flank_groups.tsv gets columns body_identity, body_alnlen and the class flank_only; copy_status.tsv gets the copy
status after the check (a twin or array copy that lost every pair falls back to unique / masked / untestable) and a fifth
column body_check (kept / flank_only / -). The flank-only versions stay as flank_groups.flank.tsv and copy_status.flank.tsv.
bodycheck.txt: background, threshold, counts. Needs bedtools and ssearch36. One ssearch36 call per copy with partners.
"""
import argparse, os, random, shutil, subprocess, sys, tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(description="element test of fs9 twin pairs")
ap.add_argument("genome"); ap.add_argument("bed"); ap.add_argument("out")
# 20000 random pairs (was 2000 until 2026-10-19): the threshold is a far-tail quantile (1 - alpha/k), and with 2000 pairs the
# seed alone changed 11-26 of ~18,000 rsi r1_r3_P1 verdicts; with 20000, 2-7 and the same threshold for every seed. Runtime
# about 3.3x (r1_r3_P1: 108 -> 353 s on therioserver).
ap.add_argument("--random", type=int, default=20000); ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--alpha", type=float, default=0.05); ap.add_argument("--min-copies", type=int, default=20)
ap.add_argument("--min-background", type=int, default=10)
# parallel ssearch36 calls (2026-10-24): one call per copy with partners, run one after another, took 8-10 h on tbr VES
# (620,807 copies, 99,793 with partners); results are keyed by pair, so they do not depend on the order of the calls
ap.add_argument("--threads", type=int, default=int(os.environ.get("BODYCHECK_THREADS") or 1))
a = ap.parse_args()
O = lambda f: os.path.join(a.out, f)
SST = os.environ.get("BODYCHECK_SSEARCH_T", "1").strip()       # default 1 since 2026-10-26; 0 = no -T (one worker per CPU)
if not SST.isdigit(): sys.exit("fs9_bodycheck: BODYCHECK_SSEARCH_T must be a number of threads (0 = no -T)")
if SST == "0": SST = ""

loc = {}
for line in open(a.bed):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 6: loc[f[3]] = (f[0], int(f[1]), int(f[2]), f[5])
fg = [l.rstrip("\n").split("\t") for l in open(O("flank_groups.tsv"))]
head, rows = fg[0][:7], [r[:7] for r in fg[1:]]
status = [l.rstrip("\n").split("\t")[:4] for l in open(O("copy_status.tsv"))]
names = sorted(loc)
if len(names) < a.min_copies or not rows:
    with open(O("bodycheck.txt"), "w") as fh:
        fh.write(f"skipped\t{'fewer than %d copies' % a.min_copies if len(names) < a.min_copies else 'no pairs'}\n")
    sys.exit(0)

# The background is drawn only from copies in no flank pair. Until 2026-10-19 it was drawn from all copies, so in a family
# made mostly of duplicated copies the background was the duplication identity itself, the duplicated pairs could not exceed
# it, and they were turned back into unique copies (GLM audit 2026-10-19).
paired = {r[0] for r in rows} | {r[1] for r in rows}
bgnames = [c for c in names if c not in paired]
if len(bgnames) < a.min_background:
    with open(O("bodycheck.txt"), "w") as fh:
        fh.write("skipped\tfewer than %d copies outside flank pairs (%d): flank result kept\n" % (a.min_background, len(bgnames)))
    sys.exit(0)
rnd = random.Random(a.seed)
rpairs = [tuple(rnd.sample(bgnames, 2)) for _ in range(a.random)]
need = sorted({r[0] for r in rows} | {r[1] for r in rows} | {x for p in rpairs for x in p})
tmp = tempfile.mkdtemp(prefix="bodycheck_", dir=a.out)
with open(os.path.join(tmp, "r.bed"), "w") as fh:
    for c in need:
        ch, s, e, st = loc[c]; fh.write(f"{ch}\t{s}\t{e}\t{c}\t0\t{st}\n")
fa = subprocess.run(["bedtools", "getfasta", "-fi", a.genome, "-bed", os.path.join(tmp, "r.bed"), "-s", "-nameOnly"],
                    capture_output=True, text=True, check=True).stdout
seq, name = {}, None
for l in fa.splitlines():
    if l.startswith(">"): name = l[1:].rsplit("(", 1)[0]; seq[name] = []
    else: seq[name].append(l.strip().upper())
seq = {k: "".join(v) for k, v in seq.items()}

def body_ids(pairs):
    """{(x, y): (identity, alnlen)} for unordered pairs, one ssearch36 call per first copy (its partners as the library)"""
    part = defaultdict(set)
    for x, y in pairs:
        x, y = (x, y) if x < y else (y, x)
        part[x].add(y)
    def one(job):
        i, x, ys = job
        q, lib = os.path.join(tmp, "q%d.fa" % i), os.path.join(tmp, "lib%d.fa" % i)
        with open(q, "w") as fh: fh.write(f">{x}\n{seq[x]}\n")
        with open(lib, "w") as fh:
            for y in sorted(ys): fh.write(f">{y}\n{seq[y]}\n")
        # BODYCHECK_SSEARCH_T (owner's decisions: opt-in 2026-10-24, default 1 since 2026-10-26; 0 = no -T = one worker per
        # CPU): a library of a few copies is scored against shuffles whose buffers are sized by the worker count, so the
        # E-values, and which weak pairs pass -E 100, depend on -T (docs/SPEEDUP_SSEARCH.md). -T 1: 35 ms CPU per one-partner
        # call instead of 250 ms, the same verdicts in every family tested, weak background pairs can differ, and the result no
        # longer depends on the machine's CPU count. Without -T, 32 parallel calls each started 128 workers on therioserver
        # (load ~1,400 during the 2026-10-26 runs)
        r = subprocess.run(["ssearch36"] + (["-T", SST] if SST else []) + ["-m", "8", "-E", "100", "-b", str(len(ys) + 5), "-d", "0", q, lib],
                           capture_output=True, text=True)
        os.remove(q); os.remove(lib)
        if r.returncode != 0: raise RuntimeError("fs9_bodycheck: ssearch36 failed: " + r.stderr[-300:])
        out = {}
        for l in r.stdout.splitlines():
            if not l.strip() or l.startswith("#"): continue
            f = l.split("\t")
            if f[1] in ys and (x, f[1]) not in out: out[(x, f[1])] = (float(f[2]), int(f[3]))   # first HSP = best
        for y in ys: out.setdefault((x, y), (0.0, 0))
        return out
    res = {}
    jobs = [(i, x, ys) for i, (x, ys) in enumerate(part.items())]
    try:
        with ThreadPoolExecutor(max_workers=max(1, a.threads)) as ex:
            for out in ex.map(one, jobs):
                res.update(out)
    except RuntimeError as e:
        sys.exit(str(e))
    return res

def get(res, x, y): return res[(x, y) if x < y else (y, x)]
rres = body_ids(rpairs)
bg = sorted(get(rres, x, y)[0] for x, y in rpairs)
def quant(q): return bg[min(len(bg) - 1, int(q * len(bg)))]
thr = quant(1 - a.alpha)
deg = defaultdict(set)
for r in rows: deg[r[0]].add(r[1]); deg[r[1]].add(r[0])
def pair_thr(x, y): return quant(1 - a.alpha / max(1, len(deg[x]), len(deg[y])))
pres = body_ids([(r[0], r[1]) for r in rows])

shutil.copy(O("flank_groups.tsv"), O("flank_groups.flank.tsv")); shutil.copy(O("copy_status.tsv"), O("copy_status.flank.tsv"))
best = {}; nfail = 0
with open(O("flank_groups.tsv"), "w") as fh:
    fh.write("\t".join(head + ["body_identity", "body_alnlen"]) + "\n")
    for r in rows:
        bi, bl = get(pres, r[0], r[1])
        cl = r[6] if bi > pair_thr(r[0], r[1]) else "flank_only"
        if cl == "flank_only": nfail += 1
        else:
            st = "array" if cl == "array" else ("twin1" if r[5] == "T1" else "twin2")
            rank = {"twin1": 4, "twin2": 3, "array": 2}[st]
            for c in (r[0], r[1]):
                if rank > best.get(c, (0, ""))[0]: best[c] = (rank, st)
        fh.write("\t".join(r[:6] + [cl, "%.1f" % bi, str(bl)]) + "\n")
fallback = {c: "unique" for c in names}
for fn, st in (("masked.txt", "masked"), ("untestable.txt", "untestable")):    # untestable wins over masked, as in fs9_twins.sh
    if os.path.exists(O(fn)):
        for l in open(O(fn)): fallback[l.strip()] = st
counts = defaultdict(int)
with open(O("copy_status.tsv"), "w") as fh:
    for s in status:
        c, st = s[0], s[1]
        if st in ("twin1", "twin2", "array"):
            new, chk = (best[c][1], "kept") if c in best else (fallback.get(c, "unique"), "flank_only")
            grp = s[3] if c in best else "-"
        else:
            new, chk, grp = st, "-", s[3]
        counts[new] += 1
        fh.write("\t".join([c, new, s[2], grp, chk]) + "\n")
with open(O("bodycheck.txt"), "w") as fh:
    fh.write(f"random_pairs\t{len(bg)}\nrandom_median\t{bg[len(bg) // 2]:.1f}\nthreshold_one_partner\t{thr:.1f}\nthreshold_max\t{bg[-1]:.1f}\n")
    fh.write(f"pairs\t{len(rows)}\nflank_only\t{nfail}\n")
    for k in sorted(counts): fh.write(f"{k}\t{counts[k]}\n")
shutil.rmtree(tmp)
print(open(O("bodycheck.txt")).read(), end="")
