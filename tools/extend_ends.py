#!/usr/bin/env python3
"""extend_ends.py - extend merged loci whose element goes on past the end of the hits (step 1, before extraction).

    extend_ends.py MERGED_BED LABELED_BED GENOME OUT_BED OUT_DIR [--win 100] [--k 10] [--step 10] [--min-n 20]
                   [--min-frac 0.02] [--share 0.2] [--min-shared 10] [--bridge 2] [--flag-frac 0.2]

Why. A locus is the union of the hits of all bank consensuses, so it ends where the longest matching consensus stops. When the
genome's element is longer than every bank consensus that hit it, the copy is cut there. Found on hla, 2026-10-29: shorties of
134 bp (r5-type head + linker + GGCAGGATGGCTCAGCT + A tail) were cut before their tail in 764 of 1,400 and 708 of 2,350 genome
copies, but kept it where r5_r3_P34 also hit (its second module matches the tail). On the SubFam plate the tail then looked like
a presence/absence variant. The continuation check (tools/continuation_check.py) does not see this: it runs after step 2, only
finds continuations that match another bank consensus, and changes nothing.

How. Loci are grouped by the label of their highest-scoring hit. For every locus the WIN bp past each end are read in the
element's orientation (never past the neighbouring locus or the contig end). Per group and side, every K-mer is counted once per
window it occurs in; K-mers made only of A/T, or with >= K-3 copies of one base, are not counted (A tails, microsatellites).
A K-mer is "shared" when it occurs in at least max(MIN_N, MIN_FRAC x windows of the group) windows: unique flank gives K-mers seen
once or twice, an untaken part of the element gives the same K-mers in hundreds of copies. A group with fewer than MIN_SHARED
shared K-mers extends nothing. K = 10 and SHARE = 0.2 because old copies keep few intact K-mers (15 % divergence: 0.85^10 = 20 %). Each end is walked outward in STEP-bp
steps; a step is element when it has at least 3 counted K-mers and at least SHARE of them are shared, a miss otherwise. The walk
stops at the (BRIDGE+1)-th miss in a row (up to BRIDGE misses are bridged: a single mutation removes every K-mer overlapping
it, a whole step). The new end is the far end of the last shared K-mer; when the walk extended the end, a following A run (>= 80 % A per 10 bp)
is added as the A tail. A walk still in element at WIN is not an end (array, dimer, longer element): that locus is left as it
is and counted (at_cap_not_extended).

Tested 2026-10-29 (docs/END_EXTENSION.md): hla shorties, rsi control.

Out: OUT_BED (same columns as MERGED_BED, ends moved), OUT_DIR/ends_extension.tsv (per group and side: loci, extended, median and
max extension, flag when >= FLAG_FRAC of the group grew by >= 10 bp at that side: the bank consensuses are shorter than the
element there), OUT_DIR/loci_extended.tsv (every locus that moved: old and new coordinates, group, bp added 5'/3').
Prints one line: EXTEND loci=.. extended5=.. extended3=.. flagged=..
"""
import argparse
import collections
import os
import statistics
import subprocess
import sys
import tempfile

import numpy as np

CODE = np.full(256, 255, dtype=np.uint8)
for i, c in enumerate("ACGT"):
    CODE[ord(c)] = i
    CODE[ord(c.lower())] = i


def kmer_ints(seq, k):
    """K-mer codes of seq (index j = K-mer starting at j); -1 where the K-mer holds a non-ACGT base or is not counted."""
    a = CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(a) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.int64)
    bad = (a == 255)
    a = np.where(bad, 0, a).astype(np.int64)
    w = np.lib.stride_tricks.sliding_window_view(a, k)
    codes = (w * (4 ** np.arange(k - 1, -1, -1, dtype=np.int64))).sum(axis=1)
    badk = np.lib.stride_tricks.sliding_window_view(bad, k).any(axis=1)
    cnt = np.stack([(w == b).sum(axis=1) for b in range(4)], axis=1)
    lowc = (cnt[:, 0] + cnt[:, 3] == k) | (cnt.max(axis=1) >= k - 3)
    codes[badk | lowc] = -1
    return codes


def walk(codes, shared, k, step, share, win, bridge=1):
    """Bases to add: codes are the K-mers of the window read outward (index 0 next to the locus end)."""
    # one step may fail between element steps: a single mutation removes all K-mers that overlap it, i.e. a whole step
    last, miss, s = 0, 0, 0
    while s * step < len(codes):
        c = codes[s * step:(s + 1) * step]
        c = c[c >= 0]
        if len(c) < 3 or np.isin(c, shared).mean() < share:
            miss += 1
            if miss > bridge:
                break
            s += 1
            continue
        miss = 0
        idx = np.nonzero(np.isin(codes[s * step:(s + 1) * step], shared))[0]
        last = s * step + int(idx[-1]) + k
        s += 1
    return min(last, win)


def a_run(seq, start, win):
    """Length of the A run at seq[start:] (>= 80 % A in each 10-bp step, then base by base), capped so start + run <= win."""
    p = start
    while p + 10 <= len(seq) and seq[p:p + 10].upper().count("A") >= 8 and p + 10 <= win:
        p += 10
    while p < len(seq) and p < win and seq[p].upper() == "A":
        p += 1
    return p - start


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("merged"); ap.add_argument("labeled"); ap.add_argument("genome"); ap.add_argument("out_bed"); ap.add_argument("out_dir")
    ap.add_argument("--win", type=int, default=100); ap.add_argument("--k", type=int, default=10); ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--min-n", type=int, default=20); ap.add_argument("--min-frac", type=float, default=0.02)
    ap.add_argument("--share", type=float, default=0.2); ap.add_argument("--flag-frac", type=float, default=0.2)
    ap.add_argument("--min-shared", type=int, default=10); ap.add_argument("--bridge", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    size = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.genome + ".fai")}
    loci = [l.rstrip("\n").split("\t") for l in open(a.merged) if l.strip()]
    if not loci:
        open(a.out_bed, "w").close(); print("EXTEND loci=0 extended5=0 extended3=0 flagged=0"); return
    # group = label of the highest-scoring labeled hit on the locus' strand
    tmp = tempfile.mkdtemp(prefix="extend_", dir=a.out_dir)
    with open(tmp + "/loci.bed", "w") as fh:
        for i, f in enumerate(loci):
            fh.write("%s\t%s\t%s\t%d\t0\t%s\n" % (f[0], f[1], f[2], i, f[5]))
    lab = {}
    p = subprocess.run("bedtools intersect -s -wa -wb -a %s/loci.bed -b %s" % (tmp, a.labeled), shell=True, check=True,
                       stdout=subprocess.PIPE, universal_newlines=True)
    for l in p.stdout.splitlines():
        f = l.split("\t")
        i, sc = int(f[3]), float(f[10]) if f[10] not in ("", ".") else 0.0
        if i not in lab or sc > lab[i][0]:
            lab[i] = (sc, f[9])
    group = [lab[i][1] if i in lab else "-" for i in range(len(loci))]
    # room on each side: up to the neighbouring locus (any strand) and the contig ends
    order = sorted(range(len(loci)), key=lambda i: (loci[i][0], int(loci[i][1])))
    left, right = {}, {}
    for j, i in enumerate(order):
        c, s, e = loci[i][0], int(loci[i][1]), int(loci[i][2])
        pv = order[j - 1] if j > 0 and loci[order[j - 1]][0] == c else None
        nx = order[j + 1] if j + 1 < len(order) and loci[order[j + 1]][0] == c else None
        left[i] = min(a.win, s - (int(loci[pv][2]) if pv is not None else 0))
        right[i] = min(a.win, (int(loci[nx][1]) if nx is not None else size[c]) - e)
    # windows in element orientation: side 3 = past the 3' end, side 5 = before the 5' end
    with open(tmp + "/win.bed", "w") as fh:
        for i, f in enumerate(loci):
            c, s, e, st = f[0], int(f[1]), int(f[2]), f[5]
            r3, r5 = (right[i], left[i]) if st == "+" else (left[i], right[i])
            if r3 > 0:
                fh.write("%s\t%d\t%d\t%d_3\t0\t%s\n" % (c, e if st == "+" else s - r3, e + r3 if st == "+" else s, i, st))
            if r5 > 0:
                fh.write("%s\t%d\t%d\t%d_5\t0\t%s\n" % (c, s - r5 if st == "+" else e, s if st == "+" else e + r5, i, st))
    p = subprocess.run("bedtools getfasta -s -tab -nameOnly -fi %s -bed %s/win.bed" % (a.genome, tmp), shell=True, check=True,
                       stdout=subprocess.PIPE, universal_newlines=True)
    win = {}
    for l in p.stdout.splitlines():
        n, sq = l.split("\t")
        n = n.rsplit("(", 1)[0]
        i, side = n.split("_")
        sq = sq if side == "3" else sq[::-1]          # read outward: index 0 next to the locus end
        win[(int(i), side)] = sq
    codes = {key: kmer_ints(sq, a.k) for key, sq in win.items()}
    # 5' windows were reversed as strings, so their K-mers are reversed sequences; counting and walking use the same strings,
    # which is all that matters for sharing between copies
    by_group = collections.defaultdict(list)
    for key in codes:
        by_group[(group[key[0]], key[1])].append(key)
    add = collections.defaultdict(lambda: [0, 0])          # i -> [bp 5', bp 3']
    stats = []
    for (g, side), keys in sorted(by_group.items()):
        n = len(keys)
        u = [np.unique(codes[k][codes[k] >= 0]) for k in keys]
        allc = np.concatenate(u) if u else np.empty(0, dtype=np.int64)
        vals, cnt = np.unique(allc, return_counts=True)
        thr = max(a.min_n, int(np.ceil(a.min_frac * n)))
        shared = vals[cnt >= thr]
        if len(shared) < a.min_shared:                  # one or two common K-mers (rsi: 1 in r7, P18, C11) are not an element part
            shared = shared[:0]
        ext, capped = [], 0
        if len(shared):
            for key in keys:
                b = walk(codes[key], shared, a.k, a.step, a.share, a.win, a.bridge)
                if b >= a.win - a.step:                   # still element at the cap: an array, a dimer or a longer element,
                    capped += 1                           # not an end the bank missed; left as it is (continuation check)
                    continue
                if b > 0:
                    b += a_run(win[key], b, a.win)
                    b = min(b, len(win[key]))
                    add[key[0]][0 if side == "5" else 1] = b
                    ext.append(b)
        ng = sum(1 for x in ext if x >= 10)
        stats.append((g, side, n, thr, len(shared), len(ext), ng, statistics.median(ext) if ext else 0, max(ext) if ext else 0,
                      capped, "flag" if n >= a.min_n and ng >= a.flag_frac * n else ""))
    nloci = collections.Counter(group)
    with open(a.out_dir + "/ends_extension.tsv", "w") as fh:
        fh.write("group\tside\twindows\tshared_threshold\tshared_kmers\textended\textended_ge10\tmedian_bp\tmax_bp\tat_cap_not_extended\tflag\n")
        for s in stats:
            fh.write("%s\t%s'\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%s\n" % s)
    new = {}
    for i, f in enumerate(loci):
        s, e = int(f[1]), int(f[2])
        b5, b3 = add[i] if i in add else (0, 0)
        new[i] = [s - b5, e + b3] if f[5] == "+" else [s - b3, e + b5]
    # two neighbours that both grew into the gap between them meet in its middle (loci never overlap)
    for j in range(len(order) - 1):
        i, n = order[j], order[j + 1]
        if loci[i][0] == loci[n][0] and new[i][1] > new[n][0]:
            mid = (int(loci[i][2]) + int(loci[n][1])) // 2
            new[i][1] = max(int(loci[i][2]), min(new[i][1], mid))
            new[n][0] = min(int(loci[n][1]), max(new[n][0], mid))
    out = open(a.out_bed, "w"); moved = open(a.out_dir + "/loci_extended.tsv", "w")
    moved.write("chrom\tstart\tend\tstrand\tnew_start\tnew_end\tgroup\tadded_5\tadded_3\n")
    for i, f in enumerate(loci):
        s, e = int(f[1]), int(f[2])
        ns, ne = new[i]
        b5, b3 = (s - ns, ne - e) if f[5] == "+" else (ne - e, s - ns)
        if (ns, ne) != (s, e):
            moved.write("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\n" % (f[0], s, e, f[5], ns, ne, group[i], b5, b3))
        out.write("\t".join([f[0], str(ns), str(ne)] + f[3:]) + "\n")
    out.close(); moved.close()
    subprocess.run(["rm", "-rf", tmp], check=True)
    e5 = sum(1 for v in add.values() if v[0] > 0); e3 = sum(1 for v in add.values() if v[1] > 0)
    print("EXTEND loci=%d extended5=%d extended3=%d flagged=%d -> %s/ends_extension.tsv" % (len(loci), e5, e3, sum(1 for s in stats if s[10]), a.out_dir))


if __name__ == "__main__":
    main()
