#!/usr/bin/env python3
"""composite_scan.py RUN_DIR OUT_DIR [--sample N] [--flank F] [--threads T] [--loci IDS.txt]

Is each assigned copy a single element, or part of a COMPOSITE - two SINE copies joined head to tail
(a tandem or a dimer), which SINEderella treats as one copy or splits into two loci?

Found by hand on rsi (2026-09-28, Tal rsi/REFINEMENT.md): r5 subgroup B is the left half of r5-r5
pairs, r6's left-extended top100 copies have an r5 head in front, r10 is the left half of an
r10 + group B element, and r8 rand100 held four head-to-tail pairs. None of these is caught by the
[array] mark (>= 3 copies within 50 kb). This scans every copy.

Per assigned locus (Status "assigned" in results/assignment_full.tsv; --sample N per subfamily,
seeded): the locus +- F bp (default 200), in the locus orientation, is searched with every
consensus of the run (ssearch36, both strands, -z 11 because the library is all homologs). Hits
>= MIN_ALN bp are kept greedily by bitscore, dropping any that overlaps a kept one by > MAX_OVL bp.
The MAIN hit is the best one overlapping the locus by >= 30 bp. Another kept hit is a PARTNER when
the gap between them is <= ADJ bp (overlap up to MAX_OVL allowed): upstream partner -> the copy is
a RIGHT half, downstream partner -> a LEFT half. A partner inside the locus means the locus itself
spans the junction. Hits 30-200 bp away are counted separately as 'near'.

Baseline: with D assigned copies per bp of genome, a random neighbour starts within ADJ bp on a
given side with probability ~ D * (ADJ + element length). Printed next to the observed rates.

Outputs: OUT_DIR/composite_loci.tsv (one row per locus), OUT_DIR/composite_summary.tsv, and the
summary on stdout.
"""
import collections
import os
import random
import re
import subprocess
import sys

ADJ = 30        # max gap (bp) between two copies to call them one composite
MAX_OVL = 20    # allowed overlap between kept hits / partner and main
MIN_ALN = 40    # shortest hit kept
SEED = 42


def arg(name, default, cast=str):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def read_fai(path):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(path)}


def cons_lengths(path):
    L, name = {}, None
    for l in open(path):
        l = l.strip()
        if l.startswith(">"):
            name = l[1:].split()[0]; L[name] = 0
        elif name:
            L[name] += len(l)
    return L


LOC = re.compile(r"^(\S+?):(\d+)-(\d+)\(([^)]*)\)")


def load_loci(run):
    loci = []
    for l in open(os.path.join(run, "results", "assignment_full.tsv")):
        f = l.rstrip("\n").split("\t")
        if len(f) < 5 or f[4] != "assigned":
            continue
        m = LOC.match(f[0])
        if not m:
            continue
        st = m.group(4) if m.group(4) in ("+", "-") else "+"
        loci.append((f[0], f[1], m.group(1), int(m.group(2)), int(m.group(3)), st))
    return loci


def main():
    run, out = sys.argv[1], sys.argv[2]
    F = arg("--flank", 200, int)
    T = arg("--threads", 32, int)
    N = arg("--sample", 0, int)
    os.makedirs(out, exist_ok=True)
    genome = os.path.join(run, "genome.clean.fa")
    cons = os.path.join(run, "consensuses.clean.fa")
    fai = read_fai(genome + ".fai")
    clen = cons_lengths(cons)
    loci = load_loci(run)
    n_all = collections.Counter(x[1] for x in loci)
    if "--loci" in sys.argv:
        want = {l.strip().replace("_", "@U@") for l in open(arg("--loci", "")) if l.strip()}
        want |= {w.replace("@U@", "_") for w in want}
        loci = [x for x in loci if x[0] in want or x[0].replace("@U@", "_") in want]
    elif N:
        by = collections.defaultdict(list)
        for x in loci:
            by[x[1]].append(x)
        rng = random.Random(SEED)
        loci = [x for sf in sorted(by) for x in (by[sf] if len(by[sf]) <= N else rng.sample(by[sf], N))]
    bed = os.path.join(out, "windows.bed")
    with open(bed, "w") as fh:
        for i, (lid, sf, c, a, b, st) in enumerate(loci):
            fh.write("%s\t%d\t%d\tw%d\t0\t%s\n" % (c, max(0, a - 1 - F), min(fai.get(c, b + F), b + F), i, st))
    wfa = os.path.join(out, "windows.fa")
    subprocess.run("bedtools getfasta -fi %s -bed %s -s -nameOnly > %s" % (genome, bed, wfa), shell=True, check=True)
    hits_path = os.path.join(out, "hits.m8")
    subprocess.run("ssearch36 -m 8 -E 1e-5 -z 11 -T %d %s %s > %s 2> /dev/null" % (T, cons, wfa, hits_path),
                   shell=True, check=True)

    # window geometry: offset of the locus inside the window (window is in locus orientation)
    geo = {}
    for i, (lid, sf, c, a, b, st) in enumerate(loci):
        ws, we = max(0, a - 1 - F), min(fai.get(c, b + F), b + F)
        up = (a - 1) - ws if st == "+" else we - b
        geo["w%d" % i] = (up, up + (b - a + 1))          # core [lo, hi) in window coords
    H = collections.defaultdict(list)
    for l in open(hits_path):
        f = l.split("\t")
        if len(f) < 12:
            continue
        w = f[1].split("(")[0]
        q1, q2, s1, s2, bits = int(f[6]), int(f[7]), int(f[8]), int(f[9]), float(f[11])
        rev = (q1 > q2) != (s1 > s2)
        lo, hi = min(s1, s2) - 1, max(s1, s2)
        if hi - lo < MIN_ALN:
            continue
        H[w].append((bits, lo, hi, f[0], min(q1, q2), max(q1, q2), "-" if rev else "+"))

    rows, S = [], collections.defaultdict(collections.Counter)
    cut = collections.defaultdict(collections.Counter)       # where a LEFT half's main hit ends (cons pos)
    partners = collections.defaultdict(collections.Counter)
    nearp = collections.defaultdict(collections.Counter)      # nearest hit 30-200 bp away: side, family
    gaps = collections.defaultdict(list)                      # adjacent-partner gaps per subfamily
    for i, (lid, sf, c, a, b, st) in enumerate(loci):
        w = "w%d" % i
        clo, chi = geo[w]
        kept = []
        for h in sorted(H.get(w, []), reverse=True):
            if all(min(h[2], k[2]) - max(h[1], k[1]) <= MAX_OVL for k in kept):
                kept.append(h)
        main = next((h for h in kept if min(h[2], chi) - max(h[1], clo) >= 30), None)
        S[sf]["n"] += 1
        if not main:
            S[sf]["no_main_hit"] += 1
            rows.append((lid, sf, "no_hit", "", "", "", "", ""))
            continue
        up = dn = None
        near = 0
        near_best = None
        for h in kept:
            if h is main:
                continue
            if h[2] <= main[1] + MAX_OVL:            # upstream of main
                gap = main[1] - h[2]
                if gap <= ADJ and (up is None or h[0] > up[0]):
                    up = h + (gap,)
                elif ADJ < gap <= 200:
                    near += 1
                    if near_best is None or gap < near_best[0]:
                        near_best = (gap, "up", h[3])
            elif h[1] >= main[2] - MAX_OVL:          # downstream
                gap = h[1] - main[2]
                if gap <= ADJ and (dn is None or h[0] > dn[0]):
                    dn = h + (gap,)
                elif ADJ < gap <= 200:
                    near += 1
                    if near_best is None or gap < near_best[0]:
                        near_best = (gap, "down", h[3])
        if near_best:
            nearp[sf]["%s %s" % (near_best[1], near_best[2])] += 1
        cls = {(False, False): "single", (True, False): "right_half", (False, True): "left_half",
               (True, True): "middle"}[(up is not None, dn is not None)]
        S[sf][cls] += 1
        if near:
            S[sf]["near_30_200"] += 1
        inside = any(p is not None and min(p[2], chi) - max(p[1], clo) >= 30 for p in (up, dn))
        if inside:
            S[sf]["partner_inside_locus"] += 1
        if dn is not None:
            cut[sf][(main[5] // 10) * 10] += 1
        for p, side in ((up, "up"), (dn, "down")):
            if p is not None:
                gaps[sf].append(p[7])
                partners[sf]["%s %s %s" % (side, p[3], "same" if p[6] == main[6] else "opposite")] += 1

        def fmt(p):
            return "" if p is None else "%s:%d-%d:%s:gap%d" % (p[3], p[4], p[5], "same" if p[6] == main[6] else "opp", p[7])
        rows.append((lid, sf, cls, "%s:%d-%d" % (main[3], main[4], main[5]), fmt(up), fmt(dn), "inside" if inside else "",
                     "" if not near_best else "%s:%s:gap%d" % (near_best[1], near_best[2], near_best[0])))

    # CHAINS: every kept hit in the window, in order - the full layout around each copy.
    # Unit label: family, then h = starts at its consensus 5' end (<= 15), e = reaches its 3' end
    # (>= length - 15); link: '+' adjacent (<= ADJ bp), '~' spacer 30-200 bp, '..' further.
    # A 'mid' unit (no h) right after another is a piecewise match of one element, not a new copy.
    chains = collections.defaultdict(collections.Counter)
    chain_ex = {}
    for i, (lid, sf, c, a, b, st) in enumerate(loci):
        w = "w%d" % i
        kept = []
        for h in sorted(H.get(w, []), reverse=True):
            if all(min(h[2], k[2]) - max(h[1], k[1]) <= MAX_OVL for k in kept):
                kept.append(h)
        kept.sort(key=lambda h: h[1])
        if not kept:
            continue
        parts = []
        for j, h in enumerate(kept):
            if j:
                g = h[1] - kept[j - 1][2]
                parts.append("+" if g <= ADJ else "~" if g <= 200 else "..")
            L = clen.get(h[3], 10 ** 6)
            tag = ("h" if h[4] <= 15 else "") + ("e" if h[5] >= L - 15 else "")
            name = re.sub(r"_\d+seqs$", "", h[3])        # r10_19seqs -> r10; r10_groupB stays itself
            parts.append("%s%s%s" % (name, "[" + tag + "]" if tag else "[mid]", "" if h[6] == "+" else "(rc)"))
        key = " ".join(parts)
        chains[sf][key] += 1
        chain_ex.setdefault((sf, key), lid)
    with open(os.path.join(out, "composite_chains.tsv"), "w") as fh:
        fh.write("subfamily\tcount\tfraction\tchain\texample_locus\n")
        for sf in sorted(chains, key=lambda k: -S[k]["n"]):
            tot = sum(chains[sf].values())
            for key, v in chains[sf].most_common(15):
                fh.write("%s\t%d\t%.3f\t%s\t%s\n" % (sf, v, v / float(tot), key, chain_ex[(sf, key)]))

    with open(os.path.join(out, "composite_loci.tsv"), "w") as fh:
        fh.write("locus\tsubfamily\tclass\tmain_hit\tupstream_partner\tdownstream_partner\tpartner_in_locus\tnearest_30_200\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    gsize = sum(fai.values())
    dens = len(load_loci(run)) / float(gsize)
    hdr = ["subfamily", "assigned", "scanned", "single", "left_half", "right_half", "middle", "partner_in_locus",
           "near_30_200", "baseline_one_side", "main_cut_mode(left_halves)", "top_partners"]
    lines = ["\t".join(hdr)]
    print("%-12s %7s %6s | %6s %6s %6s %6s | %6s %6s | %5s | %-12s %s" % (
        "subfamily", "assigned", "scan", "single", "left", "right", "middle", "inLoc", "near", "base", "cut(mode)", "top partners"))
    for sf in sorted(S, key=lambda k: -S[k]["n"]):
        s = S[sf]; n = s["n"] or 1
        pct = lambda k: 100.0 * s[k] / n
        base = 100.0 * dens * (ADJ + clen.get(sf, 250))
        cm = cut[sf].most_common(1)
        cmode = "%d-%d (%d)" % (cm[0][0], cm[0][0] + 9, cm[0][1]) if cm else "-"
        tp = ", ".join("%s %d" % (k, v) for k, v in partners[sf].most_common(3))
        lines.append("\t".join([sf, str(n_all[sf]), str(s["n"])] + ["%.1f" % pct(k) for k in
                     ("single", "left_half", "right_half", "middle", "partner_inside_locus", "near_30_200")]
                     + ["%.2f" % base, cmode, tp]))
        print("%-12s %7d %6d | %5.1f%% %5.1f%% %5.1f%% %5.1f%% | %5.1f%% %5.1f%% | %4.2f%% | %-12s %s" % (
            sf, n_all[sf], s["n"], pct("single"), pct("left_half"), pct("right_half"), pct("middle"),
            pct("partner_inside_locus"), pct("near_30_200"), base, cmode, tp))
    open(os.path.join(out, "composite_summary.tsv"), "w").write("\n".join(lines) + "\n")
    print("\nnearest other copy 30-200 bp away (top 3) and adjacent-partner gap quartiles:")
    for sf in sorted(S, key=lambda k: -S[k]["n"]):
        g = sorted(gaps[sf])
        q = "gap p25/med/p75 %d/%d/%d" % (g[len(g) // 4], g[len(g) // 2], g[3 * len(g) // 4]) if g else ""
        print("  %-12s %-26s near: %s" % (sf, q, ", ".join("%s %d" % kv for kv in nearp[sf].most_common(3))))


if __name__ == "__main__":
    main()
