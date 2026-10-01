#!/usr/bin/env python3
"""length_variants.py - are two (or more) length versions of a SINE real versions, or one element with a decayed 3' end?

The question: copies of one family end at different places on the longest consensus. Either the family has
distinct length versions (rsi r9 105 bp / r7 154 bp; MEG-RS 135 bp / MEG-RL 207 bp, Gogolevsky 2009), each
inserted as such, or the ends only scatter because the 3' end decays or is variable. Four measurements, all
from the copies, none from the consensus lengths:

  A  END MODES      histogram of the 3' end on the long consensus (5'-complete copies only). Real versions give
                    separate modes with an empty valley between them; decay gives one mode and a smooth tail.
                    valley_ratio = density at the deepest point between two modes / the lower mode's peak density;
                    close peaks (<= 40 bp) with a shallow valley are one broad mode (a version's own end scatters).
  B  LINKAGE        columns inside the shared region where a base is >= DIAG_DIFF more frequent in the short
                    mode than in the long one. Real versions are separate lineages, so those bases follow
                    the length; decay is independent of the internal sequence, so no column differs. Each copy
                    is typed by the bases it carries at those columns; linkage = P(short type | short mode) -
                    P(short type | long mode). A short mode that also holds copies of the long type (decayed
                    long copies) lowers the index but does not remove it.
  C  TSD            target-site duplication after each mode's own end (4-20 bp search of fs8, slack 45), real
                    flank pairs against shuffled ones (5' side of copy i with the 3' side of copy i+1). A real
                    insertion carries its own TSD right after its own end.
  D  RESIDUAL       when a short mode's 3' end was lost by decay, the flank past it still resembles the longer
                    consensus' extension a little; a version that really ends there does not (flank = chance).

Verdict (thresholds below, all adjustable):
  TWO_VERSIONS      >= 2 modes, valley_ratio <= VALLEY_MAX, linkage >= LINK_MIN, TSD excess >= TSD_MIN in every mode
  UNLINKED_ENDS     >= 2 modes but the internal bases do not follow the length (a variable end, one monomer)
  SINGLE_MODE       one mode: no second version; the tail is decay / a variable end
  UNRESOLVED        too few copies, or the evidence disagrees (the report lists which test failed)
It never merges anything: it is a decision input.

Usage:
  length_variants.py --genome G.fna --hits assignment_full.tsv --families MEG-RL,MEG-RS --cons bank.fa \
      --long MEG-RL --out outprefix [--threads 8]
needs samtools, bedtools, ssearch36 on PATH.
"""
import argparse, collections, concurrent.futures as cf, os, random, re, subprocess, sys

# --- thresholds ------------------------------------------------------------------------------------------
SMOOTH = 7          # box smoothing of the end histogram (bp)
MODE_MIN_FRAC = 0.05   # a mode needs this share of the 5'-complete copies ...
MODE_MIN_N = 30        # ... and at least this many
MODE_SEP = 10       # peaks closer than this are one peak
MERGE_NEAR = 40     # neighbouring peaks this close with a shallow valley are one broad mode
MODE_FLOOR = 0.2    # the mode window extends while the density is >= this x the peak
MODE_MAXHALF = 30   # ... but at most this far from the peak
MODE_PROM = 4.0     # the peak must stand this far above the median density of its surroundings (40 bp each side)
QS_MAX = 15         # a 5'-complete copy starts within this of consensus position 1
VALLEY_MAX = 0.25
DIAG_DIFF = 0.25    # a column is diagnostic when a base is this much more frequent in the short mode than in the long mode
LINK_MIN = 0.4
TSD_MIN = 10.0      # points of TSD excess over shuffled (TSD >= TSD_LEN bp)
TSD_LEN = 7
TSD_SLACK = 45
MIN_COPIES = 100
TSD_SAMPLE = 1500

COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fa(path):
    d, n = {}, None
    for line in open(path):
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = []
        elif n is not None:
            d[n].append(line.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


def parse_btop(btop, qs):
    """BTOP of a ssearch36 -m 8CB hit -> ({query column: subject base ('-' = deleted in the copy)}, last query column).
    Digits are matches; a pair is (query base, subject base); a pair with '-' as the query base is an insertion in
    the copy and uses no query column. Matches are not listed."""
    col, out, i = qs, {}, 0
    while i < len(btop):
        if btop[i].isdigit():
            j = i
            while j < len(btop) and btop[j].isdigit():
                j += 1
            col += int(btop[i:j])
            i = j
        else:
            q, s = btop[i], btop[i + 1]
            i += 2
            if q == "-":
                continue
            out[col] = s
            col += 1
    return out, col - 1


class Copy:
    """One copy aligned to the long consensus: qs..qe consensus columns, mm = {column: copy base} at non-matches,
    up = 30 bp before the element, dn = 120 bp after its aligned end, fam = assigned family."""
    def __init__(self, qs, qe, mm, up, dn, fam=""):
        self.qs, self.qe, self.mm, self.up, self.dn, self.fam = qs, qe, mm, up, dn, fam


# --- A: end modes ----------------------------------------------------------------------------------------
def _valley(sm, p, q):
    """(position, ratio) of the deepest point between peaks p < q: mean density in a 5 bp box there / the lower peak."""
    v = min(range(p, q + 1), key=lambda i: sm[i])
    box = sm[max(p, v - 2):min(q, v + 2) + 1]
    return v, (sum(box) / float(len(box))) / min(sm[p], sm[q])


def end_modes(copies, cons_len):
    """Modes of the 3' end (5'-complete copies). Neighbouring peaks <= MERGE_NEAR bp apart that no real valley separates
    (ratio > VALLEY_MAX) are one broad mode: the end of one version scatters by a few bp in the simple-repeat tail
    (rsi r9: 90-125 on the r7 consensus), and that must not look like two versions."""
    ends = [c.qe for c in copies if c.qs <= QS_MAX]
    n = len(ends)
    h = [0] * (cons_len + 2 + SMOOTH)
    for e in ends:
        h[min(e, cons_len)] += 1
    half = SMOOTH // 2
    sm = [sum(h[max(0, i - half):i + half + 1]) / float(SMOOTH) for i in range(len(h))]
    need = max(MODE_MIN_N / float(SMOOTH), MODE_MIN_FRAC * n / float(SMOOTH))
    peaks = [i for i in range(1, len(sm) - 1) if sm[i] >= need and sm[i] >= sm[i - 1] and sm[i] >= sm[i + 1]]
    peaks.sort(key=lambda i: -sm[i])
    kept = []
    for pk in peaks:
        if all(abs(pk - k) >= MODE_SEP for k in kept):
            kept.append(pk)
    kept.sort()
    # prominence: a bump in a flat tail is not a mode
    prom = []
    for pk in kept:
        around = sorted(sm[max(1, pk - 40):max(1, pk - 10)] + sm[pk + 11:pk + 41])
        floor = around[len(around) // 2] if around else 0.0
        if sm[pk] >= MODE_PROM * floor:
            prom.append(pk)
    # groups of peaks: merge near neighbours with a shallow valley, repeat until stable
    groups = [[pk] for pk in prom]
    changed = True
    while changed and len(groups) > 1:
        changed = False
        for k in range(len(groups) - 1):
            a, b = max(groups[k], key=lambda x: sm[x]), max(groups[k + 1], key=lambda x: sm[x])
            if b - a <= MERGE_NEAR and _valley(sm, a, b)[1] > VALLEY_MAX:
                groups[k] = groups[k] + groups[k + 1]
                del groups[k + 1]
                changed = True
                break
    modes = []
    for g in groups:
        pk = max(g, key=lambda x: sm[x])
        lo, hi = min(g), max(g)
        while lo > 1 and pk - lo < MODE_MAXHALF and sm[lo - 1] >= MODE_FLOOR * sm[pk]:
            lo -= 1
        while hi < len(sm) - 1 and hi - pk < MODE_MAXHALF and sm[hi + 1] >= MODE_FLOOR * sm[pk]:
            hi += 1
        modes.append({"peak": pk, "lo": lo, "hi": hi, "density": sm[pk]})
    for k in range(len(modes) - 1):
        a, b = modes[k], modes[k + 1]
        v, ratio = _valley(sm, a["peak"], b["peak"])
        a["hi"] = min(a["hi"], v)
        b["lo"] = max(b["lo"], v + 1)
        b["valley_ratio_before"] = ratio
        b["between_n"] = sum(h[a["peak"] + 1:b["peak"]])
    return modes, h, n


def assign_mode(c, modes):
    if c.qs > QS_MAX:
        return None
    for k, m in enumerate(modes):
        if m["lo"] <= c.qe <= m["hi"]:
            return k
    return None


# --- B: linkage ------------------------------------------------------------------------------------------
def column_alleles(copies, cons, col_lo, col_hi):
    """per column in [col_lo, col_hi]: Counter of the copies' bases (consensus base for matches)."""
    res = {}
    for col in range(col_lo, col_hi + 1):
        cnt, ref = collections.Counter(), cons[col - 1]
        for c in copies:
            if c.qs <= col <= c.qe:
                cnt[c.mm.get(col, ref)] += 1
        res[col] = cnt
    return res


def linkage(copies_by_mode, cons, shared_hi):
    """Short mode = the first, long mode = the last. Returns dict with the diagnostic columns and the index."""
    sh, lg = copies_by_mode[0], copies_by_mode[-1]
    a_s, a_l = column_alleles(sh, cons, 1, shared_hi), column_alleles(lg, cons, 1, shared_hi)
    diag = {}
    for col in range(1, shared_hi + 1):
        ns, nl = sum(a_s[col].values()), sum(a_l[col].values())
        if ns < 50 or nl < 50:
            continue
        # the base whose frequency is highest in the short mode relative to the long mode; a mode may be a mixture
        # (rle MEG-RS mode: 41 % G at column 25 against 6 % in the long mode), so a difference in frequency counts, not a fixed one
        bs = max("ACGT-", key=lambda b: a_s[col][b] / float(ns) - a_l[col][b] / float(nl))
        gain = a_s[col][bs] / float(ns) - a_l[col][bs] / float(nl)
        bl = a_l[col].most_common(1)[0][0]
        if gain >= DIAG_DIFF and bs != bl:      # bs == bl: the long mode's own majority base, it cannot tell the two apart
            diag[col] = (bs, bl)

    def kind(c):
        vs = vl = 0
        for col, (bs, bl) in diag.items():
            if c.qs <= col <= c.qe:
                b = c.mm.get(col, cons[col - 1])
                vs += b == bs
                vl += b == bl
        return "S" if vs > vl else "L" if vl > vs else None

    ks = collections.Counter(kind(c) for c in sh)
    kl = collections.Counter(kind(c) for c in lg)
    ps = ks["S"] / float(ks["S"] + ks["L"]) if ks["S"] + ks["L"] else None
    pl = kl["S"] / float(kl["S"] + kl["L"]) if kl["S"] + kl["L"] else None
    idx = (ps - pl) if ps is not None and pl is not None else None
    if not diag:
        idx = 0.0       # no column differs between the modes: nothing follows the length
    return {"diag_columns": diag, "short_mode_short_allele": ps, "long_mode_short_allele": pl, "index": idx,
            "n_short_informative": ks["S"] + ks["L"], "n_long_informative": kl["S"] + kl["L"]}


# --- C: TSD ----------------------------------------------------------------------------------------------
def has_tsd(up, dn, minlen=TSD_LEN, slack=TSD_SLACK):
    """fs8 rule: 5' copy ends <= 4 bp before the element, 3' copy starts <= slack bp after its end, 4-20 bp, <= 20 % mismatch."""
    for L in range(20, minlen - 1, -1):
        for d5 in range(5):
            if len(up) < L + d5:
                continue
            a = up[len(up) - L - d5:len(up) - d5]
            for d3 in range(slack + 1):
                b = dn[d3:d3 + L]
                if len(b) < L:
                    break
                if sum(x != y for x, y in zip(a, b)) <= L * 0.2:
                    return True
    return False


def tsd_rates(cs, rng):
    smp = cs if len(cs) <= TSD_SAMPLE else rng.sample(cs, TSD_SAMPLE)
    if len(smp) < 2:
        return None, None
    real = sum(has_tsd(c.up, c.dn) for c in smp) / float(len(smp))
    shuf = sum(has_tsd(smp[i].up, smp[(i + 1) % len(smp)].dn) for i in range(len(smp))) / float(len(smp))
    return 100 * real, 100 * shuf


# --- D: residual similarity past a short mode's end -------------------------------------------------------
def residual(cs, cons, end_col, rng):
    ext = cons[end_col:end_col + 40]
    if len(ext) < 30 or len(cs) < 2:
        return None

    def best(dn):
        top = 0.0
        for d in range(0, TSD_SLACK + 1):
            seg = dn[d:d + len(ext)]
            if len(seg) < len(ext):
                break
            top = max(top, sum(x == y for x, y in zip(seg, ext)) / float(len(ext)))
        return top

    smp = cs if len(cs) <= 600 else rng.sample(cs, 600)
    real = sorted(best(c.dn) for c in smp)
    shuf = sorted(best(smp[(i + 1) % len(smp)].dn) for i in range(len(smp)))
    med = lambda v: v[len(v) // 2]
    ge = lambda v: sum(x >= 0.75 for x in v) / float(len(v))
    return {"median_real": med(real), "median_shuffled": med(shuf), "frac75_real": ge(real), "frac75_shuffled": ge(shuf)}


# --- verdict ---------------------------------------------------------------------------------------------
def analyse(copies, cons, seed=7):
    rng = random.Random(seed)
    L = len(cons)
    modes, hist, n5 = end_modes(copies, L)
    rep = {"n_copies": len(copies), "n_5prime_complete": n5, "modes": modes, "hist": hist, "tests": {}, "why": []}
    if n5 < MIN_COPIES:
        rep["verdict"] = "UNRESOLVED"
        rep["why"].append("fewer than %d 5'-complete copies" % MIN_COPIES)
        return rep
    if len(modes) < 2:
        rep["verdict"] = "SINGLE_MODE"
        rep["why"].append("%d end mode(s): no second length version; copies ending away from it are decay or a variable end" % len(modes))
        if modes:
            m = modes[0]
            rep["tests"]["share_in_mode"] = sum(hist[m["lo"]:m["hi"] + 1]) / float(n5)
        return rep
    by = [[] for _ in modes]
    for c in copies:
        k = assign_mode(c, modes)
        if k is not None:
            by[k].append(c)
    rep["mode_sizes"] = [len(b) for b in by]
    ok = True
    worst = max(m.get("valley_ratio_before", 0) for m in modes[1:])
    rep["tests"]["valley_ratio_max"] = worst
    if worst > VALLEY_MAX:
        ok = False
        rep["why"].append("valley between modes not empty (ratio %.2f > %.2f)" % (worst, VALLEY_MAX))
    shared_hi = max(1, modes[0]["peak"] - 10)
    lk = linkage(by, cons, shared_hi)
    rep["tests"]["linkage"] = lk
    linked = lk["index"] is not None and lk["index"] >= LINK_MIN
    if not linked:
        rep["why"].append("internal bases do not follow the length (index %s, %d diagnostic columns)" % (
            "n/a" if lk["index"] is None else "%.2f" % lk["index"], len(lk["diag_columns"])))
    tsd = []
    for k, b in enumerate(by):
        r, s = tsd_rates(b, rng)
        tsd.append({"mode_end": modes[k]["peak"], "n": len(b), "tsd_real": r, "tsd_shuffled": s,
                    "excess": None if r is None else r - s})
    rep["tests"]["tsd"] = tsd
    tsd_ok = all(t["excess"] is not None and t["excess"] >= TSD_MIN for t in tsd)
    if not tsd_ok:
        rep["why"].append("TSD excess under %.0f points in at least one mode" % TSD_MIN)
    rep["tests"]["residual_short_mode"] = residual(by[0], cons, modes[0]["peak"], rng)
    if ok and linked and tsd_ok:
        rep["verdict"] = "TWO_VERSIONS"
    elif ok and not linked and lk["index"] is not None:
        rep["verdict"] = "UNLINKED_ENDS"
    else:
        rep["verdict"] = "UNRESOLVED"
    return rep


def format_report(rep, cons_len):
    o = ["VERDICT: %s" % rep["verdict"], "copies %d, 5'-complete %d, consensus %d bp" % (rep["n_copies"], rep["n_5prime_complete"], cons_len)]
    for m in rep["modes"]:
        o.append("  mode at %d (window %d-%d)%s" % (m["peak"], m["lo"], m["hi"], "" if "valley_ratio_before" not in m else
                 ", valley ratio to the previous mode %.3f (%d copies between the peaks)" % (m["valley_ratio_before"], m["between_n"])))
    t = rep["tests"]
    if "linkage" in t:
        lk = t["linkage"]
        o.append("  linkage: %d diagnostic columns %s; P(short allele | short mode) %s, P(short allele | long mode) %s, index %s" % (
            len(lk["diag_columns"]), sorted(lk["diag_columns"])[:12], lk["short_mode_short_allele"], lk["long_mode_short_allele"], lk["index"]))
    for x in t.get("tsd", []):
        o.append("  TSD (>=%d bp) mode %d: n=%d real %s %% shuffled %s %% excess %s" % (TSD_LEN, x["mode_end"], x["n"], x["tsd_real"], x["tsd_shuffled"], x["excess"]))
    if t.get("residual_short_mode"):
        o.append("  residual similarity past the short mode: %s" % t["residual_short_mode"])
    for w in rep["why"]:
        o.append("  note: " + w)
    return "\n".join(o)


# --- data collection -------------------------------------------------------------------------------------
def collect(genome, hits_tsv, families, cons_fa, long_name, workdir, threads, flank=100, max_copies=0):
    os.makedirs(workdir, exist_ok=True)
    cons = read_fa(cons_fa)[long_name]
    if not os.path.exists(genome + ".fai"):
        subprocess.check_call(["samtools", "faidx", genome])
    fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(genome + ".fai")}
    rows, skipped = [], 0
    for line in open(hits_tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or f[4] != "assigned" or f[1] not in families:
            continue
        m = re.match(r"(.+):(\d+)-(\d+)\(([+-])\)$", f[0])
        if not m:
            skipped += 1
            continue
        rows.append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), f[1], f[0]))
    if max_copies and len(rows) > max_copies:
        rows = random.Random(1).sample(rows, max_copies)
    with open(workdir + "/hits.bed", "w") as o:
        for i, (c, s, e, st, fam, h) in enumerate(rows):
            s0, e0 = max(0, s - 1 - flank), min(fai[c], e + flank)
            o.write("%s\t%d\t%d\tc%d\t0\t%s\n" % (c, s0, e0, i, st))   # short names: ssearch36 cuts ids at 60 characters
    subprocess.check_call(["bedtools", "getfasta", "-fi", genome, "-bed", workdir + "/hits.bed", "-s", "-name", "-fo", workdir + "/copies.fa"])
    with open(workdir + "/cons.fa", "w") as o:
        o.write(">%s\n%s\n" % (long_name, cons))
    seqs = {k.split("::")[0]: v for k, v in read_fa(workdir + "/copies.fa").items()}
    fam_of = {"c%d" % i: r[4] for i, r in enumerate(rows)}
    # chunks of 1000, searched in parallel; -z 11 because the default statistics return nothing on a large library of near-identical copies
    names = list(seqs)
    chunks = []
    for i in range(0, len(names), 1000):
        p = "%s/part_%04d.fa" % (workdir, i // 1000)
        with open(p, "w") as o:
            for n in names[i:i + 1000]:
                o.write(">%s\n%s\n" % (n, seqs[n]))
        chunks.append(p)

    def run(p):
        return subprocess.run(["ssearch36", "-m", "8CB", "-E", "1e-3", "-d", "5", "-z", "11", "-T", "1", workdir + "/cons.fa", p],
                              capture_output=True, text=True).stdout

    best = {}
    with cf.ThreadPoolExecutor(threads) as ex:
        for out in ex.map(run, chunks):
            for line in out.splitlines():
                f = line.split("\t")
                if line.startswith("#") or len(f) < 13:
                    continue
                b = float(f[11])
                if f[1] not in best or b > float(best[f[1]][11]):
                    best[f[1]] = f
    copies = []
    for name, f in best.items():
        s = seqs[name]
        qs, qe, ss, se = int(f[6]), int(f[7]), int(f[8]), int(f[9])
        if ss > se:
            s = rc(s)
            n = len(s)
            ss, se = n - ss + 1, n - se + 1
        mm, last = parse_btop(f[12], qs)
        up_end = max(0, ss - 1 - (qs - 1))
        copies.append(Copy(qs, last, mm, s[max(0, up_end - 30):up_end], s[se:se + 120], fam_of[name]))
    return copies, cons, skipped, len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--hits", required=True, help="assignment_full.tsv")
    ap.add_argument("--families", required=True, help="comma list pooled for the test (the long and the short consensus names)")
    ap.add_argument("--cons", required=True, help="consensus fasta")
    ap.add_argument("--long", required=True, help="name of the longest consensus")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-copies", type=int, default=0)
    a = ap.parse_args()
    copies, cons, skipped, nrows = collect(a.genome, a.hits, set(a.families.split(",")), a.cons, a.long, a.out + ".work", a.threads, max_copies=a.max_copies)
    rep = analyse(copies, cons)
    txt = format_report(rep, len(cons)) + "\n(%d hits used, %d with mixed strand skipped)\n" % (nrows, skipped)
    open(a.out + ".report.txt", "w").write(txt)
    with open(a.out + ".ends.tsv", "w") as o:
        o.write("end\tcopies_5prime_complete\n")
        for i, v in enumerate(rep["hist"]):
            if v:
                o.write("%d\t%d\n" % (i, v))
    with open(a.out + ".copies.tsv", "w") as o:
        o.write("fam\tqs\tqe\n")
        for c in copies:
            o.write("%s\t%d\t%d\n" % (c.fam, c.qs, c.qe))
    print(txt)


if __name__ == "__main__":
    main()
