#!/usr/bin/env python3
"""satellite_stage.py: the satellite screen of SINEderella, run inside step 1 after the searches and before extraction / SubFam.

Design: docs/SATELLITES.md (sections 4 and 5). Per consensus:
  kind A  SINE-derived satellites (monomer = a part of the SINE, or the SINE plus a bit): tools/satellite_trf_verify on the raw hits of
          sear (searches/gen-<q>.rawhits.tsv: merged hits at the homology cut, any length; falls back to the length-filtered gen-<q>.bed),
          i.e. a window around every hit, TRF on the windows, unit aligned to the consensus. Every verified locus is a satellite locus.
  kind B  arrays whose unit is longer than the SINE (rsi MEG-RS): tools/satellite_screen regular-spacing runs with a chance null on the
          full-length hits (gen-<q>.bed); the family is SAT_B when the excess over chance is >= 20 % of its hits.
Writes OUT/indication.tsv (per consensus), OUT/loci.bed (kind, consensus, locus, monomers, SINE part), OUT/units.fa, OUT/<q>.kindA.loci.tsv.
Exclusion (default on; --no-exclude writes the tables only): hits of gen-<q>.bed that overlap a kind-A locus, or a kind-B run that counts
as a satellite locus, are removed. --exclude-b decides which kind-B runs count: `verified` (default) = runs whose units are near-identical
(satellite_kindB_verify: units cut from hit to hit, neighbours and next-but-one aligned, median identity >= 85 %; rsi: 18 of 21 MEG-RS runs,
the 26-unit 2.6 kb array at NC_142508.1, and 1 850 of 1 862 five-to-nine-unit runs of the r-families are NOT arrays); `flagged` = runs of
SAT_B consensuses only;
`long` = those plus, in any consensus, runs with at least --long-min-units (10) units and at least the chance-calibrated minimum (the smallest
length chance cannot explain: rle MEG-RS 5, tbr VES 50). Kind-B runs are geometric only (no sequence check of the units yet), and in rsi the
calibrated minimum was 5 for most r-families, which would have counted hundreds of 5-unit clusters of ordinary copies: hence the 10-unit floor
and `flagged` as default. `all` = every run; `none`. They are removed from gen-<q>.bed (the original is kept as gen-<q>.bed.before_satellites and the
removed hits in OUT/excluded_hits.bed), so that extraction, SubFam (the peel input) and assignment never see them. Nothing is deleted.

Usage: satellite_stage.py --searches DIR --genome G.fa --cons CONSENSUSES.fa --out OUT [--threads 16] [--no-exclude]
       [--exclude-b verified|flagged|long|all|none] [--max-verify-runs 2000] [--only NAME,...]
"""
import argparse
import collections
import glob
import os
import re
import shutil
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import array_order as ao  # noqa: E402
import satellite_screen as ss  # noqa: E402
import satellite_kindB_verify as kb  # noqa: E402

FLAG_PCT = 20.0


def read_fa(path):
    d, n = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = []
        elif n is not None:
            d[n].append(line)
    return {k: "".join(v) for k, v in d.items()}


def query_name(bed):
    """searches/gen-consensuses.clean.part_<NAME>.bed -> NAME"""
    b = os.path.basename(bed)
    b = re.sub(r"\.(bed|rawhits\.tsv)$", "", b)
    return b.split(".part_", 1)[1] if ".part_" in b else re.sub(r"^gen-", "", b)


def overlaps(lst, s, e):
    """lst sorted by start: any interval overlapping [s, e)"""
    import bisect
    i = bisect.bisect_left(lst, (s, -1))
    j = i - 1
    while j >= 0 and lst[j][1] > s:
        if lst[j][0] < e:
            return True
        j -= 1
    while i < len(lst) and lst[i][0] < e:
        if lst[i][1] > s:
            return True
        i += 1
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--searches", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--cons", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--no-exclude", action="store_true")
    ap.add_argument("--exclude-b", choices=["verified", "flagged", "long", "all", "none"], default="verified")
    ap.add_argument("--long-min-units", type=int, default=10, help="with --exclude-b long: a run counts on its own only with at least this many units and more than the chance-calibrated minimum")
    ap.add_argument("--max-verify-runs", type=int, default=2000, help="per consensus, the unit check (ssearch36 per run) is done on at most this many regularly spaced runs, the longest first; the others stay 'untested' and are not excluded")
    ap.add_argument("--only", default="", help="comma list of consensus names: screen only these (SINEderella --add: the new consensuses) and merge their rows into the existing tables of OUT")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cons = read_fa(a.cons)
    only = set(x for x in a.only.split(",") if x)
    beds = sorted(glob.glob(os.path.join(a.searches, "gen-*.bed")))
    beds = [b for b in beds if not b.endswith(".before_satellites")]
    if only:
        beds = [b for b in beds if query_name(b) in only or query_name(b).replace("|", "_") in {x.replace("|", "_") for x in only}]
    if not beds:
        print("satellite_stage: no gen-*.bed in %s%s, nothing to do" % (a.searches, " for %s" % ",".join(sorted(only)) if only else ""))
        return 0
    verify = os.path.join(HERE, "satellite_trf_verify.py")
    rows, loci, units, excluded, allA, regions = [], [], [], [], [], []
    for bed in beds:
        q = query_name(bed)
        if q not in cons:
            # the bed may carry the file-part name while the bank carries the id; try both
            cand = [k for k in cons if k == q or k.replace("|", "_") == q]
            if not cand:
                print("satellite_stage: %s: consensus %s not in %s, skipped" % (os.path.basename(bed), q, a.cons))
                continue
            q = cand[0]
        raw = re.sub(r"\.bed$", ".rawhits.tsv", bed)
        hits_src = raw if os.path.exists(raw) and os.path.getsize(raw) else bed
        cfa = os.path.join(a.out, q + ".cons.fa")
        open(cfa, "w").write(">%s\n%s\n" % (q, cons[q]))
        pref = os.path.join(a.out, q + ".kindA")
        subprocess.run([sys.executable, verify, hits_src, "--genome", a.genome, "--cons", cfa, "--out", pref, "--threads", str(a.threads)],
                       check=False, stdout=open(pref + ".log", "w"), stderr=subprocess.STDOUT)
        A = []
        if os.path.exists(pref + ".loci.tsv"):
            for l in list(open(pref + ".loci.tsv"))[1:]:
                f = l.rstrip("\n").split("\t")
                # contig start end period copies cons_start cons_end score identity monomer_cov cons_cov
                per_, aln_, pid_, cs_, ce_ = int(f[3]), int(f[7]), float(f[8]), int(f[9]), int(f[10])
                A.append((f[0], int(f[1]), int(f[2]), per_, float(f[4]), cs_, ce_, aln_ * pid_, pid_, min(1.0, aln_ / float(per_)),
                          (ce_ - cs_ + 1) / float(len(cons[q]))))
        allA.extend((q,) + x for x in A)
        if os.path.exists(pref + ".units.fa"):
            units.append(open(pref + ".units.fa").read())
        # kind B on the full-length hits
        full = ss.read_bed(bed)
        brows, bsum = ss.screen(full)
        flagB = bsum["excess_B"] >= FLAG_PCT
        B = [(r[1], r[2], r[3], r[4], r[7]) for r in brows if r[0] == "B"]      # contig start end hits median_gap
        # The wide tier (>= 10 copies, gaps <= 30 kb) replaces the narrow runs it contains. When it joins a real array with
        # dispersed copies beside it, the median unit identity of the joined run can fall under the cut and the whole array
        # stays in (rsi MEG-RS, mini genome 2026-10-05: a 2 168 bp array verified at 89 % on its own became one 144-hit run with
        # the copies 2-35 kb beyond its end, 82 %, COPIES). So the narrow runs (6 kb rule) are verified as well, on their own, and
        # a verified narrow run is excluded even when the wide run that contains it is not.
        nar = ao.regular_runs([(c, s) for c, s, e, _st in full])
        ngroups = collections.defaultdict(list)
        for i, rid in nar.items():
            ngroups[rid].append(full[i])
        have = {(x[0], x[1], x[2]) for x in B}
        Bn = []
        for v in ngroups.values():
            v.sort(key=lambda h: h[1])
            gaps = [v[t + 1][1] - v[t][1] for t in range(len(v) - 1)]
            rec = (v[0][0], v[0][1], v[-1][2], len(v), int(statistics.median(gaps)) if gaps else 0)
            if (rec[0], rec[1], rec[2]) not in have:
                Bn.append(rec)
        # the unit check of every run (sequence, not geometry): verified arrays are satellite loci whatever the family share
        verdicts = {}
        if B and a.exclude_b != "none":
            # one ssearch36 call per run: a hit-dense family (tbr VES, one copy per 3 kb) has tens of thousands of chance runs of
            # 5-9 copies, so only the --max-verify-runs longest runs are checked; the rest keep the geometric verdict "untested"
            # and are never excluded (2026-10-05). Real arrays are long and come first.
            ordered = sorted(B + Bn, key=lambda x: -x[3])
            todo, rest = ordered[:a.max_verify_runs], ordered[a.max_verify_runs:]
            if rest:
                print("satellite_stage: %s: %d regularly spaced runs, unit check on the %d longest (--max-verify-runs); %d runs of <= %d hits untested"
                      % (q, len(B), len(todo), len(rest), rest[0][3]), flush=True)
            runs = [(q, c, s, e, gap, n) for c, s, e, n, gap in todo]
            st = {q: collections.defaultdict(list)}
            for c, s, e, _strand in full:
                st[q][c].append((s, e))
            for c in st[q]:
                st[q][c].sort()
            res = kb.verify_runs(runs, st, a.genome, a.out, threads=a.threads)
            verdicts = {(runs[i][1], runs[i][2]): res[i] for i in res}
            for c, s, e, n, gap in rest:
                verdicts[(c, s)] = (0, None, None, "untested")
        arrB = [x for x in B if verdicts.get((x[0], x[1]), (0, None, None, "-"))[3] == "ARRAY"]
        # narrow runs verified inside a wide run that failed: excluded too, listed in loci.bed with sine_part "narrow"
        arrN = [x for x in Bn if verdicts.get((x[0], x[1]), (0, None, None, "-"))[3] == "ARRAY"
                and not any(w[0] == x[0] and w[1] <= x[1] and x[2] <= w[2] for w in arrB)]
        arrB = arrB + arrN
        lmin = bsum.get("long_min")
        longB = [x for x in B if lmin is not None and x[3] >= max(lmin, a.long_min_units)]
        mono = int(sum(x[4] for x in A))
        rows.append((q, len(cons[q]), len(full), len(A), mono, max((x[4] for x in A), default=0), len(B), bsum["pct_B"], bsum["excess_B"],
                     lmin if lmin is not None else "-", len(longB), len(arrB), "SAT_A" if A else "-", "SAT_B" if flagB else "-"))   # flag_A re-set below
        for c, s, e, n, gap in B:
            v = verdicts.get((c, s), (0, None, None, "-"))
            loci.append(("B", q, c, s, e, gap, n, "-", v[3] + ("" if v[1] is None else " id%.0f" % v[1])))
        for c, s, e, n, gap in arrN:
            v = verdicts.get((c, s), (0, None, None, "-"))
            loci.append(("B", q, c, s, e, gap, n, "narrow", v[3] + ("" if v[1] is None else " id%.0f" % v[1])))
        print("satellite_stage: %s: %d full hits; kind A %d loci (%d monomers); kind B %d runs (%d verified arrays, %d long), excess %.1f %% %s%s" % (
            q, len(full), len(A), mono, len(B), len(arrB), len(longB), bsum["excess_B"],
            "SAT_A " if A else "", "SAT_B" if flagB else ""), flush=True)
        # exclusion
        if a.no_exclude:
            continue
        ex = collections.defaultdict(list)
        for rec in A:
            ex[rec[0]].append((rec[1], rec[2]))
        if a.exclude_b == "all" or (a.exclude_b in ("flagged", "long") and flagB):
            for c, s, e, n, gap in B:
                ex[c].append((s, e))
        elif a.exclude_b == "long":
            for c, s, e, n, gap in longB:
                ex[c].append((s, e))
        elif a.exclude_b == "verified":
            for c, s, e, n, gap in arrB:
                ex[c].append((s, e))
        if not ex:
            continue
        for c in ex:
            ex[c].sort()
            # the same spans go to exclude_regions.bed: step 1 removes every MERGED locus inside them, whatever consensus found
            # it. Filtering only this consensus' own hits let array units re-enter through another consensus' hit at the same
            # place (rsi MEG-RS 2026-10-05: 6 of the 129 copies left after the screen lay inside verified array spans).
            kindA = {(s, e) for s, e in ((rec[1], rec[2]) for rec in A if rec[0] == c)}
            regions += ["%s\t%d\t%d\t%s\t%s\n" % (c, s, e, "A" if (s, e) in kindA else "B", q) for s, e in ex[c]]
        keep, drop = [], []
        for line in open(bed):
            f = line.rstrip("\n").split("\t")
            try:
                c, s, e = f[0], int(f[1]), int(f[2])
            except (ValueError, IndexError):
                keep.append(line)
                continue
            (drop if overlaps(ex.get(c, []), s, e) else keep).append(line)
        if drop:
            if not os.path.exists(bed + ".before_satellites"):
                shutil.copy2(bed, bed + ".before_satellites")
            open(bed, "w").writelines(keep)
            excluded += [q + "\t" + l for l in drop]
            print("satellite_stage: %s: %d of %d full hits removed from %s (kept in .before_satellites)" % (q, len(drop), len(keep) + len(drop), os.path.basename(bed)), flush=True)
    # one locus, one consensus: the same array is found through every consensus that shares its head (rsi: a 140 bp monomer
    # matched all 17 r-families). Derivation rule (tested on the rsi shared loci, docs/SATELLITES.md 5h): group the records of one
    # locus; among the consensuses within 3 identity points of the best and covered by the monomer over >= 90 % of the alignment
    # (if none, all of them), take the consensus of which the monomer covers the largest share (a whole r2 beats the r2 part of a
    # composite); ties -> the shorter consensus. The others are listed as also_matches.
    groups = []      # [contig, start, end, [records]]
    for rec in sorted(allA, key=lambda x: (x[1], x[2])):
        c, s, e = rec[1], rec[2], rec[3]
        for g in groups:
            if g[0] == c and min(g[2], e) - max(g[1], s) > 0.5 * min(g[2] - g[1], e - s):
                g[3].append(rec)
                g[1], g[2] = min(g[1], s), max(g[2], e)
                break
        else:
            groups.append([c, s, e, [rec]])
    attributed, also = [], collections.defaultdict(set)
    for c, gs, ge, recs in groups:
        # record = (q, contig, start, end, period, copies, cons_start, cons_end, score, identity, monomer_cov, cons_cov)
        best_id = max(r[9] for r in recs)
        cand = [r for r in recs if r[9] >= best_id - 3.0 and r[10] >= 0.9] or [r for r in recs if r[9] >= best_id - 3.0]
        pick = max(cand, key=lambda r: (r[11], -len(cons[r[0]])))
        q, _, s, e, per, cop, cs, ce = pick[:8]
        attributed.append((q, c, s, e, per, cop, cs, ce))
        for r in recs:
            if r[0] != q:
                also[(c, s)].add(r[0])
    per_q = collections.Counter(x[0] for x in attributed)
    mono_q = collections.defaultdict(float)
    big_q = collections.defaultdict(float)
    for q, c, s, e, per, cop, cs, ce in attributed:
        mono_q[q] += cop
        big_q[q] = max(big_q[q], cop)
        loci.append(("A", q, c, s, e, per, cop, "%d-%d" % (cs, ce), ",".join(sorted(also.get((c, s), ()))) or "-"))
    rows = [(r[0], r[1], r[2], per_q.get(r[0], 0), int(mono_q.get(r[0], 0)), big_q.get(r[0], 0)) + r[6:12] + ("SAT_A" if per_q.get(r[0], 0) else "-", r[13]) for r in rows]
    HDR = "consensus\tcons_len\tfull_hits\tkindA_loci\tkindA_monomers\tkindA_largest\tkindB_runs\tkindB_pct\tkindB_excess_pct\tkindB_long_min\tkindB_long_runs\tkindB_verified_arrays\tflag_A\tflag_B\n"
    LHDR = "#kind\tconsensus\tcontig\tstart\tend\tperiod_or_unit\tmonomers_or_hits\tsine_part\talso_matches_or_unit_check\n"
    done = {r[0] for r in rows}

    def kept_rows(path, header, col):
        """with --only: the rows of the existing table for the consensuses not screened now (same header only; another format is
        set aside as .prev and reported)"""
        if not only or not os.path.exists(path):
            return []
        lines = open(path).read().splitlines(True)
        if not lines or lines[0] != header:
            os.replace(path, path + ".prev")
            print("satellite_stage: %s had another format; kept as .prev, rewritten with the screened consensuses only" % os.path.basename(path))
            return []
        return [l for l in lines[1:] if l.split("\t")[col] not in done]

    old_ind = kept_rows(os.path.join(a.out, "indication.tsv"), HDR, 0)
    old_loci = kept_rows(os.path.join(a.out, "loci.bed"), LHDR, 1)
    with open(os.path.join(a.out, "indication.tsv"), "w") as o:
        o.write(HDR)
        o.writelines(old_ind)
        for r in rows:
            o.write("%s\t%d\t%d\t%d\t%d\t%.0f\t%d\t%.1f\t%.1f\t%s\t%d\t%d\t%s\t%s\n" % r)
    with open(os.path.join(a.out, "loci.bed"), "w") as o:
        o.write(LHDR)
        o.writelines(old_loci)
        for l in sorted(loci, key=lambda x: (x[2], x[3])):
            l = l if len(l) == 9 else l + ("-",)
            o.write("\t".join(str(x) for x in l) + "\n")
    if only:      # units of every consensus screened so far: rebuilt from the per-consensus files
        units = [open(f).read() for f in sorted(glob.glob(os.path.join(a.out, "*.kindA.units.fa")))]
    open(os.path.join(a.out, "units.fa"), "w").write("".join(units))
    exf = os.path.join(a.out, "excluded_hits.bed")
    old_ex = [l for l in open(exf)] if (only and os.path.exists(exf)) else []
    old_ex = [l for l in old_ex if l.split("\t")[0] not in done]
    with open(exf, "w") as o:
        o.writelines(old_ex)
        o.writelines(excluded)
    rgf = os.path.join(a.out, "exclude_regions.bed")       # contig start end kind(A|B) consensus; read by step 1 after the merge
    old_rg = [l for l in open(rgf)] if (only and os.path.exists(rgf)) else []
    old_rg = [l for l in old_rg if l.rstrip("\n").split("\t")[4] not in done]
    with open(rgf, "w") as o:
        o.writelines(sorted(old_rg + regions, key=lambda l: (l.split("\t")[0], int(l.split("\t")[1]))))
    na = sum(r[3] for r in rows)
    print("satellite_stage: %d consensuses, %d kind-A loci (%d distinct; %d found through more than one consensus), %d consensuses SAT_B, %d hits excluded -> %s" % (
        len(rows), na, len(attributed), len(also), sum(1 for r in rows if r[13] == "SAT_B"), len(excluded), a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
