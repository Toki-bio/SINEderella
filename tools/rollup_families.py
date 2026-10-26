#!/usr/bin/env python3
"""Roll a finished SINEderella run up to families: write <run>/family_view/, a run directory in which every subfamily label is replaced
by its family, so the existing steps (3 postprocess, 4 plots, 6 report, 7 boundaries, 8a plates) can be run on it unchanged.

    rollup_families.py --run RUN_ROOT --families families.clean.tsv --family-cons family_consensuses.fa [--out DIR] [--votes vote_summary.tsv[.gz]]

(families.clean.tsv and family_consensuses.fa come from tools/families.py.)

What changes in family_view/
    consensuses.clean.fa                          the family (parent) consensuses
    step2/step2_output/assignment_full.tsv        Subfamily column = family; see 'status' below
    step2/step2_output/all_votes.tsv              best voted consensus -> its family
    genome.clean_step1/searches/all_hits.labeled.bed      column 4 = family
    genome.clean_step1/searches/regions.by_subfam.bed     regions rebuilt from the relabelled hits (a CONFLICT is now two FAMILIES)
    genome.clean_step1/extracted.fasta            copy of the original (same copies)
    family_view.tsv, family_view.summary.tsv      per-copy and per-family tables

Status of a copy in the family view
    * subfamily-firm in the run (status = assigned)  -> assigned to its family (the vote and score as they were);
    * otherwise, WHEN a vote table is given (--votes; step 2 writes vote_summary.tsv.gz): the family votes of its ten cycles are summed;
      ten of ten in one family and a score at or above the lowest threshold of that family's subfamilies  -> assigned to the family with
      `layer2=unresolved` in the runner tag (family-firm, subfamily between its members);
    * all other copies keep their status; their label becomes the family of their best vote.
    A runner-up that belongs to the copy's own family is dropped from the runner tag (a LEAK inside a family is expected; across families
    it stays the alarm).
Without --votes only the first rule applies. Approximation: the runner of a family call is the best runner of the winning consensus rows.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import os
import re
import shutil
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_family_map(path):
    m = {}
    for i, line in enumerate(open(path)):
        f = line.rstrip("\r\n").split("\t")
        if i == 0 and f[0] == "subfamily":
            continue
        if len(f) >= 2:
            m[f[0]] = f[1]
    return m


def step2_dir(run):
    cands = sorted(glob.glob(os.path.join(run, "step2", "step2_output*")), key=os.path.getmtime)
    if not cands:
        raise SystemExit(f"no step2/step2_output* in {run}")
    return cands[-1]


def parse_runner(tag):
    d = {}
    for part in (tag or "").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def read_votes(path):
    """seq -> list of (cons, votes, bs_sum, rbs_sum, rcons, span, runner_span); span = (start, end) of the copy on cons in its best
    cycle, runner_span the same on rcons (step 2 since 2026-10-26, columns 7-10; None in older vote tables)"""
    out = defaultdict(list)
    sp = lambda a, b: (int(a), int(b)) if a.isdigit() and b.isdigit() and int(b) > 0 else None
    for line in opener(path):
        f = line.rstrip("\r\n").split("\t")
        if len(f) < 4:
            continue
        rb = float(f[4]) if len(f) > 4 and f[4] not in ("", "NA") else 0.0
        s = sp(f[6], f[7]) if len(f) > 7 else None
        r = sp(f[8], f[9]) if len(f) > 9 else None
        out[f[0]].append((f[1], int(f[2]), float(f[3]), rb, f[5] if len(f) > 5 else "NA", s, r))
    return out


PARTIAL_COVER = 0.60    # a copy covering less than this share of every consensus it competes for is a partial copy
PARTIAL_MIN_BP = 50     # shorter matches are not labelled
PARTIAL_RUNNER = 0.90   # a runner-up scoring at least this share of the winner competes too


def partial_call(rows, lens, thr, f_of, fam_name, fam_bs):
    """Partial copy (2026-10-26): a family-firm copy whose subfamily is open because it matches only a stretch of its
    consensuses, typically a unit that several composites share (rsi: the r1 head of r1_r3_P1, r1_r3_P2 and C11 without r1 in the
    bank). Competitors = the consensuses of the family that won a cycle, plus the best one's runner-up when it scores >= 90 % of it.
    Partial when the copy covers < 60 % of every competitor and >= 50 bp of the best. It must pass the best consensus's threshold
    scaled to the stretch (threshold x stretch length / consensus length; the threshold is 0.45 x the 10th-best full-length sum).
    Returns None (not partial, or spans unknown) or (label, best, coverage, stretch threshold, passes, {consensus: span})."""
    own = [r for r in rows if f_of(r[0]) == fam_name and r[1] > 0]
    if not own or any(r[5] is None for r in own):
        return None
    best = max(own, key=lambda r: (r[1], r[2]))
    comp = {r[0]: r[5] for r in own}
    if best[4] not in ("", "NA") and f_of(best[4]) == fam_name and best[6] is not None and best[2] > 0 \
            and best[3] / best[2] >= PARTIAL_RUNNER and best[4] not in comp:
        comp[best[4]] = best[6]
    cov = {c: (s[1] - s[0] + 1) / float(lens[c]) for c, s in comp.items() if lens.get(c)}
    bspan = best[5][1] - best[5][0] + 1
    if not cov or max(cov.values()) >= PARTIAL_COVER or bspan < PARTIAL_MIN_BP or best[0] not in thr:
        return None
    sthr = thr[best[0]] * bspan / float(lens[best[0]])
    order = [best[0]] + sorted(c for c in comp if c != best[0])
    label = ",".join("%s:%d-%d" % (c, comp[c][0], comp[c][1]) for c in order)
    return label, best[0], cov[best[0]], sthr, fam_bs >= sthr, comp


def group_partials(items):
    """Group partial copies into stretches. items: list of (family, {consensus: (start, end)}). Two copies of one family are
    linked when, on a consensus both match, their starts and their ends each differ by at most max(10, 10 % of the shorter
    span); a stretch = a connected group. Returns (group index per item, groups), a group = (family, [item indexes],
    {consensus: (median start, median end)} over the consensuses matched by at least half of its copies), groups sorted by
    family and size (largest first)."""
    par = list(range(len(items)))
    def find(i):
        while par[i] != i:
            par[i] = par[par[i]]
            i = par[i]
        return i
    by = defaultdict(list)
    for i, (F, comp) in enumerate(items):
        for c, (s, e) in comp.items():
            by[(F, c)].append((s, e, i))
    for key, v in by.items():
        v.sort()
        for a in range(len(v)):
            sa, ea, ia = v[a]
            for b in range(a + 1, len(v)):
                sb, eb, ib = v[b]
                tol = max(10, 0.1 * min(ea - sa + 1, eb - sb + 1))
                if sb - sa > tol:
                    break
                if abs(eb - ea) <= tol:
                    par[find(ia)] = find(ib)
    comps = defaultdict(list)
    for i in range(len(items)):
        comps[find(i)].append(i)
    med = lambda x: sorted(x)[len(x) // 2]
    groups = []
    for idx in comps.values():
        spans = defaultdict(list)
        for i in idx:
            for c, se in items[i][1].items():
                spans[c].append(se)
        keep = {c: (med([s for s, _ in v]), med([e for _, e in v])) for c, v in spans.items() if 2 * len(v) >= len(idx)}
        groups.append((items[idx[0]][0], sorted(idx), keep))
    groups.sort(key=lambda g: (g[0], -len(g[1]), min(g[1])))
    gi = [0] * len(items)
    for k, g in enumerate(groups):
        for i in g[1]:
            gi[i] = k
    return gi, groups


def stretch_label(spans):
    """'cons:start-end' of a stretch on each of its consensuses, the longest-matched one first"""
    order = sorted(spans, key=lambda c: (-(spans[c][1] - spans[c][0]), c))
    return ",".join("%s:%d-%d" % (c, spans[c][0], spans[c][1]) for c in order)


def merge_regions(rows):
    """bedtools merge -c 4 -o distinct (book-ended and overlapping features merged); rows = (chrom, start, end, name)"""
    rows = sorted(rows, key=lambda r: (r[0], r[1], r[2]))
    out = []
    for c, s, e, n in rows:
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
            if n not in out[-1][3]:
                out[-1][3].append(n)
        else:
            out.append([c, s, e, [n]])
    return [(c, s, e, ",".join(n)) for c, s, e, n in out]


TIE_CAP = 0.10          # the relative margin never grows beyond 10 % of the score


class ConsDist:
    """divergence (share of differing columns among the columns where both consensuses have a base, pairwise MAFFT alignment) of two bank consensuses, cached"""

    def __init__(self, bank, mafft="mafft"):
        import families as fm
        self.fm, self.mafft = fm, mafft
        self.seq = fm.read_fasta(bank)
        self.cache = {}

    def get(self, a, b):
        key = (a, b) if a <= b else (b, a)
        if key not in self.cache:
            if a not in self.seq or b not in self.seq:
                self.cache[key] = 0.0
            else:
                x, y = self.fm.align_mafft({a: self.seq[a], b: self.seq[b]}, self.mafft)
                cols = [(p, q) for p, q in zip(x.upper(), y.upper()) if p in "ACGT" and q in "ACGT"]     # shared columns only: an inserted element is not a difference in score per column
                self.cache[key] = sum(1 for p, q in cols if p != q) / float(len(cols)) if cols else 0.0
        return self.cache[key]


def rollup(run, fam_map_path, fam_cons, out, votes_path=None, bank=None, near_tie=0.02, tie_rel=0.5):
    fam = read_family_map(fam_map_path)
    f_of = lambda s: fam.get(s, s)
    s2 = step2_dir(run)
    os.makedirs(os.path.join(out, "step2", "step2_output"), exist_ok=True)
    os.makedirs(os.path.join(out, "genome.clean_step1", "searches"), exist_ok=True)
    shutil.copyfile(fam_cons, os.path.join(out, "consensuses.clean.fa"))
    src_ext = os.path.join(run, "genome.clean_step1", "extracted.fasta")
    dst_ext = os.path.join(out, "genome.clean_step1", "extracted.fasta")
    if os.path.exists(dst_ext):
        os.remove(dst_ext)
    try:
        os.symlink(os.path.abspath(src_ext), dst_ext)
    except (OSError, NotImplementedError):
        shutil.copyfile(src_ext, dst_ext)

    # hits and regions
    hits_in = os.path.join(run, "genome.clean_step1", "searches", "all_hits.labeled.bed")
    hits_out = os.path.join(out, "genome.clean_step1", "searches", "all_hits.labeled.bed")
    rows = []
    with open(hits_out, "w") as fo:
        for line in open(hits_in):
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 4:
                continue
            f[3] = f_of(f[3])
            fo.write("\t".join(f) + "\n")
            rows.append((f[0], int(f[1]), int(f[2]), f[3]))
    with open(os.path.join(out, "genome.clean_step1", "searches", "regions.by_subfam.bed"), "w") as fo:
        for c, s, e, n in merge_regions(rows):
            fo.write(f"{c}\t{s}\t{e}\t{n}\n")

    # thresholds per subfamily from the run itself
    thr = {}
    rows_full = []
    for i, line in enumerate(open(os.path.join(s2, "assignment_full.tsv"))):
        f = line.rstrip("\r\n").split("\t")
        if i == 0 and f[0] == "Sequence":
            continue
        if len(f) < 7:
            f += ["NA"] * (7 - len(f))
        rows_full.append(f)
        try:
            t = float(f[5])
            thr[f[1]] = max(thr.get(f[1], 0.0), t)
        except ValueError:
            pass
    fam_thr = defaultdict(list)
    for s, t in thr.items():
        if t > 0:
            fam_thr[f_of(s)].append(t)
    votes = read_votes(votes_path) if votes_path else {}
    dist = ConsDist(bank) if bank and tie_rel > 0 else None
    lens = {}
    if bank and os.path.exists(bank):
        import families as fm
        lens = {k: len(v) for k, v in fm.read_fasta(bank).items()}
    partial_rows = []

    stats = defaultdict(lambda: defaultdict(int))
    per_copy = []
    out_rows = []                 # assignment rows, written after the partial copies are grouped into stretches
    pending = []                  # (row index, partial_rows index) of the copies labelled partial
    for seq, sub, bs, nv, status, th, tag in rows_full:
        F = f_of(sub)
        new_status, new_bs, new_nv, new_th = status, bs, nv, th
        layer2 = ""
        rd = parse_runner(tag)
        runner = rd.get("runner", "NA")
        rbs = rd.get("runner_bs", "0")
        rratio = rd.get("runner_ratio", "0")
        if status == "assigned":
            stats[F]["sub_firm"] += 1
            # layer 2: a subfamily call whose runner-up is another subfamily of the same family is a near-tie when its score advantage is
            # below max(near_tie, tie_rel x divergence of the two consensuses)
            if runner != "NA" and runner != sub and f_of(runner) == F and f_of(sub) == F:
                try:
                    margin = 1.0 - float(rratio)
                except ValueError:
                    margin = 1.0
                dd = dist.get(sub, runner) if dist else 0.0
                if margin < min(TIE_CAP, max(near_tie, tie_rel * dd)):
                    layer2 = ";layer2=tie"
                    stats[F]["tie"] += 1
        elif seq in votes:
            fv = defaultdict(int)
            fb = defaultdict(float)
            fr = defaultdict(float)
            fc = {}
            for cons, v, b, r, rc, _s, _rs in votes[seq]:
                g = f_of(cons)
                fv[g] += v
                fb[g] += b
                if f_of(rc) != g and r >= fr[g]:
                    fr[g] = r
                    fc[g] = rc
            best = max(fv, key=lambda g: (fv[g], fb[g]))
            F = best
            new_nv = str(fv[best])
            new_bs = str(int(fb[best] + 0.5))
            lowest = min(fam_thr[best]) if fam_thr.get(best) else 0
            if fv[best] == 10 and fb[best] >= lowest:
                new_status = "assigned"
                new_th = str(int(lowest))
                layer2 = ";layer2=unresolved"
                pc = partial_call(votes[seq], lens, thr, f_of, best, fb[best]) if lens else None
                if pc:
                    label, pbest, pcov, sthr, ok, comp = pc
                    partial_rows.append([seq, best, pbest, label, "%.2f" % pcov, str(int(fb[best] + 0.5)), str(int(sthr + 0.5)),
                                         "partial" if ok else "below_stretch_threshold", "-", comp])
                    if ok:
                        layer2 = ";layer2=partial:"          # + stretch id, filled in after grouping
                        pending.append((len(out_rows), len(partial_rows) - 1))
                        stats[F]["partial"] += 1
                rc = fc.get(best, "NA")
                runner = f_of(rc) if rc != "NA" else "NA"
                rbs = str(int(fr[best] + 0.5))
                rratio = f"{(fr[best] / fb[best]) if fb[best] else 0:.4f}"
                stats[F]["family_only"] += 1
            elif fv[best] == 10:
                new_status = "rejected_low_bitscore"
        # runner of the same family does not count at family level
        if runner != "NA" and f_of(runner) == F:
            runner, rbs, rratio = "NA", "0", "0.0000"
        elif runner != "NA":
            runner = f_of(runner)
        if layer2 == ";layer2=tie":                       # the runner is in the same family: keep it visible as a tie, not as a LEAK
            new_tag = f"runner=NA;runner_bs=0;runner_ratio=0.0000{layer2};tie_with={rd.get('runner', 'NA')}"
        else:
            new_tag = "NA" if tag == "NA" and not layer2 else f"runner={runner};runner_bs={rbs};runner_ratio={rratio}{layer2}"
        out_rows.append([seq, F, str(new_bs), str(new_nv), new_status, str(new_th), new_tag])
        per_copy.append([seq, sub, F, status, new_status, new_nv, layer2.lstrip(";")])
        stats[F]["total"] += 1
        if new_status != "assigned":
            stats[F]["not_firm"] += 1
    # partial copies grouped into stretches; each copy's layer2 gets its stretch id (<family>_partial<n>, largest first)
    gi, groups = group_partials([(partial_rows[j][1], partial_rows[j][9]) for _, j in pending])
    sid, nth = [], defaultdict(int)
    for F_, idx, sp in groups:
        nth[F_] += 1
        sid.append("%s_partial%d" % (F_, nth[F_]))
    for k, (ri, j) in enumerate(pending):
        out_rows[ri][6] = out_rows[ri][6].replace("layer2=partial:", "layer2=partial:" + sid[gi[k]], 1)
        per_copy[ri][6] = per_copy[ri][6].replace("layer2=partial:", "layer2=partial:" + sid[gi[k]], 1)
        partial_rows[j][8] = sid[gi[k]]
    with open(os.path.join(out, "step2", "step2_output", "assignment_full.tsv"), "w") as full_out, \
            open(os.path.join(out, "step2", "step2_output", "all_votes.tsv"), "w") as votes_out:
        full_out.write("Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo\n")
        for r in out_rows:
            full_out.write("\t".join(r) + "\n")
            votes_out.write(f"{r[0]}\t{r[1]}\t{r[3]}\t{r[2]}\n")
    with open(os.path.join(out, "partial_stretches.tsv"), "w") as fo:
        fo.write("Stretch\tFamily\tCopies\tMedian_length\tSpans\n")
        for k, (F_, idx, sp) in enumerate(groups):
            L = []
            for i in idx:
                r = partial_rows[pending[i][1]]
                b0, b1 = r[9][r[2]]
                L.append(b1 - b0 + 1)
            L.sort()
            fo.write("%s\t%s\t%d\t%d\t%s\n" % (sid[k], F_, len(idx), L[len(L) // 2], stretch_label(sp)))
    with open(os.path.join(out, "family_view.tsv"), "w") as fo:
        fo.write("Sequence\tSubfamily\tFamily\tStatus_subfamily\tStatus_family\tFamily_votes\tLayer2\n")
        for r in per_copy:
            fo.write("\t".join(r) + "\n")
    # partial copies: family-firm copies that match only a stretch of their consensuses (partial_call)
    with open(os.path.join(out, "partial_copies.tsv"), "w") as fo:
        fo.write("Sequence\tFamily\tBest\tStretch\tCoverage_of_best\tFamily_score\tStretch_threshold\tStatus\tStretch_id\n")
        for r in partial_rows:
            fo.write("\t".join(r[:9]) + "\n")
    members = defaultdict(list)
    for s in set(list(fam.keys()) + list(thr.keys())):
        members[f_of(s)].append(s)
    with open(os.path.join(out, "family_view.summary.tsv"), "w") as fo:
        fo.write("family\tsubfamilies\tcopies\tsubfamily_firm\tsubfamily_near_tie\tfamily_firm_only\tfamily_firm_total\tnot_firm\n")
        for F in sorted(stats):
            st = stats[F]
            fo.write(f"{F}\t{','.join(sorted(members[F]))}\t{st['total']}\t{st['sub_firm']}\t{st['tie']}\t{st['family_only']}\t"
                     f"{st['sub_firm'] + st['family_only']}\t{st['not_firm']}\n")
    return stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--families", required=True, help="families.clean.tsv from tools/families.py")
    ap.add_argument("--family-cons", required=True, help="family_consensuses.fa from tools/families.py")
    ap.add_argument("--out", default=None)
    ap.add_argument("--votes", default=None, help="vote_summary.tsv[.gz] (seq, cons, votes, bs_sum, runner_bs_sum, runner_cons)")
    ap.add_argument("--bank", default=None, help="consensuses.clean.fa of the run; with it the near-tie margin also scales with the divergence of the two subfamily consensuses")
    ap.add_argument("--near-tie", type=float, default=0.02, help="absolute score margin below which two subfamilies of one family are a tie (default 0.02)")
    ap.add_argument("--tie-rel", type=float, default=0.5, help="x divergence of the two consensuses, also a margin, capped at 10 % (default 0.5; 0 = absolute only)")
    a = ap.parse_args(argv)
    out = a.out or os.path.join(a.run, "family_view")
    stats = rollup(a.run, a.families, a.family_cons, out, a.votes, a.bank, a.near_tie, a.tie_rel)
    tot = sum(s["total"] for s in stats.values())
    firm = sum(s["sub_firm"] + s["family_only"] for s in stats.values())
    print(f"family_view: {len(stats)} families, {tot} copies, {firm} family-firm ({sum(s['family_only'] for s in stats.values())} of them "
          f"subfamily-unresolved) -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
