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
    """seq -> list of (cons, votes, bs_sum, rbs_sum, rcons)"""
    out = defaultdict(list)
    for line in opener(path):
        f = line.rstrip("\r\n").split("\t")
        if len(f) < 4:
            continue
        rb = float(f[4]) if len(f) > 4 and f[4] not in ("", "NA") else 0.0
        out[f[0]].append((f[1], int(f[2]), float(f[3]), rb, f[5] if len(f) > 5 else "NA"))
    return out


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

    stats = defaultdict(lambda: defaultdict(int))
    per_copy = []
    full_out = open(os.path.join(out, "step2", "step2_output", "assignment_full.tsv"), "w")
    full_out.write("Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo\n")
    votes_out = open(os.path.join(out, "step2", "step2_output", "all_votes.tsv"), "w")
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
            for cons, v, b, r, rc in votes[seq]:
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
        full_out.write("\t".join([seq, F, str(new_bs), str(new_nv), new_status, str(new_th), new_tag]) + "\n")
        votes_out.write(f"{seq}\t{F}\t{new_nv}\t{new_bs}\n")
        per_copy.append((seq, sub, F, status, new_status, new_nv, layer2.lstrip(";")))
        stats[F]["total"] += 1
        if new_status != "assigned":
            stats[F]["not_firm"] += 1
    full_out.close()
    votes_out.close()
    with open(os.path.join(out, "family_view.tsv"), "w") as fo:
        fo.write("Sequence\tSubfamily\tFamily\tStatus_subfamily\tStatus_family\tFamily_votes\tLayer2\n")
        for r in per_copy:
            fo.write("\t".join(r) + "\n")
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
