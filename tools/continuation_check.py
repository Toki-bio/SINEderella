#!/usr/bin/env python3
"""continuation_check.py RUN_DIR: do the copies of a consensus continue into another consensus past their ends?

Found on rsi (2026-10-22): 535 of the 761 r10_19seqs copies did not end where r10 ends. Their "right flank" was the r8 module
of the composite r10_r8_P18 (P18 positions 293-423), because those copies are a composite the bank did not have: r10 joined
directly to the r8 module, without P18's middle piece. The vote gave them to the r10 unit, their copy coordinates stopped at
r10's end, the flank tests read SINE sequence as flank, and 216 more such copies matched neither consensus well enough and
were left unassigned. Nothing in the workflow reported it; the plate builder only noticed that "copies stay similar past the
flank end" and widened the plate.

This step asks, per consensus, separately for its assigned copies and for the unassigned copies whose best (too weak) match
it was: read each copy with WINDOW bp on both sides, in its own orientation; find where its own consensus matches inside it;
search what lies past each end of that match (the rest of the locus plus the margin) against the bank (ssearch36, forward
strand only), and keep the best hit that starts at the junction (within JUNC bp of the end of the own match) and is at least
MINLEN bp long. Reading past the own match, not past the locus, matters: in rsi the long r10 loci already covered r10 and the
r8 module (about 280 bp), so reading past the locus end found only genomic flank (first version of this tool, 2026-10-22).
inside_locus tells the two cases apart: a continuation mostly inside the locus = the copy was labelled with a shorter
consensus; mostly outside = the locus coordinates stop early. When at least MIN_SHARE of the sampled copies (and
at least MIN_N copies) continue into the same consensus, the row is flagged. For a flag into another consensus a candidate
consensus of "copy + continuation" is built from those copies with tools/vendor/sine_consensus.sh (two seeds, as the
consensus audit), and compared with the bank: if the bank already holds it, the problem is the assignment; if not, a
consensus is missing. A continuation into the same consensus (dimers, tandem copies) is reported, not built.

Nothing is changed: the run's copies, verdicts and bank stay as they are. The table goes on the report.

Out (RUN_DIR/results/continuation/):
  continuation.tsv   one row per (consensus, group, side, target) with at least MIN_N copies continuing:
                     consensus group side sampled target copies share junction_bp target_start target_end identity
                     inside_locus flag candidate candidate_bp seeds_differ_mm bank_match
                     (junction_bp, target_start/end, identity, inside_locus: medians over the continuing copies)
  candidates.fa      the candidate consensuses (seed 1)
  <candidate>.copies.fa  the copy + continuation sequences each candidate was built from
Prints one line: CONTINUATION groups=.. flagged=.. candidates=.. -> results/continuation/continuation.tsv
"""
import argparse
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "vendor", "sine_consensus.sh")
sys.path.insert(0, HERE)
from consensus_audit import align, read_fa  # noqa: E402  (same global alignment the audit uses)

WINDOW = 150        # bp read past each end
JUNC = 20           # a continuation starts within this many bp of the copy's end
MINLEN = 50         # and is at least this long
MIN_COPIES = 20     # groups with fewer copies are not tested
SAMPLE = 400        # copies per group (seeded)
MIN_SHARE = 0.20    # flag when at least this share of the sampled copies continue into one consensus
MIN_N = 10          # ... and at least this many copies
BUILD_MAX = 200     # copies per candidate build
SUBSAMPLE = 100     # sine_consensus.sh -n
EVALUE = "1e-5"
UNASSIGNED = ("rejected_low_bitscore", "no_unanimous")
LOCUS = re.compile(r"^(.+):(\d+)-(\d+)\(([+-])\)$")


def read_lengths(fai):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}


def groups_from(assign_path):
    """{(consensus, group): [(contig, start, end, strand, locus)]}; group = assigned | unassigned (best match too weak).
    Loci with a mixed strand ('+,-') have no single orientation and are skipped."""
    g = defaultdict(list)
    with open(assign_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        if head[:2] != ["Sequence", "Subfamily"] or len(head) < 5:
            raise SystemExit("continuation_check: %s: unexpected header %s" % (assign_path, head[:5]))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            m = LOCUS.match(f[0])
            if not m:
                continue
            st = f[4].strip()
            grp = "assigned" if st == "assigned" else ("unassigned" if st in UNASSIGNED else None)
            if grp is None:
                continue
            g[(f[1], grp)].append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), f[0]))
    return g


def getfasta(genome, bed_rows, tmp):
    """bed_rows: [(ctg, s, e, name, strand)] -> {name: seq}, strand-aware (bedtools getfasta -s)"""
    if not bed_rows:
        return {}
    bed = os.path.join(tmp, "q.bed")
    with open(bed, "w") as fh:
        for r in bed_rows:
            fh.write("%s\t%d\t%d\t%s\t0\t%s\n" % r)
    out = subprocess.run(["bedtools", "getfasta", "-fi", genome, "-bed", bed, "-s", "-nameOnly", "-tab"],
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("bedtools getfasta failed: " + out.stderr.strip()[:300])
    seqs = {}
    for line in out.stdout.splitlines():
        n, sq = line.split("\t")
        seqs[re.sub(r"\([+-]\)$", "", n)] = sq.upper()
    return seqs


def ssearch(query_fa, lib_fa, threads, nbest=5):
    """ssearch36 -m 8, forward strand of the queries only; a failed call is an error, not 'no hits' (ssearch36 can exit 0
    with '***' on stderr when it cannot read a file)."""
    out = subprocess.run(["ssearch36", "-m", "8", "-3", "-E", EVALUE, "-b", str(nbest), "-d", "0", "-T", str(threads), query_fa, lib_fa],
                         capture_output=True, text=True)
    if out.returncode != 0 or "***" in out.stderr:
        raise RuntimeError("ssearch36 failed: " + out.stderr.strip()[:300])
    hits = defaultdict(list)
    for line in out.stdout.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        hits[f[0]].append((float(f[11]), f[1], float(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8]), int(f[9])))
    return hits   # name -> [(bits, subject, identity, alnlen, qstart, qend, sstart, send)]


def junction_hit(side, hits, wlen):
    """best hit that starts at the junction: 3' window -> qstart <= JUNC; 5' window -> qend >= wlen - JUNC + 1"""
    ok = []
    for h in hits:
        bits, subj, ident, alen, qs, qe, ss, se = h
        if alen < MINLEN or se < ss:
            continue
        if (side == "3" and qs <= JUNC) or (side == "5" and qe >= wlen - JUNC + 1):
            ok.append(h)
    return max(ok) if ok else None


def build(name, seqs, workdir, seed):
    d = os.path.join(workdir, "%s.s%d" % (name, seed))
    os.makedirs(d, exist_ok=True)
    src = os.path.join(workdir, name + ".fa")
    r = subprocess.run(["bash", SCRIPT, "-n", str(SUBSAMPLE), "-r", str(seed), src], cwd=d, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    p = os.path.join(d, name + "_consensus.fasta")
    if r.returncode != 0 or not os.path.exists(p):
        return None
    s = read_fa(p)
    return "".join(s.values()) if s else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-build", action="store_true", help="report only, build no candidate")
    a = ap.parse_args()
    run = os.path.abspath(a.run_dir)
    genome = os.path.join(run, "genome.clean.fa")
    bank_path = os.path.join(run, "consensuses.clean.fa")
    assign = os.path.join(run, "results", "assignment_full.tsv")
    for f in (genome, genome + ".fai", bank_path, assign):
        if not os.path.exists(f):
            raise SystemExit("continuation_check: missing %s" % f)
    for tool in ("ssearch36", "bedtools") + (() if a.no_build else ("mafft",)):
        if shutil.which(tool) is None:
            raise SystemExit("continuation_check: %s not on PATH" % tool)
    out = os.path.join(run, "results", "continuation")
    work = os.path.join(out, "work")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(work)
    lens = read_lengths(genome + ".fai")
    bank = read_fa(bank_path)
    groups = groups_from(assign)
    rnd = random.Random(a.seed)

    # 1) a seeded sample of every group; each copy is read with WINDOW bp on both sides, in its own orientation
    sample, padded, left = {}, {}, {}
    beds = []
    for key in sorted(groups):
        cps = [c for c in groups[key] if c[0] in lens]
        if len(cps) < MIN_COPIES:
            continue
        cps = sorted(cps, key=lambda c: c[4])
        if len(cps) > SAMPLE:
            cps = rnd.sample(cps, SAMPLE)
        gi = len(sample)
        sample[key] = cps
        for i, (ctg, st, en, strand, _) in enumerate(cps):
            ps, pe = max(0, st - WINDOW), min(lens[ctg], en + WINDOW)
            n = "g%d_c%d" % (gi, i)
            beds.append((ctg, ps, pe, n, strand))
            left[n] = (st - ps) if strand == "+" else (pe - en)      # bp before the locus in the copy's orientation
    keys = list(sample)
    padded = getfasta(genome, beds, work)

    # 2) where does each copy's own consensus match inside it? (one search per group, forward strand)
    own = {}                                                         # n -> (sstart, send) in the padded copy, 1-based
    for gi, key in enumerate(keys):
        if key[0] not in bank:
            continue
        q, lib = os.path.join(work, "own_q.fa"), os.path.join(work, "own_lib.fa")
        with open(q, "w") as fh:
            fh.write(">%s\n%s\n" % (key[0], bank[key[0]]))
        with open(lib, "w") as fh:
            for i in range(len(sample[key])):
                n = "g%d_c%d" % (gi, i)
                if n in padded:
                    fh.write(">%s\n%s\n" % (n, padded[n]))
        for h in ssearch(q, lib, a.threads, nbest=len(sample[key]) + 10).get(key[0], []):
            bits, subj, ident, alen, qs, qe, ss, se = h
            if se > ss and (subj not in own or bits > own[subj][2]):
                own[subj] = (ss, se, bits, qs, qe)   # qs-qe: the matched part of the own consensus

    # 3) the sequence past each end of that match (the rest of the locus plus the margin), searched against the bank
    tails = {}
    for n, (ss, se, *_) in own.items():
        sq = padded[n]
        tails[n + "_3"] = sq[se:]
        tails[n + "_5"] = sq[:ss - 1]
    qfa = os.path.join(work, "tails.fa")
    with open(qfa, "w") as fh:
        for n, sq in tails.items():
            if len(sq) >= MINLEN:
                fh.write(">%s\n%s\n" % (n, sq))
    # every bank consensus may be reported (-b): with only the top five, a junction hit to a sixth was never seen when
    # related consensuses scored higher further from the junction (GLM audit 2026-10-22)
    hits = ssearch(qfa, bank_path, a.threads, nbest=len(bank)) if os.path.getsize(qfa) else {}

    # 4) per group, side and target: how many copies continue, and how much of the continuation lies inside the locus
    rows = []
    cont = defaultdict(list)   # (key, side, target) -> [(copy index, hit, bp of the continuation inside the locus)]
    for tn, sq in tails.items():
        if len(sq) < MINLEN:
            continue
        g, c, side = tn.split("_")
        n = g + "_" + c
        h = junction_hit(side, hits.get(tn, []), len(sq))
        if not h:
            continue
        cp = sample[keys[int(g[1:])]][int(c[1:])]
        locus_len = cp[2] - cp[1]
        ss, se = own[n][:2]
        if side == "3":   # continuation spans padded positions se+qs .. se+qe; the locus ends at left + locus_len
            inside = max(0, min(se + h[5], left[n] + locus_len) - max(se + h[4] - 1, left[n]))
        else:             # continuation spans padded positions qs .. qe; the locus starts after left[n]
            inside = max(0, min(h[5], left[n] + locus_len) - max(h[4] - 1, left[n]))
        cont[(keys[int(g[1:])], side, h[1])].append((int(c[1:]), h, inside))
    for (key, side, target), lst in sorted(cont.items()):
        n = len(lst)
        if n < MIN_N:
            continue
        gi = keys.index(key)
        tested = sum(1 for x in own if x.startswith("g%d_" % gi))   # copies in which the own consensus was found
        share = n / tested
        hs = [h for _, h, _ in lst]
        jb = statistics.median((h[4] - 1) if side == "3" else (len(tails["g%d_c%d_5" % (gi, i)]) - h[5]) for i, h, _ in lst)
        span = statistics.median(h[5] - h[4] + 1 for h in hs)
        ins = statistics.median(x for _, _, x in lst)
        # where the own match stops on the own consensus (3') or starts (5'): a continuation into the own consensus that
        # resumes well past that point is an internal deletion of it (P without its middle piece), not a dimer
        oe = statistics.median(own["g%d_c%d" % (gi, i)][4 if side == "3" else 3] for i, _, _ in lst)
        rows.append({"consensus": key[0], "group": key[1], "side": side + "'", "sampled": len(sample[key]), "tested": tested, "target": target,
                     "copies": n, "share": round(share, 3), "junction_bp": int(jb),
                     "target_start": int(statistics.median(h[6] for h in hs)), "target_end": int(statistics.median(h[7] for h in hs)),
                     "identity": round(statistics.median(h[2] for h in hs), 1),
                     "inside_locus": "%d of %d bp" % (ins, span),
                     "flag": "yes" if share >= MIN_SHARE else "no", "lst": lst, "key": key, "side0": side, "own_edge": int(oe)})
    nf = sum(len(v) for v in sample.values()) - len(own)
    if nf:
        print("continuation_check: own consensus not found in %d sampled copies (skipped)" % nf)

    # 3) candidates for flags into another consensus
    ncand = 0
    with open(os.path.join(out, "candidates.fa"), "w") as cfa:
        for r in rows:
            r.update(candidate="", candidate_bp="", seeds_differ_mm="", bank_match="")
            if r["flag"] != "yes":
                continue
            deletion = ""
            if r["target"] == r["consensus"]:
                if r["side0"] == "3" and r["target_start"] > r["own_edge"] + JUNC:
                    deletion = "deletion inside %s (own match ends at %d, resumes at %d); " % (r["consensus"], r["own_edge"], r["target_start"])
                elif r["side0"] == "5" and r["target_end"] < r["own_edge"] - JUNC:
                    deletion = "deletion inside %s (continuation ends at %d, own match resumes at %d); " % (r["consensus"], r["target_end"], r["own_edge"])
                else:
                    r["bank_match"] = "continues into itself (dimer or tandem copies): not built"
                    continue
            if a.no_build:
                continue
            name = ("%s__3p__%s_%d-%d" if r["side0"] == "3" else "%s__5p__%s_%d-%d") % (
                r["consensus"], r["target"], r["target_start"], r["target_end"])
            name = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
            cps = sample[r["key"]]
            pick = sorted(r["lst"], key=lambda x: cps[x[0]][4])[:BUILD_MAX]
            gi = keys.index(r["key"])
            seqs = []
            for i, h, _ in pick:     # own-consensus match + the continuation, cut from the padded copy
                n = "g%d_c%d" % (gi, i)
                ss, se = own[n][:2]
                seqs.append(padded[n][ss - 1:se + h[5]] if r["side0"] == "3" else padded[n][h[4] - 1:se])
            with open(os.path.join(work, name + ".fa"), "w") as fh, open(os.path.join(out, name + ".copies.fa"), "w") as cf:
                for i, s in enumerate(seqs):
                    fh.write(">s%d\n%s\n" % (i, s))
                    cf.write(">%s\n%s\n" % (cps[pick[i][0]][4], s))
            c1 = build(name, seqs, work, 1)
            c2 = build(name, seqs, work, 2)
            if not c1:
                r["bank_match"] = "candidate build failed"
                continue
            ncand += 1
            r["candidate"], r["candidate_bp"] = name, len(c1)
            r["seeds_differ_mm"] = align(c1, c2)[0] if c2 else "second seed failed"
            cfa.write(">%s\n%s\n" % (name, c1))
            cq = os.path.join(work, "candidate_query.fa")   # ssearch36 misreads file names like ..._291-420.cand.fa
            with open(cq, "w") as fh:
                fh.write(">%s\n%s\n" % (name, c1))
            bh = ssearch(cq, bank_path, a.threads, nbest=len(bank)).get(name, [])
            best = max(bh) if bh else None
            # in the bank = the whole candidate matches one consensus end to end: global alignment (as the consensus audit) with
            # few mismatches and gap columns. ssearch36's alignment length counts gaps, so "covers 90 % of both" passed A+B
            # against A+M+B although the whole middle piece M (140 bp) was missing (test_continuation_check).
            mm, gp = align(c1, bank[best[1]]) if best else (0, 0)
            for h in bh:   # the end-to-end match need not be the best local scorer: every consensus with a hit is aligned
                m2, g2 = align(c1, bank[h[1]])
                if m2 + g2 < mm + gp:
                    best, mm, gp = h, m2, g2
            if best and mm + gp <= max(8, 0.05 * len(c1)):
                r["bank_match"] = deletion + "already in the bank as %s (%d mismatches, %d gap columns): the copies were assigned elsewhere" % (best[1], mm, gp)
            elif best:
                r["bank_match"] = deletion + "not in the bank (closest %s: %d mismatches and %d gap columns end to end; local match at its positions %d-%d)" % (
                    best[1], mm, gp, best[6], best[7])
            else:
                r["bank_match"] = deletion + "not in the bank (no hit)"
    cols = ["consensus", "group", "side", "sampled", "tested", "target", "copies", "share", "junction_bp", "target_start", "target_end",
            "identity", "inside_locus", "flag", "candidate", "candidate_bp", "seeds_differ_mm", "bank_match"]
    with open(os.path.join(out, "continuation.tsv"), "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["flag"] != "yes", -r["share"])):
            o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    shutil.rmtree(work, ignore_errors=True)
    print("CONTINUATION groups=%d flagged=%d candidates=%d -> results/continuation/continuation.tsv" % (
        len(sample), sum(r["flag"] == "yes" for r in rows), ncand))
    return 0


if __name__ == "__main__":
    sys.exit(main())
