#!/usr/bin/env python3
"""continuation_check.py RUN_DIR: do the copies of a consensus continue into another consensus past their ends?

Found on rsi (2026-10-22): 535 of the 761 r10_19seqs copies did not end where r10 ends. Their "right flank" was the r8 module
of the composite r10_r8_P18 (P18 positions 293-423), because those copies are a composite the bank did not have: r10 joined
directly to the r8 module, without P18's middle piece. The vote gave them to the r10 unit, their copy coordinates stopped at
r10's end, the flank tests read SINE sequence as flank, and 216 more such copies matched neither consensus well enough and
were left unassigned. Nothing in the workflow reported it; the plate builder only noticed that "copies stay similar past the
flank end" and widened the plate.

This step asks, per consensus, separately for its assigned copies and for the unassigned copies whose best (too weak) match
it was: take the WINDOW bp right past each copy's 3' end and right before its 5' end (strand-aware, read in the copy's
orientation), search them against the bank (ssearch36, forward strand only), and keep, per window, the best hit that starts
at the junction (within JUNC bp of it) and is at least MINLEN bp long. When at least MIN_SHARE of the sampled copies (and
at least MIN_N copies) continue into the same consensus, the row is flagged. For a flag into another consensus a candidate
consensus of "copy + continuation" is built from those copies with tools/vendor/sine_consensus.sh (two seeds, as the
consensus audit), and compared with the bank: if the bank already holds it, the problem is the assignment; if not, a
consensus is missing. A continuation into the same consensus (dimers, tandem copies) is reported, not built.

Nothing is changed: the run's copies, verdicts and bank stay as they are. The table goes on the report.

Out (RUN_DIR/results/continuation/):
  continuation.tsv   one row per (consensus, group, side, target) with at least MIN_N copies continuing:
                     consensus group side sampled target copies share junction_bp target_start target_end identity flag
                     candidate candidate_bp seeds_differ_mm bank_match
                     (junction_bp, target_start/end, identity: medians over the continuing copies)
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


def windows(copy, clen):
    """BED intervals (start, end) of the 3' and 5' windows in genome coordinates; reading them with the copy's strand gives
    sequence in the copy's orientation: 3' window starts at the junction, 5' window ends at it."""
    ctg, s, e, strand, _ = copy
    right = (e, min(clen, e + WINDOW))
    left = (max(0, s - WINDOW), s)
    return {"3": right, "5": left} if strand == "+" else {"3": left, "5": right}


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


def ssearch(query_fa, lib_fa, threads):
    """ssearch36 -m 8, forward strand of the queries only; a failed call is an error, not 'no hits' (ssearch36 can exit 0
    with '***' on stderr when it cannot read a file)."""
    out = subprocess.run(["ssearch36", "-m", "8", "-3", "-E", EVALUE, "-b", "5", "-T", str(threads), query_fa, lib_fa],
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

    # 1) windows of a seeded sample of every group, one search for all
    sample = {}
    beds = []
    for key in sorted(groups):
        cps = [c for c in groups[key] if c[0] in lens]
        if len(cps) < MIN_COPIES:
            continue
        cps = sorted(cps, key=lambda c: c[4])
        if len(cps) > SAMPLE:
            cps = rnd.sample(cps, SAMPLE)
        sample[key] = cps
        for i, c in enumerate(cps):
            for side, (ws, we) in windows(c, lens[c[0]]).items():
                if we - ws >= MINLEN:
                    beds.append((c[0], ws, we, "g%d_c%d_%s" % (len(sample) - 1, i, side), c[3]))
    keys = list(sample)
    wseq = getfasta(genome, beds, work)
    qfa = os.path.join(work, "windows.fa")
    with open(qfa, "w") as fh:
        for n, s in wseq.items():
            fh.write(">%s\n%s\n" % (n, s))
    hits = ssearch(qfa, bank_path, a.threads) if wseq else {}

    # 2) per group, side and target: how many copies continue
    rows = []
    cont = defaultdict(list)   # (key, side, target) -> [(copy index, hit)]
    for n, s in wseq.items():
        gi, ci, side = n.split("_")
        h = junction_hit(side, hits.get(n, []), len(s))
        if h:
            cont[(keys[int(gi[1:])], side, h[1])].append((int(ci[1:]), h))
    for (key, side, target), lst in sorted(cont.items()):
        n = len(lst)
        if n < MIN_N:
            continue
        share = n / len(sample[key])
        hs = [h for _, h in lst]
        jb = statistics.median((h[4] - 1) if side == "3" else (len(wseq["g%d_c%d_%s" % (keys.index(key), i, side)]) - h[5])
                               for i, h in lst)
        rows.append({"consensus": key[0], "group": key[1], "side": side + "'", "sampled": len(sample[key]), "target": target,
                     "copies": n, "share": round(share, 3), "junction_bp": int(jb),
                     "target_start": int(statistics.median(h[6] for h in hs)), "target_end": int(statistics.median(h[7] for h in hs)),
                     "identity": round(statistics.median(h[2] for h in hs), 1),
                     "flag": "yes" if share >= MIN_SHARE else "no", "lst": lst, "key": key, "side0": side})

    # 3) candidates for flags into another consensus
    ncand = 0
    with open(os.path.join(out, "candidates.fa"), "w") as cfa:
        for r in rows:
            r.update(candidate="", candidate_bp="", seeds_differ_mm="", bank_match="")
            if r["flag"] != "yes":
                continue
            if r["target"] == r["consensus"]:
                r["bank_match"] = "continues into itself (dimer or tandem copies): not built"
                continue
            if a.no_build:
                continue
            name = ("%s__3p__%s_%d-%d" if r["side0"] == "3" else "%s__5p__%s_%d-%d") % (
                r["consensus"], r["target"], r["target_start"], r["target_end"])
            name = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
            cps = sample[r["key"]]
            pick = sorted(r["lst"], key=lambda x: cps[x[0]][4])[:BUILD_MAX]
            el = getfasta(genome, [(cps[i][0], cps[i][1], cps[i][2], "e%d" % i, cps[i][3]) for i, _ in pick], work)
            seqs = []
            for i, h in pick:
                w = wseq["g%d_c%d_%s" % (keys.index(r["key"]), i, r["side0"])]
                e = el.get("e%d" % i, "")
                seqs.append(e + w[:h[5]] if r["side0"] == "3" else w[h[4] - 1:] + e)
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
            bh = ssearch(cq, bank_path, a.threads).get(name, [])
            best = max(bh) if bh else None
            # in the bank = the whole candidate matches one consensus end to end: global alignment (as the consensus audit) with
            # few mismatches and gap columns. ssearch36's alignment length counts gaps, so "covers 90 % of both" passed A+B
            # against A+M+B although the whole middle piece M (140 bp) was missing (test_continuation_check).
            mm, gp = align(c1, bank[best[1]]) if best else (0, 0)
            if best and mm + gp <= max(8, 0.05 * len(c1)):
                r["bank_match"] = "already in the bank as %s (%d mismatches, %d gap columns): the copies were assigned elsewhere" % (best[1], mm, gp)
            elif best:
                r["bank_match"] = "not in the bank (closest %s: %d mismatches and %d gap columns end to end; local match at its positions %d-%d)" % (
                    best[1], mm, gp, best[6], best[7])
            else:
                r["bank_match"] = "not in the bank (no hit)"
    cols = ["consensus", "group", "side", "sampled", "target", "copies", "share", "junction_bp", "target_start", "target_end",
            "identity", "flag", "candidate", "candidate_bp", "seeds_differ_mm", "bank_match"]
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
