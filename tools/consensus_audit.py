#!/usr/bin/env python3
"""Consensus audit: rebuild every consensus of the bank from its own assigned copies and compare with the bank sequence.

The rebuild is an independent check on a consensus, not a replacement for it: tools/vendor/sine_consensus.sh (the bootstrap
consensus builder of github.com/Toki-bio/SINE_consensus: seeded random subsamples, mini-consensuses accumulated and re-aligned
until stable; a base is called when at least 50 % of the copies carry it and no other base ties it - columns with a base
in fewer than 30 % of the copies are skipped first) is run on the copies assigned to each consensus, twice with
different seeds, and the result is aligned globally to the bank sequence.

The same copies also answer a second question that nothing else in the workflow asks: do they resemble EACH OTHER? Every
gate (the score cut, the flank test, the per-copy verdict) compares a copy with a consensus, so a consensus whose
composition alone draws hits can collect copies that each pass on their own and are unrelated to one another. The audit
therefore aligns a seeded sample of the assigned copies with each other (MAFFT, the same binary the rebuild uses) and
reports the median and the 10th-90th percentile of pairwise identity over the columns where both copies have a base,
next to the GC of the bank sequence, of the copies and of the rebuild. Numbers only: no threshold, no verdict change.

Per family (results/consensus_audit/summary.tsv):
  copies, bank_bp and bank_gc, median GC of the copies, copies_compared / copy_pairs and the pairwise identity between
  the copies (median, p10, p90), rebuilt_bp and GC (run 1, run 2), mismatches and gap columns against the bank,
  run1_vs_run2 mismatches, verdict
  MATCH     both rebuilds within MAX_MM mismatches of the bank and within the length tolerance
  SHORTER   a rebuild is shorter than the bank by more than the tolerance: the copies are fragmentary or the bank has a stretch
            (seed, longer copies) that fewer than half of the assigned copies carry
  LONGER    a rebuild is longer: a stretch is carried by at least half of the copies but is not in the bank (a longer version in the pool,
            or a tail the bank consensus lacks)
  DIVERGED  more than MAX_MM mismatches against the bank
  UNSTABLE  the two rebuilds differ from each other by more than MAX_MM mismatches
  SKIPPED   fewer than MIN_COPIES assigned copies (test on a species with more)
The verdict is about the rebuild only; the copy-vs-copy identity is reported beside it and never enters it.
Nothing is changed in the bank: the table is a decision input.

Usage: consensus_audit.py RUN_DIR [--jobs 8] [--min-copies 20] [--seeds 1,2] [--subsample 100]
                                  [--identity-sample 60] [--identity-seed 20261123]
RUN_DIR needs consensuses.clean.fa and results/assigned.fasta (ids "locus|family|bits").
"""
import argparse
import os
import random
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "vendor", "sine_consensus.sh")
MIN_COPIES = 20
MAX_MM = 4
LEN_TOL_FRAC = 0.05
LEN_TOL_BP = 8
IDENT_SAMPLE = 60        # copies drawn for the copy-against-copy alignment (60 copies = 1770 pairs)
IDENT_SEED = 20261123    # fixed, so the same assigned set always draws the same sample
IDENT_MIN_COPIES = 4     # fewer copies than this: no identity (6 pairs is the floor of a median)
IDENT_MIN_COLUMNS = 30   # a pair sharing fewer columns than this is not counted


def read_fa(path):
    d, n = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                n = line[1:].split()[0]
                d[n] = []
            elif n is not None:
                d[n].append(line)
    return {k: "".join(v).upper() for k, v in d.items()}


def split_assigned(path, outdir, families):
    """results/assigned.fasta -> outdir/<family>.fa; returns {family: copies}"""
    n = {f: 0 for f in families}
    out = {}
    cur = None
    try:
        for line in open(path):
            if line.startswith(">"):
                parts = line[1:].strip().split("|")
                fam = parts[1] if len(parts) > 1 else None
                cur = None
                if fam in n:
                    if fam not in out:
                        out[fam] = open(os.path.join(outdir, fam + ".fa"), "w")
                    cur = out[fam]
                    n[fam] += 1
            if cur is not None:
                cur.write(line)
    finally:
        for fh in out.values():
            fh.close()
    return n


def read_rows(path):
    """sequences of a FASTA in file order (a list, so two records with the same name are both kept)"""
    rows, cur = [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                cur = []
                rows.append(cur)
            elif cur is not None:
                cur.append(line)
    return ["".join(r).upper() for r in rows]


def gc_pct(seq):
    """GC % over the A/C/G/T of one sequence (gaps, N and everything else ignored); None when there is no base"""
    s = seq.upper()
    acgt = sum(s.count(c) for c in "ACGT")
    return 100.0 * (s.count("G") + s.count("C")) / acgt if acgt else None


def percentile(values, q):
    """nearest-rank percentile of a sorted list (no interpolation, so the number is exactly reproducible)"""
    if not values:
        return None
    k = int(round(q / 100.0 * (len(values) - 1)))
    return values[max(0, min(len(values) - 1, k))]


def pairwise_identities(rows):
    """identity % of every pair of aligned rows, over the columns where BOTH rows have a base; sorted ascending.

    A column with a gap or an N in either row is not counted in either the numerator or the denominator, so an
    insertion in one copy neither helps nor hurts the pair. Pairs sharing fewer than IDENT_MIN_COLUMNS columns
    are dropped (nothing to measure)."""
    keep = {"A": "A", "C": "C", "G": "G", "T": "T"}
    seqs = ["".join(keep.get(c, ".") for c in r.upper()) for r in rows]
    out = []
    for i, a in enumerate(seqs):
        for b in seqs[i + 1:]:
            shared = sum(1 for x, y in zip(a, b) if x != "." and y != ".")
            if shared < IDENT_MIN_COLUMNS:
                continue
            same = sum(1 for x, y in zip(a, b) if x == y and x != ".")
            out.append(100.0 * same / shared)
    out.sort()
    return out


def copy_identity(fam, workdir, outdir, sample, seed):
    """Do the copies assigned to one family resemble each other?

    Aligns a seeded sample of at most `sample` of them with each other (MAFFT, the same binary and the same
    --threadit 0 the rebuild uses, so the alignment is deterministic) and returns
    {copies_compared, copy_pairs, copy_identity_median, copy_identity_p10, copy_identity_p90, copy_gc_median},
    or None when there are too few copies or MAFFT does not finish. The alignment is kept, so the numbers can be
    read off a plate: results/consensus_audit/copies/<family>.sample.aln.fa"""
    d = read_fa(os.path.join(workdir, fam + ".fa"))
    names = sorted(d)                                        # sorted first: the draw must not depend on file order
    if len(names) < IDENT_MIN_COPIES:
        return None
    if len(names) > sample:
        names = sorted(random.Random(seed).sample(names, sample))
    gcs = sorted(g for g in (gc_pct(d[n]) for n in names) if g is not None)
    src = os.path.join(workdir, fam + ".sample.fa")
    with open(src, "w") as fh:
        for n in names:
            fh.write(">%s\n%s\n" % (n, d[n]))
    os.makedirs(outdir, exist_ok=True)
    aln = os.path.join(outdir, fam + ".sample.aln.fa")
    try:
        with open(aln, "w") as out:
            subprocess.run(["mafft", "--threadit", "0", "--auto", "--quiet", "--preservecase", src],
                           check=True, stdout=out, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, OSError):
        return None
    ids = pairwise_identities(read_rows(aln))
    if not ids:
        return None
    return {"copies_compared": len(names), "copy_pairs": len(ids),
            "copy_identity_median": "%.1f" % percentile(ids, 50),
            "copy_identity_p10": "%.1f" % percentile(ids, 10),
            "copy_identity_p90": "%.1f" % percentile(ids, 90),
            "copy_gc_median": "%.1f" % percentile(gcs, 50) if gcs else ""}


def align(a, b):
    """global alignment (match +1, mismatch -1, gap -1); returns (mismatches, gap_columns). a, b <= ~1500 bp"""
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0] = -i
    for j in range(1, m + 1):
        S[0][j] = -j
    for i in range(1, n + 1):
        ai = a[i - 1]
        Si, Sp = S[i], S[i - 1]
        for j in range(1, m + 1):
            Si[j] = max(Sp[j - 1] + (1 if ai == b[j - 1] else -1), Sp[j] - 1, Si[j - 1] - 1)
    i, j, mm, gp = n, m, 0, 0
    while i or j:
        if i and j and S[i][j] == S[i - 1][j - 1] + (1 if a[i - 1] == b[j - 1] else -1):
            mm += a[i - 1] != b[j - 1]
            i -= 1
            j -= 1
        elif i and S[i][j] == S[i - 1][j] - 1:
            gp += 1
            i -= 1
        else:
            gp += 1
            j -= 1
    return mm, gp


def verdict(bank_len, rebuilt, mms, mm12):
    """rebuilt: list of lengths of the successful rebuilds; mms: mismatches against the bank; mm12: between the two rebuilds"""
    if not rebuilt:
        return "FAILED"
    if max(mms) > MAX_MM:
        return "DIVERGED"
    if mm12 is not None and mm12 > MAX_MM:
        return "UNSTABLE"
    tol = max(LEN_TOL_BP, LEN_TOL_FRAC * bank_len)
    if min(rebuilt) < bank_len - tol:
        return "SHORTER"
    if max(rebuilt) > bank_len + tol:
        return "LONGER"
    return "MATCH"


def rebuild(fam, workdir, seed, subsample):
    d = os.path.join(workdir, "%s.s%d" % (fam, seed))
    os.makedirs(d, exist_ok=True)
    src = os.path.abspath(os.path.join(workdir, fam + ".fa"))
    try:
        subprocess.run(["bash", SCRIPT, "-n", str(subsample), "-r", str(seed), src], cwd=d, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, OSError):
        return None
    p = os.path.join(d, fam + "_consensus.fasta")
    if not os.path.exists(p):
        return None
    seqs = read_fa(p)
    return "".join(seqs.values()) if seqs else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--min-copies", type=int, default=MIN_COPIES)
    ap.add_argument("--seeds", default="1,2")
    ap.add_argument("--subsample", type=int, default=100)
    ap.add_argument("--identity-sample", type=int, default=IDENT_SAMPLE,
                    help="copies drawn per family for the copy-against-copy alignment (default %d)" % IDENT_SAMPLE)
    ap.add_argument("--identity-seed", type=int, default=IDENT_SEED,
                    help="seed of that draw; fixed so the numbers repeat (default %d)" % IDENT_SEED)
    ap.add_argument("--bank", default=None)
    a = ap.parse_args()
    run = os.path.abspath(a.run_dir)
    bank_path = a.bank or os.path.join(run, "consensuses.clean.fa")
    assigned = os.path.join(run, "results", "assigned.fasta")
    for f in (bank_path, assigned, SCRIPT):
        if not os.path.exists(f):
            print("consensus_audit: missing %s, skipped" % f)
            return 0
    if shutil.which("mafft") is None or shutil.which("gawk") is None:
        print("consensus_audit: mafft or gawk not on PATH, skipped")
        return 0
    seeds = [int(x) for x in a.seeds.split(",")]
    out = os.path.join(run, "results", "consensus_audit")
    work = os.path.join(out, "work")
    copies_dir = os.path.join(out, "copies")
    shutil.rmtree(work, ignore_errors=True)
    shutil.rmtree(copies_dir, ignore_errors=True)
    os.makedirs(work)
    bank = read_fa(bank_path)
    counts = split_assigned(assigned, work, set(bank))
    todo = [f for f in bank if counts.get(f, 0) >= a.min_copies]
    # the copies are compared with each other for every family that has any, including the ones too small to rebuild
    ident_todo = [f for f in bank if counts.get(f, 0) >= IDENT_MIN_COPIES]
    jobs = [(f, s) for f in todo for s in seeds]
    print("consensus_audit: %d consensuses, %d with >= %d copies, %d rebuilds, %d compared copy against copy"
          % (len(bank), len(todo), a.min_copies, len(jobs), len(ident_todo)), flush=True)
    res, ident = {}, {}
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {(f, s): ex.submit(rebuild, f, work, s, a.subsample) for f, s in jobs}
        ifuts = {f: ex.submit(copy_identity, f, work, copies_dir, a.identity_sample, a.identity_seed) for f in ident_todo}
        for k, fu in futs.items():
            res[k] = fu.result()
        for f, fu in ifuts.items():
            ident[f] = fu.result()
    cols = (["family", "copies", "bank_bp", "bank_gc", "copy_gc_median", "copies_compared", "copy_pairs",
             "copy_identity_median", "copy_identity_p10", "copy_identity_p90"]
            + sum([["rebuilt_bp_s%d" % s, "rebuilt_gc_s%d" % s, "mismatches_s%d" % s, "gap_columns_s%d" % s] for s in seeds], [])
            + ["rebuilds_differ_mm", "verdict"])
    rows = []
    with open(os.path.join(out, "rebuilt.fa"), "w") as fa:
        for fam, seq in bank.items():
            row = {"family": fam, "copies": counts.get(fam, 0), "bank_bp": len(seq)}
            g = gc_pct(seq)
            row["bank_gc"] = "" if g is None else "%.1f" % g
            row.update(ident.get(fam) or {})
            if fam not in todo:
                row["verdict"] = "SKIPPED"
                rows.append(row)
                continue
            got = []
            for s in seeds:
                r = res.get((fam, s))
                if r:
                    mm, gp = align(seq, r)
                    row["rebuilt_bp_s%d" % s], row["mismatches_s%d" % s], row["gap_columns_s%d" % s] = len(r), mm, gp
                    rg = gc_pct(r)
                    row["rebuilt_gc_s%d" % s] = "" if rg is None else "%.1f" % rg
                    got.append((r, mm))
                    fa.write(">%s_seed%d\n%s\n" % (fam, s, r))
            mm12 = align(got[0][0], got[1][0])[0] if len(got) >= 2 else None
            row["rebuilds_differ_mm"] = "" if mm12 is None else mm12
            row["verdict"] = verdict(len(seq), [len(r) for r, _ in got], [m for _, m in got], mm12)
            rows.append(row)
    with open(os.path.join(out, "summary.tsv"), "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    shutil.rmtree(work, ignore_errors=True)
    from collections import Counter
    print("consensus_audit: " + ", ".join("%s %d" % kv for kv in sorted(Counter(r["verdict"] for r in rows).items())) + " -> results/consensus_audit/summary.tsv")
    meds = sorted(float(r["copy_identity_median"]) for r in rows if r.get("copy_identity_median"))
    if meds:
        print("consensus_audit: identity between the copies of one family, median per family: %.1f - %.1f %% (%d families)"
              % (meds[0], meds[-1], len(meds)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
