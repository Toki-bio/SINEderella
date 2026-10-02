#!/usr/bin/env python3
"""Consensus audit: rebuild every consensus of the bank from its own assigned copies and compare with the bank sequence.

The rebuild is an independent check on a consensus, not a replacement for it: tools/vendor/sine_consensus.sh (the bootstrap
consensus builder of github.com/Toki-bio/SINE_consensus: seeded random subsamples, mini-consensuses accumulated and re-aligned
until stable, a base is called when >= 30 % of the copies carry it) is run on the copies assigned to each consensus, twice with
different seeds, and the result is aligned globally to the bank sequence.

Per family (results/consensus_audit/summary.tsv):
  copies, bank_bp, rebuilt_bp (run 1, run 2), mismatches and gap columns against the bank, run1_vs_run2 mismatches, verdict
  MATCH     both rebuilds within MAX_MM mismatches of the bank and within the length tolerance
  SHORTER   a rebuild is shorter than the bank by more than the tolerance: the copies are fragmentary or the bank has a stretch
            (seed, longer copies) that fewer than 30 % of the assigned copies carry
  LONGER    a rebuild is longer: a stretch is carried by >= 30 % of the copies but is not in the bank (a longer version in the pool,
            or a tail the bank consensus lacks)
  DIVERGED  more than MAX_MM mismatches against the bank
  UNSTABLE  the two rebuilds differ from each other by more than MAX_MM mismatches
  SKIPPED   fewer than MIN_COPIES assigned copies (test on a species with more)
Nothing is changed in the bank: the table is a decision input.

Usage: consensus_audit.py RUN_DIR [--jobs 8] [--min-copies 20] [--seeds 1,2] [--subsample 100]
RUN_DIR needs consensuses.clean.fa and results/assigned.fasta (ids "locus|family|bits").
"""
import argparse
import os
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
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    bank = read_fa(bank_path)
    counts = split_assigned(assigned, work, set(bank))
    todo = [f for f in bank if counts.get(f, 0) >= a.min_copies]
    jobs = [(f, s) for f in todo for s in seeds]
    print("consensus_audit: %d consensuses, %d with >= %d copies, %d rebuilds" % (len(bank), len(todo), a.min_copies, len(jobs)), flush=True)
    res = {}
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {(f, s): ex.submit(rebuild, f, work, s, a.subsample) for f, s in jobs}
        for k, fu in futs.items():
            res[k] = fu.result()
    cols = ["family", "copies", "bank_bp"] + sum([["rebuilt_bp_s%d" % s, "mismatches_s%d" % s, "gap_columns_s%d" % s] for s in seeds], []) + ["rebuilds_differ_mm", "verdict"]
    rows = []
    with open(os.path.join(out, "rebuilt.fa"), "w") as fa:
        for fam, seq in bank.items():
            row = {"family": fam, "copies": counts.get(fam, 0), "bank_bp": len(seq)}
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
