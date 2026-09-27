#!/usr/bin/env python3
"""Put the consensus exactly as searched onto each published plate, as row 2.

Usage: add_seed_row.py BANK.fa SPECIES_CODE ALIGNMENT.aln.fa [...] [--threads N]

The publish stage rebuilds row 1 of every plate from the copies and moves its edge; the seed that
the genome was actually searched with is otherwise gone from the plate. This adds it back, named
<subfamily>_seed_as_searched, directly under row 1, so both can be read together.

The seed is added with `mafft --add` into the finished alignment: existing rows are not realigned,
only gap columns are inserted where the seed needs them, so packed flanks stay packed. The file is
rewritten only if every original row is unchanged once those inserted all-gap columns are removed;
otherwise it is left as it was and a warning is printed. A plate that already has the row is
skipped, so the step can be re-run. If the plate is in the reverse-complement orientation of the
seed (it should not be after correct_published_aln.py), the reverse complement is added and the row
name says so. Run this after the report and its verdict columns are built: the verdicts are computed
on the plates without this row.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
KINDS = ("_top100.aln.fa", "_rand100.aln.fa", "_subfam.aln.fa")


def read_fa(path):
    names, seqs = [], []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            names.append(line[1:].strip())
            seqs.append([])
        elif seqs:
            seqs[-1].append(line)
    return names, ["".join(s) for s in seqs]


def write_fa(path, names, seqs):
    with open(path, "w", encoding="utf-8") as fh:
        for n, s in zip(names, seqs):
            fh.write(">%s\n" % n)
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")


def kshare(a, b, k=8):
    ka = {a[i:i + k] for i in range(len(a) - k + 1)}
    if not ka:
        return 0.0
    kb = {b[i:i + k] for i in range(len(b) - k + 1)}
    return len(ka & kb) / float(len(ka))


def ungap(s):
    return s.replace("-", "").replace(".", "")


def subfamily_of(path, code):
    b = os.path.basename(path)
    for k in KINDS:
        if b.endswith(k):
            b = b[:-len(k)]
            break
    else:
        return None
    return b[len(code) + 1:] if b.startswith(code + "_") else b


def add_one(path, bank, code, threads):
    sf = subfamily_of(path, code)
    if sf is None or sf not in bank:
        return "skip (no seed %r in bank)" % sf
    names, seqs = read_fa(path)
    if any(n.startswith(sf + "_seed_as_searched") for n in names):
        return "skip (already has the seed row)"
    if len(seqs) < 2:
        return "skip (fewer than 2 rows)"
    seed = bank[sf].upper()
    body = ungap(seqs[0]).upper()
    rc = seed.translate(COMP)[::-1]
    label = sf + "_seed_as_searched"
    if kshare(rc, body) > 2 * kshare(seed, body) and kshare(rc, body) >= 0.10:
        seed, label = rc, label + "_revcomp"

    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as td:
        aln, add = os.path.join(td, "aln.fa"), os.path.join(td, "seed.fa")
        # mafft names must be unique and simple; map back afterwards
        write_fa(aln, ["r%d" % i for i in range(len(seqs))], seqs)
        write_fa(add, ["seedrow"], [seed])
        out = subprocess.run(["mafft", "--add", add, "--preservecase", "--quiet", "--thread", str(threads), aln],
                             capture_output=True, text=True)
        if out.returncode != 0:
            return "WARN mafft failed: %s" % out.stderr.strip()[:200]
        n2, s2 = read_fa_text(out.stdout)
    if len(s2) != len(seqs) + 1 or n2[-1] != "seedrow":
        return "WARN unexpected mafft output; left unchanged"
    width = len(s2[0])
    old = s2[:-1]
    # columns where every original row is a gap were inserted for the seed
    keep = [j for j in range(width) if any(r[j] not in "-." for r in old)]
    # an all-gap column already present in the plate is dropped by `keep` too; restore by comparison
    # against the original with its own all-gap columns removed
    own = [j for j in range(len(seqs[0])) if any(r[j] not in "-." for r in seqs)]
    for orig, new in zip(seqs, old):
        if "".join(new[j] for j in keep) != "".join(orig[j] for j in own):
            return "WARN mafft moved bases or gaps in an original row; left unchanged"
    names_out = [names[0], label] + names[1:]
    seqs_out = [old[0], s2[-1]] + old[1:]
    write_fa(path, names_out, seqs_out)
    return "added %s (columns %d -> %d; all-gap columns dropped, others unchanged)" % (label, len(seqs[0]), width)


def read_fa_text(text):
    names, seqs = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            names.append(line[1:].strip())
            seqs.append([])
        elif seqs:
            seqs[-1].append(line.strip())
    return names, ["".join(s) for s in seqs]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bank")
    ap.add_argument("code")
    ap.add_argument("alignments", nargs="+")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    bn, bs = read_fa(a.bank)
    bank = {n.split()[0]: re.sub(r"[^A-Za-z]", "", s) for n, s in zip(bn, bs)}
    for p in a.alignments:
        if not os.path.isfile(p):  # an unmatched shell glob (no subfam plate, say)
            continue
        print("%-45s %s" % (os.path.basename(p), add_one(p, bank, a.code, a.threads)))


if __name__ == "__main__":
    main()
