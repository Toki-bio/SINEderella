#!/usr/bin/env python3
"""Put the consensus exactly as searched back on each published plate, and mark every automatic
change to it as a proposal.

Usage: add_seed_row.py BANK.fa SPECIES_CODE ALIGNMENT.aln.fa [...] [--threads N]

The publish chain rebuilds row 1 of every plate from the copies: it can add bases past the ends
of the consensus the genome was searched with, and it can drop bases from them. Neither is safe to
apply automatically (his review, 2026-09-28: on rsi MEG-RL the column walk put TGGGGGAAATA of
flank in front of the real start, and there is no rule for trimming that he trusts). So the plate
states the proposal and leaves the decision to him:

  row 1  <subfamily>_extended   the consensus rebuilt from the copies. LOWERCASE = a proposed
                                change, never applied - compare with row 2 directly below:
                                  outside the original's span: a proposed addition;
                                  at an end inside it: a proposed trim - an original base the
                                  copies do not carry, put back in lowercase.
                                UPPERCASE = what the copies support within the original: the
                                element the verdict judges.
  row 2  <subfamily>            the consensus exactly as searched, under its own name.
  copies                        uppercase over row 1's uppercase span, lowercase outside. Only
                                the case changes; no row is realigned or repacked.

The original is added with `mafft --add` into the finished alignment: existing rows are not
realigned, only gap columns are inserted where it needs them. The file is rewritten only if every
existing row is unchanged once those inserted all-gap columns are removed; otherwise it is left as
it was and a warning is printed. Plates published before this (row 2 named
<subfamily>_seed_as_searched, row 1 named <subfamily>) are converted without realigning.
If the plate is in the reverse-complement orientation of the original (it should not be after
correct_published_aln.py), the reverse complement is added and its name says so.

Each plate's proposals go to proposals.tsv next to it (one row per plate, replaced on re-run):
bp added / trimmed at each end, and how well the copies support them - the median per-copy
identity to the proposed bases, measured on each copy's own ungapped sequence next to the
original's edge for additions (unrelated DNA gives ~0.25) and on the alignment columns for trims.
"""
import argparse
import csv
import os
import random
import re
import statistics
import subprocess
import tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
KINDS = ("_top100.aln.fa", "_rand100.aln.fa", "_subfam.aln.fa")
OLD_TAG = "_seed_as_searched"
EXT = "_extended"
GAPS = "-."
STRAY_GAP = 5       # columns without an original letter that cut an end block off the main one
STRAY_MAX = 12      # an end block this small is a stray, not the original's span,
STRAY_OCC = 0.50    # ... unless at least this share of the copies have letters in its columns
STRAY_ID = 0.40     # ... and match its letters at this median identity (flank background ~0.25)
FIELDS = ["plate", "subfamily", "orig_len",
          "add5_bp", "add5_ungapped", "add5_support", "add5_null",
          "add3_bp", "add3_ungapped", "add3_support", "add3_null",
          "trim5_bp", "trim5_support", "trim3_bp", "trim3_support", "copies", "rebuilt_vs_orig"]
TRIM_FIELDS = ("trim5_bp", "trim5_support", "trim3_bp", "trim3_support")


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


def read_fa_text(text):
    names, seqs = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            names.append(line[1:].strip())
            seqs.append([])
        elif seqs:
            seqs[-1].append(line.strip())
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


def add_original(names, seqs, seed, sf, threads):
    """mafft --add the original; returns (names, seqs) with it as row 2, or a warning string."""
    body = ungap(seqs[0]).upper()
    rc = seed.translate(COMP)[::-1]
    label = sf
    if kshare(rc, body) > 2 * kshare(seed, body) and kshare(rc, body) >= 0.10:
        seed, label = rc, sf + "_revcomp"
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as td:
        aln, add = os.path.join(td, "aln.fa"), os.path.join(td, "seed.fa")
        write_fa(aln, ["r%d" % i for i in range(len(seqs))], seqs)   # mafft wants simple unique names
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
    keep = [j for j in range(width) if any(r[j] not in GAPS for r in old)]
    own = [j for j in range(len(seqs[0])) if any(r[j] not in GAPS for r in seqs)]
    for orig, new in zip(seqs, old):
        if "".join(new[j] for j in keep) != "".join(orig[j] for j in own):
            return "WARN mafft moved bases or gaps in an original row; left unchanged"
    return [names[0], label] + names[1:], [old[0], s2[-1]] + old[1:]


def anchored_identity(p, f):
    """Identity of proposal p to flank f, both read outward from the element edge: p is aligned in
    full from f's first base, gaps allowed (so an indel in a copy does not wreck the rest), and may
    end anywhere in f. (len(p) - edits) / len(p)."""
    if not p:
        return 0.0
    prev = [j for j in range(len(f) + 1)]
    for i in range(1, len(p) + 1):
        cur = [i] + [0] * len(f)
        pc = p[i - 1]
        for j in range(1, len(f) + 1):
            cur[j] = min(prev[j - 1] + (pc != f[j - 1]), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return max(0.0, (len(p) - min(prev)) / float(len(p)))


def _shuffled(p, k):
    x = list(p)
    random.Random(1000 + k).shuffle(x)
    return "".join(x)


def median_or_blank(xs):
    return "%.2f" % statistics.median(xs) if xs else ""


def mark(seqs):
    """Case marking + trim restoration on a plate whose row 2 is the original. Returns (seqs, stats)."""
    ext, orig = list(seqs[0]), seqs[1]
    oc = [j for j, c in enumerate(orig) if c not in GAPS]
    if not oc:
        return seqs, None
    n_orig = len(oc)
    # The original's span is its MAIN block. MAFFT sometimes aligns a few of its end bases on their
    # own, far from the rest (SINEbase Rhin-1's first GGGG sat ~40 columns left of the body on cth,
    # cse and fho; VES's first 5 bases on cth). First-to-last letter then stretched the span and made
    # junk rebuilt letters in between uppercase. Stray end blocks (<= STRAY_MAX letters, > STRAY_GAP
    # columns from the next original letter) are shown as trim proposals instead.
    # A block is stray only where the copies are NOT: cth Rhin-1's original head GGGGCGGCCGGT sits 25
    # columns from its body, but so do the copies' heads (occupancy 1.0, support 9) - treating it as
    # stray made every copy's head "shared flank" (85 of 100) and the call went SINE 96.8 -> Not SINE.
    ncop = max(1, len(seqs) - 2)

    # Occupancy alone is not enough on a published plate: flanks are packed against the element, so
    # flank letters fill a stray block's columns (mtu Rhin-1 GGGG). The copies must also MATCH it:
    # cth head 0.58, rsi r4 0.83; stray GGGG on rna/vmu/tni/cse/fho 0.00-0.25. The gap was 15; mtu's
    # GGGG sat 14 columns from the body and was never split off.
    def occupied(cols):
        if sum(sum(1 for i in range(2, len(seqs)) if seqs[i][j] not in GAPS) for j in cols) \
                / float(ncop * len(cols)) < STRAY_OCC:
            return False
        ids = []
        for i in range(2, len(seqs)):
            pr = [(seqs[i][j].upper(), orig[j].upper()) for j in cols if seqs[i][j] not in GAPS]
            if pr and len(pr) >= len(cols) / 2.0:
                ids.append(sum(a == b for a, b in pr) / float(len(pr)))
        return bool(ids) and sorted(ids)[len(ids) // 2] >= STRAY_ID

    stray = {"5": [], "3": []}
    while True:
        runs, cur = [], [oc[0]]
        for j in oc[1:]:
            if j - cur[-1] > STRAY_GAP:
                runs.append(cur)
                cur = [j]
            else:
                cur.append(j)
        runs.append(cur)
        if len(runs) > 1 and len(runs[0]) <= STRAY_MAX and not occupied(runs[0]):
            stray["5"] += runs[0]
            oc = [j for r in runs[1:] for j in r]
            continue
        if len(runs) > 1 and len(runs[-1]) <= STRAY_MAX and not occupied(runs[-1]):
            stray["3"] += runs[-1]
            oc = [j for r in runs[:-1] for j in r]
            continue
        break
    s_lo, s_hi = oc[0], oc[-1]
    copies = list(range(2, len(seqs)))
    stats = {"orig_len": n_orig, "copies": len(copies)}

    # proposed trims: original bases the rebuild dropped at an end - put back. Only where row 1 does
    # not extend past the original on that side; otherwise the columns are ambiguous and nothing is
    # restored (the dropped bases stay visible in row 2 directly below).
    letters = [j for j, c in enumerate(ext) if c not in GAPS]
    restored_all = []
    for side in ("5", "3"):
        # stray original bases on this side: a trim proposal wherever row 1 has no letter there
        restored = []
        for j in stray[side]:
            if ext[j] in GAPS:
                ext[j] = orig[j].lower()
            restored.append(j)
        if letters and side == "5" and letters[0] >= s_lo:
            j = s_lo
            while j <= s_hi and ext[j] in GAPS:
                if orig[j] not in GAPS:
                    ext[j] = orig[j].lower()
                    restored.append(j)
                j += 1
        elif letters and side == "3" and letters[-1] <= s_hi:
            j = s_hi
            while j >= s_lo and ext[j] in GAPS:
                if orig[j] not in GAPS:
                    ext[j] = orig[j].lower()
                    restored.append(j)
                j -= 1
        sup = []
        for i in copies:
            r = seqs[i]
            pairs = [(r[j].upper(), orig[j].upper()) for j in restored if r[j] not in GAPS]
            if pairs and len(pairs) >= len(restored) / 2.0:
                sup.append(sum(a == b for a, b in pairs) / float(len(pairs)))
        restored_all += restored
        stats["trim%s_bp" % side] = len(restored)
        stats["trim%s_support" % side] = median_or_blank(sup) if restored else ""

    # proposed additions: row-1 letters outside the original's span, lowercase (not the stray trims)
    trims = set(restored_all)
    add5 = [j for j, c in enumerate(ext) if c not in GAPS and j < s_lo and j not in trims]
    add3 = [j for j, c in enumerate(ext) if c not in GAPS and j > s_hi and j not in trims]
    for j in range(len(ext)):
        if ext[j] not in GAPS:
            ext[j] = ext[j].lower() if (j < s_lo or j > s_hi or j in trims) else ext[j].upper()
    # support: each copy's own ungapped flank next to the original's edge, read outward, aligned to
    # the proposal with gaps allowed; null = the same with the proposal shuffled (gapped alignment
    # of unrelated DNA scores well above 0.25, so the null is what support is read against)
    for side, cols in (("5", add5), ("3", add3)):
        prop = "".join(ext[j] for j in cols).upper()
        sup, null, ung = [], [], []
        for k, i in enumerate(copies if prop else ()):
            r = seqs[i]
            if side == "5":
                fl = ungap(r[:s_lo]).upper()[::-1][:len(prop) + 10]
                pp = prop[::-1]
            else:
                fl = ungap(r[s_hi + 1:]).upper()[:len(prop) + 10]
                pp = prop
            if len(fl) >= len(pp) / 2.0:
                a = fl[:len(pp)]
                ung.append(sum(x == y for x, y in zip(a, pp)) / float(len(a)))   # position by position
                sup.append(anchored_identity(pp, fl))
                null.append(anchored_identity(_shuffled(pp, k), fl))
        stats["add%s_bp" % side] = len(cols)
        stats["add%s_ungapped" % side] = median_or_blank(ung) if prop else ""
        stats["add%s_support" % side] = median_or_blank(sup) if prop else ""
        stats["add%s_null" % side] = median_or_blank(null) if prop else ""

    # identity of the rebuilt consensus to the original where both have a base: low means the copies
    # are a different family than the query found them with (lly "Rhin-1" copies start
    # GGGTTCCCTGGTGGTGTAGTGGC against SINEbase Rhin-1 GGGGGGCCGGTTGCTCAGTTGGT)
    both = [j for j in range(s_lo, s_hi + 1) if ext[j] not in GAPS and orig[j] not in GAPS]
    stats["rebuilt_vs_orig"] = ("%.2f" % (sum(ext[j].upper() == orig[j].upper() for j in both) / float(len(both)))
                                if both else "")
    up = [j for j, c in enumerate(ext) if c.isupper()]
    e_lo, e_hi = (up[0], up[-1]) if up else (s_lo, s_hi)
    out = ["".join(ext), orig.upper()]
    for i in copies:
        out.append("".join(c if c in GAPS else (c.upper() if e_lo <= j <= e_hi else c.lower())
                           for j, c in enumerate(seqs[i])))
    return out, stats


def process(path, bank, code, threads):
    sf = subfamily_of(path, code)
    if sf is None or sf not in bank:
        return "skip (no original %r in bank)" % sf, None
    names, seqs = read_fa(path)
    if len(seqs) < 2:
        return "skip (fewer than 2 rows)", None
    if len(names) > 1 and names[1] in (sf, sf + "_revcomp") and names[0].split()[0] == sf + EXT:
        msg = "re-marked"
    elif len(names) > 1 and names[1].startswith(sf + OLD_TAG):
        names = [sf + EXT, names[1].replace(sf + OLD_TAG, sf, 1)] + names[2:]
        msg = "converted from the %s naming" % OLD_TAG
    else:
        r = add_original(names, seqs, bank[sf].upper(), sf, threads)
        if isinstance(r, str):
            return r, None
        names, seqs = r
        names[0] = sf + EXT
        msg = "added the original as row 2"
    seqs, stats = mark(seqs)
    write_fa(path, names, seqs)
    if stats:
        stats.update(plate=os.path.basename(path), subfamily=sf, remarked=msg == "re-marked")
        msg += "; proposals +%s/+%s bp, trim %s/%s bp (5'/3')" % (
            stats["add5_bp"], stats["add3_bp"], stats["trim5_bp"], stats["trim3_bp"])
    return msg, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bank")
    ap.add_argument("code")
    ap.add_argument("alignments", nargs="+")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    bn, bs = read_fa(a.bank)
    bank = {n.split()[0]: re.sub(r"[^A-Za-z]", "", s) for n, s in zip(bn, bs)}
    by_dir = {}
    for p in a.alignments:
        if not os.path.isfile(p):  # an unmatched shell glob (no subfam plate, say)
            continue
        msg, stats = process(p, bank, a.code, a.threads)
        print("%-45s %s" % (os.path.basename(p), msg))
        if stats:
            by_dir.setdefault(os.path.dirname(os.path.abspath(p)), []).append(stats)
    for d, rows in by_dir.items():
        tsv = os.path.join(d, "proposals.tsv")
        old = {}
        if os.path.isfile(tsv):
            with open(tsv, newline="") as fh:
                old = {r["plate"]: r for r in csv.DictReader(fh, delimiter="\t")}
        for r in rows:
            prev = old.get(r["plate"])
            if prev and r.get("remarked"):
                # a re-run cannot see trims that the first run already put back: keep what it recorded
                for f in TRIM_FIELDS:
                    r[f] = prev.get(f, r.get(f, ""))
            old[r["plate"]] = r
        with open(tsv + ".tmp", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            for k in sorted(old):
                w.writerow(old[k])
        os.replace(tsv + ".tmp", tsv)


if __name__ == "__main__":
    main()
