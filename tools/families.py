#!/usr/bin/env python3
"""Family table for the families+subfamilies mode (MANUAL 6.1.6, docs/FAMILIES.md).

A bank (consensuses.clean.fa) is a flat list of names. This reads the optional family table that says which subfamily belongs to which
family, checks it against the bank, and gives every family a PARENT consensus:

    families.tsv   (tab separated, '#' comments, header line optional)
        subfamily   family   [parent_consensus]

    subfamily          name exactly as in the bank (the bank sanitiser keeps the first word of the header)
    family             family name (free text, no spaces)
    parent_consensus   optional: name of a sequence in --parents (or, failing that, in the bank) used as the parent consensus of the family;
                       give it on any row of the family; two different values for one family are an error

Rules
    * a bank consensus that is not in the table is its own single-member family (named like itself); no table at all = every consensus is
      its own family, which is today's flat behaviour;
    * a name in the table that is not in the bank is an error (a typo would silently drop a subfamily);
    * a subfamily may appear once; a family name may equal the name of its own single subfamily but not the name of a subfamily of another
      family;
    * parent consensus: SUPPLIED when the table names one; otherwise a single-member family uses the member itself, and a family of several
      members gets the DERIVED consensus of its members' consensuses, weighted by copy number when --weights is given (subfamily<TAB>n,
      e.g. the firm counts of an earlier run), else equally. Columns where the weighted share of non-gap bases is below --min-occupancy
      (default 0.5) are dropped.

Outputs (in --out): families.clean.tsv (subfamily, family, parent_source = supplied | derived | member), family_consensuses.fa
(one record per family, header = family name), families.summary.tsv (family, n_subfamilies, parent_source, parent_length, members).
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict

BASES = "ACGT"


class FamilyTableError(ValueError):
    pass


def read_fasta(path: str) -> "OrderedDict[str, str]":
    seqs: "OrderedDict[str, str]" = OrderedDict()
    name = None
    for line in open(path):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in seqs:
                raise FamilyTableError(f"{path}: duplicate sequence name {name}")
            seqs[name] = ""
        elif name is not None:
            seqs[name] += line.strip()
    return seqs


def parse_table(path: str) -> "OrderedDict[str, tuple[str, str]]":
    """subfamily -> (family, parent_name or '')"""
    rows: "OrderedDict[str, tuple[str, str]]" = OrderedDict()
    for ln, line in enumerate(open(path), 1):
        line = line.rstrip("\r\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 2:
            raise FamilyTableError(f"{path}:{ln}: need at least subfamily<TAB>family, got {line!r}")
        sub, fam = f[0].strip(), f[1].strip()
        par = f[2].strip() if len(f) > 2 else ""
        if ln == 1 or not rows:
            if sub.lower() == "subfamily" and fam.lower() == "family":
                continue  # header line
        if not sub or not fam or " " in fam or " " in sub:
            raise FamilyTableError(f"{path}:{ln}: empty or space-containing name in {line!r}")
        if sub in rows:
            raise FamilyTableError(f"{path}:{ln}: subfamily {sub} listed twice")
        rows[sub] = (fam, par)
    return rows


def resolve(bank_names, table, parents_names=()):
    """Return OrderedDict subfamily -> family for every bank name, and family -> supplied parent name ('' if none).
    Checks everything described in the module docstring."""
    bank = list(bank_names)
    bankset = set(bank)
    missing = [s for s in table if s not in bankset]
    if missing:
        raise FamilyTableError("in the family table but not in the bank: " + ", ".join(missing))
    sub2fam: "OrderedDict[str, str]" = OrderedDict()
    for s in bank:
        sub2fam[s] = table[s][0] if s in table else s
    fam_members = defaultdict(list)
    for s, fam in sub2fam.items():
        fam_members[fam].append(s)
    for fam, mem in fam_members.items():
        if fam in bankset and fam not in mem:
            raise FamilyTableError(f"family name {fam} is also the name of a subfamily that belongs to family {sub2fam[fam]}")
    supplied = {}
    for s, (fam, par) in table.items():
        if not par:
            continue
        if fam in supplied and supplied[fam] != par:
            raise FamilyTableError(f"family {fam}: two different parent consensuses ({supplied[fam]}, {par})")
        supplied[fam] = par
    for fam, par in supplied.items():
        if par not in parents_names and par not in bankset:
            raise FamilyTableError(f"family {fam}: parent consensus {par} is in neither the parents file nor the bank")
    return sub2fam, fam_members, supplied


def weighted_consensus(aligned, weights, min_occupancy=0.5):
    """aligned: list of equal-length aligned strings; weights: same length (> 0). Per column: if the weighted share of non-gap ACGT is below
    min_occupancy the column is dropped, else the weighted plurality base is written (ties broken by the order ACGT)."""
    if not aligned:
        return ""
    n = len(aligned[0])
    if any(len(a) != n for a in aligned):
        raise FamilyTableError("sequences are not aligned to one length")
    tot = float(sum(weights))
    out = []
    for i in range(n):
        cnt = dict.fromkeys(BASES, 0.0)
        for a, w in zip(aligned, weights):
            c = a[i].upper()
            if c in cnt:
                cnt[c] += w
        occ = sum(cnt.values()) / tot
        if occ < min_occupancy:
            continue
        out.append(max(BASES, key=lambda b: cnt[b]))
    return "".join(out)


def align_mafft(seqs, mafft="mafft"):
    """seqs: OrderedDict name -> sequence; returns aligned strings in the same order."""
    if shutil.which(mafft) is None:
        raise FamilyTableError(f"{mafft} not found; install it or pass --aligned with pre-aligned members")
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "in.fa")
        with open(inp, "w") as fh:
            for k, v in seqs.items():
                fh.write(f">{k}\n{v}\n")
        res = subprocess.run([mafft, "--auto", "--quiet", "--thread", "1", inp], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            raise FamilyTableError("mafft failed: " + res.stderr[-300:])
    out: "OrderedDict[str, str]" = OrderedDict()
    name = None
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif name:
            out[name] += line.strip()
    return [out[k] for k in seqs]


def read_weights(path):
    w = {}
    for line in open(path):
        f = line.rstrip("\r\n").split("\t")
        if len(f) >= 2 and f[0] and not f[0].startswith("#"):
            try:
                w[f[0]] = float(f[1])
            except ValueError:
                continue  # header
    return w


def build(bank_path, table_path, out_dir, parents_path=None, weights_path=None, aligned=False, mafft="mafft", min_occupancy=0.5):
    bank = read_fasta(bank_path)
    table = parse_table(table_path) if table_path else OrderedDict()
    parents = read_fasta(parents_path) if parents_path else OrderedDict()
    sub2fam, fam_members, supplied = resolve(bank.keys(), table, parents.keys())
    weights = read_weights(weights_path) if weights_path else {}
    os.makedirs(out_dir, exist_ok=True)
    rows, fam_cons, summary = [], OrderedDict(), []
    for fam, mem in fam_members.items():
        if fam in supplied:
            par = supplied[fam]
            seq = parents[par] if par in parents else bank[par]
            src = "supplied"
        elif len(mem) == 1:
            seq, src = bank[mem[0]], "member"
        else:
            seqs = OrderedDict((m, bank[m]) for m in mem)
            al = list(seqs.values()) if aligned else align_mafft(seqs, mafft)
            w = [max(weights.get(m, 1.0), 1e-9) if weights else 1.0 for m in mem]
            seq, src = weighted_consensus(al, w, min_occupancy), "derived"
            if not seq:
                raise FamilyTableError(f"family {fam}: derived consensus is empty (members do not align)")
        fam_cons[fam] = seq
        for m in mem:
            rows.append((m, fam, src))
        summary.append((fam, len(mem), src, len(seq), ",".join(mem)))
    with open(os.path.join(out_dir, "families.clean.tsv"), "w") as fh:
        fh.write("subfamily\tfamily\tparent_source\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    with open(os.path.join(out_dir, "family_consensuses.fa"), "w") as fh:
        for fam, seq in fam_cons.items():
            fh.write(f">{fam}\n{seq}\n")
    with open(os.path.join(out_dir, "families.summary.tsv"), "w") as fh:
        fh.write("family\tn_subfamilies\tparent_source\tparent_length\tmembers\n")
        for s in summary:
            fh.write("\t".join(str(x) for x in s) + "\n")
    return rows, fam_cons, summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bank", required=True, help="consensus bank FASTA (subfamily consensuses)")
    ap.add_argument("--table", default=None, help="families.tsv (omit: every consensus is its own family)")
    ap.add_argument("--parents", default=None, help="FASTA of supplied parent consensuses")
    ap.add_argument("--weights", default=None, help="subfamily<TAB>copy number, for the derived parent consensus")
    ap.add_argument("--out", required=True)
    ap.add_argument("--aligned", action="store_true", help="the bank members of each family are already aligned to one length (tests)")
    ap.add_argument("--mafft", default="mafft")
    ap.add_argument("--min-occupancy", type=float, default=0.5)
    a = ap.parse_args(argv)
    try:
        rows, fam_cons, summary = build(a.bank, a.table, a.out, a.parents, a.weights, a.aligned, a.mafft, a.min_occupancy)
    except FamilyTableError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    multi = [s for s in summary if s[1] > 1]
    print(f"{len(rows)} subfamilies, {len(fam_cons)} families ({len(multi)} with several subfamilies); "
          f"parents: " + ", ".join(f"{k}={sum(1 for s in summary if s[2] == k)}" for k in ("supplied", "derived", "member")), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
