#!/usr/bin/env python3
"""Flip published gapped MSAs when row-1 consensus element is backwards.

Secondary guard after rebuild_consensus_row (MAFFT --adjustdirection drift).
Primary orientation must be set on the consensus bank before step1 sear.
"""
from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from consensus_bank_lib import orient_by_simple_repeat_tail, ungap

COMP = str.maketrans(
    "ACGTacgt",
    "TGCAtgca",
)


def read_aln(path: Path) -> tuple[list[str], list[str]]:
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:].split()[0]
            names.append(cur)
            buf = []
        else:
            buf.append(line)
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


def write_aln(path: Path, names: list[str], seqs: list[str]) -> None:
    with io.open(path, "w", encoding="utf-8") as fh:
        for n, s in zip(names, seqs):
            fh.write(f">{n}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i : i + 80] + "\n")


def rc_gapped_row(s: str) -> str:
    return s.translate(COMP)[::-1]


def element_from_row(seq: str) -> str:
    return ungap("".join(c for c in seq if c.isupper()))


def apply(path: Path, dry_run: bool) -> str | None:
    names, seqs = read_aln(path)
    if len(seqs) < 2:
        return None
    elem = element_from_row(seqs[0])
    if len(elem) < 30:
        print(f"SKIP {path.name}: consensus element too short ({len(elem)} bp)")
        return "broken"
    r = orient_by_simple_repeat_tail(elem)
    if r.action != "flipped":
        return r.action
    if dry_run:
        print(f"WOULD FLIP {path.name}  fwd={r.fwd_score:.1f} rev={r.rev_score:.1f}")
        return "would_flip"
    seqs = [rc_gapped_row(s) for s in seqs]
    write_aln(path, names, seqs)
    print(f"FLIPPED {path.name}  fwd={r.fwd_score:.1f} rev={r.rev_score:.1f}")
    return "flipped"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("paths", nargs="+", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    counts: dict[str, int] = {}
    for p in args.paths:
        if not p.is_file():
            print(f"MISSING {p}", file=sys.stderr)
            continue
        act = apply(p, args.dry_run)
        if act:
            counts[act] = counts.get(act, 0) + 1
    if counts:
        print("summary:", " ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
