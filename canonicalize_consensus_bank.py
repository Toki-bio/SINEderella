#!/usr/bin/env python3
"""Merge RC-duplicate consensus seeds and canonicalize orientation.

See consensus_bank_lib.py and PUBLISH_WORKFLOW.md § Consensus bank.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from consensus_bank_lib import (
    find_rc_clusters,
    identity,
    orient_at_rich_3prime,
    pick_canonical,
    rc,
    read_fa,
    ungap,
    write_fa,
)


def merge_cluster(names: list[str], cons: dict[str, str]) -> tuple[str, str, list[tuple]]:
    kept = sorted(names)[0]
    seq = ungap(cons[kept])
    aliases = []
    for other in names:
        if other == kept:
            continue
        canon, _, reason = pick_canonical(seq, cons[other])
        id_d = identity(seq, cons[other])
        id_r = identity(seq, rc(cons[other]))
        seq = canon
        aliases.append((other, kept, reason, max(id_d, id_r)))
    seq = orient_at_rich_3prime(seq)
    note = "simple-repeat tail at 3 prime"
    aliases = [(a, b, f"{r}; {note}", p) for a, b, r, p in aliases]
    return kept, seq, aliases


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_fa", type=Path)
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument("--aliases", type=Path, default=None)
    ap.add_argument("--min-id", type=float, default=80.0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cons = read_fa(args.input_fa)
    if not cons:
        print("ERROR: empty input", file=sys.stderr)
        return 1

    clusters = find_rc_clusters(cons, args.min_id)
    drop: set[str] = set()
    alias_rows: list[tuple] = []
    out = dict(cons)

    for cl in clusters:
        kept, seq, aliases = merge_cluster(cl, cons)
        out[kept] = seq
        for dropped, k, reason, pct in aliases:
            drop.add(dropped)
            alias_rows.append((dropped, k, reason, f"{pct:.1f}"))
            out.pop(dropped, None)
        print(f"MERGE: keep {kept} <- {[a[0] for a in aliases]}")

    for name in list(out):
        out[name] = orient_at_rich_3prime(out[name])

    print(f"{len(cons)} -> {len(out)} ({len(drop)} merged)")
    if args.dry_run:
        return 0

    write_fa(args.output, out)
    if args.aliases:
        with open(args.aliases, "w") as fh:
            fh.write("dropped\tkept\treason\tidentity_pct\n")
            for row in sorted(alias_rows):
                fh.write("\t".join(row) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
