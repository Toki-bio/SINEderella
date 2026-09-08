#!/usr/bin/env python3
"""Orient every consensus in a bank by simple tandem repeat tail at 3'.

Run after RC merge, before step1 sear. Idempotent when already correct.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from consensus_bank_lib import orient_by_simple_repeat_tail, read_fa, write_fa


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_fa", type=Path)
    ap.add_argument("-o", "--output", type=Path, default=None,
                    help="Default: overwrite input")
    ap.add_argument("--report", type=Path, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cons = read_fa(args.input_fa)
    if not cons:
        print("ERROR: empty input", file=sys.stderr)
        return 1

    out = {}
    rows = []
    n_flip = n_undec = 0
    for name in sorted(cons):
        r = orient_by_simple_repeat_tail(cons[name])
        out[name] = r.seq if r.seq else cons[name]
        rows.append((name, r.action, f"{r.fwd_score:.1f}", f"{r.rev_score:.1f}"))
        if r.action == "flipped":
            n_flip += 1
            print(f"FLIP {name}  fwd={r.fwd_score:.1f} rev={r.rev_score:.1f}")
        elif r.action == "undecided":
            n_undec += 1
            print(f"UNDECIDED {name}  fwd={r.fwd_score:.1f} rev={r.rev_score:.1f}")

    print(f"orient: {len(cons)} seqs, {n_flip} flipped, {n_undec} undecided")
    if args.dry_run:
        return 0

    dest = args.output or args.input_fa
    write_fa(dest, out)
    if args.report:
        with open(args.report, "w") as fh:
            fh.write("name\taction\tfwd_score\trev_score\n")
            for row in rows:
                fh.write("\t".join(row) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
