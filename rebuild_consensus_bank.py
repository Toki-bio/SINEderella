#!/usr/bin/env python3
"""Rebuild consensus bank from step2 assigned copies (majority, no N ties).

Writes consensuses.rebuilt.fa for step4 plots. Does not mutate consensuses.clean.fa
unless --inplace is passed (backs up to .pre_rebuild.bak).
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from consensus_bank_lib import orient_at_rich_3prime, read_fa, rebuild_from_subfams, write_fa


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_root", type=Path)
    ap.add_argument("-o", "--output", type=Path, default=None,
                    help="Default: RUN_ROOT/consensuses.rebuilt.fa")
    ap.add_argument("--inplace", action="store_true",
                    help="Replace consensuses.clean.fa (backup .pre_rebuild.bak)")
    ap.add_argument("--max-seqs", type=int, default=500)
    args = ap.parse_args()

    run_root = args.run_root.resolve()
    rebuilt = rebuild_from_subfams(run_root, max_seqs=args.max_seqs)
    if not rebuilt:
        print("ERROR: no subfamily FASTAs under step2_output/subfamilies/", file=sys.stderr)
        return 1

    for name in list(rebuilt):
        rebuilt[name] = orient_at_rich_3prime(rebuilt[name])

    out = args.output or (run_root / "consensuses.rebuilt.fa")
    write_fa(out, rebuilt)
    print(f"Wrote {len(rebuilt)} consensuses -> {out}")

    clean = run_root / "consensuses.clean.fa"
    if args.inplace and clean.is_file():
        bak = run_root / "consensuses.clean.fa.pre_rebuild.bak"
        shutil.copy2(clean, bak)
        write_fa(clean, rebuilt)
        print(f"Inplaced consensuses.clean.fa (backup {bak})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
