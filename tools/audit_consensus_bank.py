#!/usr/bin/env python3
"""Audit consensus bank for RC duplicates, N content, step4 divergence."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from consensus_bank_lib import (
    find_rc_clusters,
    identity,
    pctid_stats,
    pick_canonical,
    rc,
    read_fa,
    rebuild_from_subfams,
    ungap,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_root", type=Path)
    ap.add_argument("--min-id", type=float, default=90.0)
    ap.add_argument("--focus", nargs="*", default=[])
    args = ap.parse_args()
    run_root = args.run_root.resolve()
    cons_path = run_root / "consensuses.clean.fa"
    if not cons_path.is_file():
        print(f"ERROR: missing {cons_path}", file=sys.stderr)
        return 1

    cons = read_fa(cons_path)
    print(f"consensus bank: {len(cons)} sequences")

    clusters = find_rc_clusters(cons, args.min_id)
    print(f"RC/similarity clusters (>={args.min_id}%): {len(clusters)}")
    for cl in clusters:
        print(f"  {cl}")
        rep = cl[0]
        for other in cl[1:]:
            d = identity(cons[rep], cons[other])
            r = identity(cons[rep], rc(cons[other]))
            _, _, reason = pick_canonical(cons[rep], cons[other])
            print(f"    {rep} vs {other}: direct={d:.1f}% RC={r:.1f}% -> {reason}")

    for f in args.focus:
        s = cons.get(f)
        if not s:
            print(f"  {f}: MISSING")
            continue
        u = ungap(s)
        st = pctid_stats(run_root, f)
        print(f"  {f}: len={len(u)} N={u.count('N')}")
        if st:
            print(f"    divergence median={st['median']:.1f}% n={st['n']}")

    rebuilt = rebuild_from_subfams(run_root)
    for f in args.focus:
        if f in cons and f in rebuilt:
            print(f"  {f} seed vs rebuilt: direct={identity(cons[f], rebuilt[f]):.1f}% "
                  f"RC={identity(cons[f], rc(rebuilt[f])):.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
