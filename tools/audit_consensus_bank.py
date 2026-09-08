#!/usr/bin/env python3
"""Audit consensus bank for RC duplicates, N content, step4 divergence."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from consensus_bank_lib import (
    find_rc_clusters,
    identity,
    orient_by_simple_repeat_tail,
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
    ap.add_argument("--orient-check", action="store_true",
                    help="Report tail-orient status; exit 1 if any would flip")
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
        o = orient_by_simple_repeat_tail(u)
        print(f"  {f}: len={len(u)} N={u.count('N')} orient={o.action} "
              f"fwd={o.fwd_score:.1f} rev={o.rev_score:.1f}")
        if st:
            print(f"    divergence median={st['median']:.1f}% n={st['n']}")

    if args.orient_check:
        bad = []
        for name, s in cons.items():
            o = orient_by_simple_repeat_tail(s)
            if o.action == "flipped":
                bad.append(name)
                print(f"ORIENT FAIL {name}: would flip (fwd={o.fwd_score:.1f} "
                      f"rev={o.rev_score:.1f})")
        if bad:
            print(f"orient-check: {len(bad)} consensus(es) backwards — fix before step1 sear")
            return 1
        print("orient-check: OK")

    rebuilt = rebuild_from_subfams(run_root)
    for f in args.focus:
        if f in cons and f in rebuilt:
            print(f"  {f} seed vs rebuilt: direct={identity(cons[f], rebuilt[f]):.1f}% "
                  f"RC={identity(cons[f], rc(rebuilt[f])):.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
