#!/usr/bin/env python3
"""Extract total copy counts per subfamily from oma_report composition table."""
import re
import sys
from pathlib import Path


def from_report(html_text: str) -> dict:
    m = re.search(r"id=['\"]composition['\"].*?</section>", html_text, re.DOTALL)
    if not m:
        return {}
    counts = {}
    pat = re.compile(
        r"<tr><td>(oma_[^<]+)</td>"
        r"<td>\d+</td><td>\d+</td>"
        r'<td class="hl">(\d+)</td>'
    )
    for row in pat.finditer(m.group(0)):
        counts[row.group(1)] = int(row.group(2))
    return counts


def from_tsv(path: Path) -> dict:
    counts = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        sf = parts[0].strip()
        if not sf.startswith("oma_"):
            continue
        # firm + soft = total (cols 1+2) or col 3 if present
        try:
            if parts[1].strip().isdigit() and parts[2].strip().isdigit():
                counts[sf] = int(parts[1]) + int(parts[2])
            else:
                counts[sf] = int(parts[3])
        except (ValueError, IndexError):
            pass
    return counts


def main():
    p = Path(sys.argv[1])
    text = p.read_text(encoding="utf-8")
    c = from_report(text)
    if not c and p.suffix == ".tsv":
        c = from_tsv(p)
    for sf in sorted(c):
        print("%s\t%d" % (sf, c[sf]))


if __name__ == "__main__":
    main()
