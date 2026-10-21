"""add_seed_row.mark keeps flanks packed (2026-10-21).

Where the copies' element reaches past the original consensus, the columns between are aligned element columns. mark()
relabels them flank (lowercase); until 2026-10-21 it changed only the case, so their gaps sat inside the flank and the
plate showed unpacked flanks (rsi_main_20261021 r10_19seqs, r6, r7: every copy).
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import add_seed_row as A  # noqa: E402

EXT = "------gg-gTTTTTTTTTTgg-g------"   # row 1: rebuilt, proposes 3 bp on each side of the original
ORIG = "----------TTTTTTTTTT----------"  # row 2: the consensus as searched
COPIES = [
    "--acgtGG-GTTTTTTTTTTCC-Ctgca--",    # flank packed against an element that reaches past the original
    "---cgtGGAGTTTATTTTTTCCACtg----",
    "-----tG--GTTTTTTTTTTC--Ctgcatt",
]


def test_flanks_packed_after_mark():
    out, stats = A.mark([EXT, ORIG] + COPIES)
    up = [j for j, c in enumerate(out[0]) if c.isupper()]
    lo, hi = up[0], up[-1]
    for before, after in zip(COPIES, out[2:]):
        assert len(after) == len(before)
        assert A.ungap(after).upper() == A.ungap(before).upper()        # no letter lost or moved between rows
        assert not re.search(r"[A-Za-z]-+[A-Za-z]", after[:lo]), after  # 5' flank: gaps, then letters to the edge
        assert not re.search(r"[A-Za-z]-+[A-Za-z]", after[hi + 1:]), after
        assert after[:lo] == after[:lo].lower() and after[hi + 1:] == after[hi + 1:].lower()
        assert after[lo:hi + 1] == after[lo:hi + 1].upper()


def test_mark_is_idempotent():
    once, _ = A.mark([EXT, ORIG] + COPIES)
    twice, _ = A.mark(once)
    assert once[2:] == twice[2:]
