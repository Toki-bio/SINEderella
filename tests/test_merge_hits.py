"""tools/merge_hits.py: strand-aware merge of the hits of all consensuses (no '+,-' loci, no reverse-complemented chimeras)."""
import io
import os
import subprocess
import sys

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools", "merge_hits.py")


def run(rows, *args):
    text = "".join("\t".join(map(str, r)) + "\n" for r in rows)
    r = subprocess.run([sys.executable, TOOL, *args], input=text, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return [l.split("\t") for l in r.stdout.splitlines()], r.stderr


def test_pure_strand_clusters_as_before():
    out, _ = run([("c1", 0, 100, "80.5", 150, "+"), ("c1", 90, 220, "85.25", 160, "+"), ("c1", 220, 300, "70", 90, "+"), ("c1", 500, 600, "90", 200, "-")])
    assert out == [["c1", "0", "300", "85.25", "160", "+"], ["c1", "500", "600", "90", "200", "-"]]


def test_old_mode_joins_strands():
    out, _ = run([("c1", 0, 100, "80", 150, "+"), ("c1", 50, 200, "90", 160, "-")], "--old")
    assert out == [["c1", "0", "200", "90", "160", "+,-"]]


def test_adjacent_elements_on_opposite_strands_are_split():
    # DIP on "-" (6755806-6756140) next to B1 on "+" (6756090-6756230): the real case ptg000073l, overlap 50 bp of 140
    rows = [("c", 6755806, 6756040, "80", 206, "-"), ("c", 6755823, 6756140, "81", 301, "-"), ("c", 6755835, 6756040, "82", 204, "-"),
            ("c", 6756090, 6756230, "83", 126, "+"), ("c", 6756090, 6756230, "84", 144, "+")]
    out, err = run(rows)
    assert out == [["c", "6755806", "6756140", "82", "301", "-"], ["c", "6756090", "6756230", "84", "144", "+"]]
    assert "1 split into 2 pieces" in err


def test_same_locus_seen_from_both_strands_is_one_interval_with_the_best_strand():
    out, err = run([("c", 100, 300, "88", 400, "+"), ("c", 110, 290, "70", 120, "-")])
    assert out == [["c", "100", "300", "88", "400", "+"]]
    assert "1 merged as one locus" in err
    out, _ = run([("c", 100, 300, "75", 90, "+"), ("c", 110, 290, "90", 300, "-")])
    assert out == [["c", "100", "300", "90", "300", "-"]]


def test_book_ended_opposite_strands_are_two_loci():
    out, _ = run([("c", 0, 200, "80", 100, "+"), ("c", 200, 400, "80", 120, "-")])
    assert [(o[1], o[2], o[5]) for o in out] == [("0", "200", "+"), ("200", "400", "-")]


def test_three_elements_plus_minus_plus_and_sorting():
    rows = [("c", 0, 150, "80", 100, "+"), ("c", 140, 300, "80", 110, "-"), ("c", 290, 450, "80", 120, "+"), ("d", 5, 50, "70", 60, "-")]
    out, _ = run(rows)
    assert [(o[0], o[1], o[2], o[5]) for o in out] == [("c", "0", "150", "+"), ("c", "140", "300", "-"), ("c", "290", "450", "+"), ("d", "5", "50", "-")]


def test_no_plus_minus_strand_is_ever_written():
    import random
    r = random.Random(3)
    rows = []
    for _ in range(300):
        s = r.randint(0, 5000)
        rows.append(("c", s, s + r.randint(30, 300), "80", r.randint(50, 300), r.choice("+-")))
    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    out, _ = run(rows)
    assert all(o[5] in ("+", "-") for o in out)
    starts = [int(o[1]) for o in out]
    assert starts == sorted(starts)
    # no two pieces of the same strand overlap
    for st in "+-":
        iv = sorted((int(o[1]), int(o[2])) for o in out if o[5] == st)
        assert all(a[1] <= b[0] for a, b in zip(iv, iv[1:]))
