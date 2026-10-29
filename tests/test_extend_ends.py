"""tools/extend_ends.py: loci whose element goes on past the hits are extended to the element's end, nothing else moves.
Needs bedtools on PATH."""
import os
import random
import subprocess
import sys

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools", "extend_ends.py")


def mutate(s, rate, rnd):
    return "".join(rnd.choice("ACGT".replace(c, "")) if rnd.random() < rate else c for c in s)


def test_extend_ends(tmp_path):
    rnd = random.Random(7)
    body = "".join(rnd.choice("ACGT") for _ in range(150))
    tail = "".join(rnd.choice("ACGT") for _ in range(40))       # the part the bank consensus E lacks
    other = "".join(rnd.choice("ACGT") for _ in range(190))     # F: a complete consensus
    seq, merged, labeled, truth = [], [], [], []
    pos = 0

    def add(s):
        nonlocal pos
        seq.append(s); pos += len(s)

    add("".join(rnd.choice("ACGT") for _ in range(500)))
    for i in range(300):
        kind = "cut" if i < 100 else ("full" if i < 200 else "other")
        el = mutate(body + tail if kind != "other" else other, 0.08, rnd) + "A" * 15
        minus = i % 2 == 1
        start = pos
        add(el if not minus else el[::-1].translate(str.maketrans("ACGT", "TGCA")))
        hit_len = 150 if kind == "cut" else len(el) - 15            # hits stop before the A tail
        if not minus:
            s, e = start, start + hit_len
        else:
            s, e = start + len(el) - hit_len, start + len(el)
        st = "-" if minus else "+"
        merged.append("chr1\t%d\t%d\t80\t100\t%s" % (s, e, st))
        labeled.append("chr1\t%d\t%d\t%s\t100\t%s" % (s, e, "F" if kind == "other" else "E", st))
        truth.append((kind, len(el) - hit_len))
        add("".join(rnd.choice("ACGT") for _ in range(400)))        # unique flank
    (tmp_path / "g.fa").write_text(">chr1\n" + "".join(seq) + "\n")
    (tmp_path / "g.fa.fai").write_text("chr1\t%d\t6\t%d\t%d\n" % (pos, pos, pos + 1))
    (tmp_path / "m.bed").write_text("\n".join(merged) + "\n")
    (tmp_path / "l.bed").write_text("\n".join(labeled) + "\n")
    r = subprocess.run([sys.executable, TOOL, str(tmp_path / "m.bed"), str(tmp_path / "l.bed"), str(tmp_path / "g.fa"),
                        str(tmp_path / "o.bed"), str(tmp_path / "x")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    out = [l.split("\t") for l in (tmp_path / "o.bed").read_text().splitlines()]
    assert len(out) == 300 and all(len(f) == 6 for f in out)
    good = 0
    for (kind, need), old, new in zip(truth, merged, out):
        o = old.split("\t")
        added = (int(new[2]) - int(o[2])) if o[5] == "+" else (int(o[1]) - int(new[1]))
        moved5 = (int(o[1]) - int(new[1])) if o[5] == "+" else (int(new[2]) - int(o[2]))
        assert moved5 == 0
        if kind == "cut":
            assert added <= need + 5                                 # never into the flank (a flank that starts with A's
            #                                                          is read as more A tail: a few bp)
            good += added >= need - 12
        else:
            assert added == 0, (kind, added)                         # complete copies do not move
    assert good >= 75                                                # most cut copies regain tail + A tail (8 % divergence)
    tsv = (tmp_path / "x" / "ends_extension.tsv").read_text().splitlines()
    assert any(l.startswith("E\t3'") and l.endswith("flag") for l in tsv)
    assert not any(l.startswith("F\t") and l.endswith("flag") for l in tsv)
    assert "EXTEND loci=300" in r.stdout
