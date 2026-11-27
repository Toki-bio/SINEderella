"""tools/composite_pieces.sh: composite stretches as step-1 queries whose header names the composite."""
import os
import subprocess
import sys

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools", "composite_pieces.sh")


def test_pieces(tmp_path):
    bank = tmp_path / "bank.fa"
    bank.write_text(">U1\n" + "A" * 80 + "\n>C2\n" + "G" * 60 + "\n>C1\n" + "C" * 60 + "G" * 40 + "T" * 30 + "\n")
    units = tmp_path / "units.tsv"
    units.write_text("consensus\tkind\tstart\tend\tunit\tidentity\talternatives\n"
                     "U1\tunit\t1\t80\tU1\t100.0\t-\n"
                     "C1\tcomposite\t1\t60\tU1\t91.0\t-\n"       # 60 bp: a piece
                     "C1\tcomposite\t61\t100\t-\t-\t-\n"         # 40 bp: too short
                     "C1\tcomposite\t101\t130\t-\t-\t-\n"         # 30 bp: too short
                     "C2\tcomposite\t1\t60\t-\t-\t-\n")           # 60 bp, but no unit explains it: never a piece
    out = tmp_path / "s"
    r = subprocess.run(["bash", TOOL, str(units), str(bank), str(out)], capture_output=True, text=True, check=True)
    assert r.stdout.strip() == "1"
    f = out / "pieces.part_C1__piece_1-60.fa"
    assert f.read_text() == ">C1 piece_1-60 unit=U1\n" + "C" * 60 + "\n"      # first word = the composite: hits get its label
    assert (out / "pieces.named.fa").read_text() == ">C1__piece_1-60\n" + "C" * 60 + "\n"
    r = subprocess.run(["bash", TOOL, str(units), str(bank), str(out), "30"], capture_output=True, text=True, check=True)
    assert r.stdout.strip() == "1"                                          # MIN_LEN 30: the "-" stretches still are not pieces
    assert sorted(os.listdir(out)) == ["pieces.named.fa", "pieces.part_C1__piece_1-60.fa"]
