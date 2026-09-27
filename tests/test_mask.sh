#!/usr/bin/env bash
# Masked search, end to end (docs/MASKING.md). Needs the SINEderella toolchain on PATH.
#
# A synthetic genome carries planted copies of two unrelated families: A (masked) and B (not).
# 1. mask_genome.py refuses a BED naming an unknown sequence and one running past a sequence end,
#    and keeps headers and lengths.
# 2. An unmasked run finds A and B.
# 3. A run with --mask-bed on A's intervals finds no A copy and every B copy at its planted
#    coordinates, and records the mask in manifest.txt.
# Sequence names contain '_' on purpose: SINEderella rewrites '_' as '@U@' in genome headers, and a
# BED written against the original assembly must still match.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SD="$(dirname "$HERE")"
W="${1:-$(mktemp -d "${TMPDIR:-$HOME/tmp}/test_mask.XXXXXX")}"
mkdir -p "$W"; cd "$W"
export THREADS="${THREADS:-4}" SKIP_ASSEMBLY_QC=1 SKIP_CANONICALIZE=1
fail(){ echo "FAIL: $*" >&2; exit 1; }

python3 - <<'EOF'
import random
random.seed(7)
B = "ACGT"
def rnd(n): return "".join(random.choice(B) for _ in range(n))
def mut(s, rate=0.08):
    return "".join(random.choice(B.replace(c, "")) if random.random() < rate else c for c in s)
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
famA, famB = rnd(220) + "A" * 12, rnd(220) + "A" * 12
seqs = {"chr_one": list(rnd(400000)), "chr_two": list(rnd(300000))}
truth = []
used = {k: [] for k in seqs}
for fam, cons, n in (("A", famA, 80), ("B", famB, 60)):
    for i in range(n):
        c = random.choice(list(seqs))
        while True:
            p = random.randrange(1000, len(seqs[c]) - 1000)
            if all(abs(p - q) > 600 for q in used[c]):
                break
        used[c].append(p)
        cp = mut(cons)
        st = random.choice("+-")
        if st == "-":
            cp = rc(cp)
        seqs[c][p:p + len(cp)] = list(cp)
        truth.append((c, p, p + len(cp), fam, st))
with open("genome.fa", "w") as fh:
    for c, s in seqs.items():
        s = "".join(s)
        fh.write(">%s test sequence\n" % c)
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")
with open("bank.fa", "w") as fh:
    fh.write(">famA\n%s\n>famB\n%s\n" % (famA, famB))
with open("truth.bed", "w") as fh:
    for t in sorted(truth):
        fh.write("%s\t%d\t%d\t%s\t0\t%s\n" % t)
with open("maskA.bed", "w") as fh:
    for c, s, e, f, st in sorted(truth):
        if f == "A":
            fh.write("%s\t%d\t%d\n" % (c, s, e))
with open("bad_name.bed", "w") as fh:
    fh.write("chr_three\t10\t20\n")
with open("bad_end.bed", "w") as fh:
    fh.write("chr_two\t299990\t300010\n")
EOF

echo "== 1. mask_genome.py guards"
python3 "$SD/tools/mask_genome.py" genome.fa bad_name.bed -o x.fa 2>/dev/null && fail "accepted an unknown sequence name"
python3 "$SD/tools/mask_genome.py" genome.fa bad_end.bed -o x.fa 2>/dev/null && fail "accepted an interval past the end"
[[ ! -e x.fa ]] || fail "left an output file after refusing"
python3 "$SD/tools/mask_genome.py" genome.fa maskA.bed -o masked.fa --stats stats.tsv
cmp <(grep '>' genome.fa) <(grep '>' masked.fa) || fail "headers changed"
cmp <(awk '/^>/{if(n)print n; n=0; next}{n+=length($0)} END{print n}' genome.fa) \
    <(awk '/^>/{if(n)print n; n=0; next}{n+=length($0)} END{print n}' masked.fa) || fail "lengths changed"
echo "ok"

count(){  # count(run_dir, family) -> hits in step1 labelled BED
  awk -v f="$2" -F'\t' '$4==f' "$1"/genome.clean_step1/searches/all_hits.labeled.bed | wc -l
}

echo "== 2. unmasked control"
# later steps can stumble on a tiny synthetic genome; this test is about step1, so judge step1 only
mkdir -p plain && ( cd plain && SINEderella ../genome.fa ../bank.fa > log.txt 2>&1 ) || true
P=$(ls -d plain/run_* | tail -1)
[[ -s "$P"/genome.clean_step1/searches/all_hits.labeled.bed ]] || fail "unmasked step1 produced no hits (plain/log.txt)"
a0=$(count "$P" famA); b0=$(count "$P" famB)
echo "famA $a0 / 80, famB $b0 / 60"
(( a0 >= 76 && b0 >= 57 )) || fail "control did not find the planted copies"

echo "== 3. masked run"
mkdir -p masked && ( cd masked && SINEderella --mask-bed ../maskA.bed ../genome.fa ../bank.fa > log.txt 2>&1 ) || true
M=$(ls -d masked/run_* | tail -1)
[[ -s "$M"/genome.clean_step1/searches/all_hits.labeled.bed ]] || fail "masked step1 produced no hits (masked/log.txt)"
a1=$(count "$M" famA); b1=$(count "$M" famB)
echo "famA $a1 / 80 (masked), famB $b1 / 60"
(( a1 == 0 )) || fail "masked family still found: $a1 hits"
(( b1 == b0 )) || fail "unmasked family changed: $b1 vs $b0 in the control"
grep -q '^MASK_BP' "$M/manifest.txt" || fail "manifest has no MASK_BP"
# B hits at planted coordinates (names sanitised: _ -> @U@)
awk -F'\t' '$4=="B"{gsub(/_/,"@U@",$1); print $1"\t"$2"\t"$3}' truth.bed | sort -k1,1 -k2,2n > truthB.bed
awk -F'\t' '$4=="famB"{print $1"\t"$2"\t"$3}' "$M"/genome.clean_step1/searches/all_hits.labeled.bed | sort -k1,1 -k2,2n > hitsB.bed
hit=$(bedtools intersect -u -f 0.9 -r -a truthB.bed -b hitsB.bed | wc -l)
echo "famB planted copies recovered at >= 90% reciprocal overlap: $hit / 60"
(( hit == b0 )) || fail "famB coordinates do not match the planted copies"
# extraction is from the unmasked genome: no N in any extracted copy
! grep -v '^>' "$M"/genome.clean_step1/extracted.fasta | grep -q N || fail "N in extracted sequence"
echo "PASS: masked search"
