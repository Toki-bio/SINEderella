#!/usr/bin/env bash
# test_fs9_engines.sh [WORKDIR] - fs9_twins.sh gives the same files with its default engines (fs9_repmask.py, fs9_seeds.py,
# confirmation in FS9_THREADS ranges) as with the awk reference (FS9_AWK=1, one awk for the confirmation).
#  toy 1: flankscan/tests/make_toy_flank.py (every flank case planted: twins, arrays, hot spots, soft-masked duplications),
#         run with SOFTMASK=both so both masking modes are compared
#  toy 2: a stress toy for the seed buckets: 400 copies whose flanks are drawn from a few shared blocks with point mutations,
#         some blocks at several offsets, so buckets exceed CAP, pairs have seeds on several diagonals (tie-breaks of bd),
#         copy names of different lengths (byte order), a KCOUNT file with repeated 20-mers, flanks cut by contig ends
# Compared byte for byte: rep20.tsv, seqs.tsv, cand.p1/p2, conf.p1/p2, pairs*.body, flank_groups.tsv, copy_status.tsv.
# Needs bedtools, samtools, gawk, python3 + numpy (BODYCHECK=0: the body check reads the same files in both modes).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD="$(dirname "$HERE")"
W=${1:-${TMPDIR:-/tmp}/fs9_engines}; rm -rf "$W"; mkdir -p "$W"
export LC_ALL=C BODYCHECK=0
FAILS=0
cmp_dirs() {   # NAME DIR_AWK DIR_NEW
    local f n=0
    for f in rep20.tsv seqs.tsv cand.p1.tsv cand.p2.tsv conf.p1.tsv conf.p2.tsv pairs.body pairs.strict.body pairs.masked95.body \
             flank_groups.tsv copy_status.tsv masked.txt untestable.txt; do
        [[ -e "$2/$f" || -e "$3/$f" ]] || continue
        if cmp -s "$2/$f" "$3/$f"; then n=$((n + 1)); else echo "FAIL $1: $f differs ($(wc -l < "$2/$f") vs $(wc -l < "$3/$f") rows)"; FAILS=$((FAILS + 1)); fi
    done
    echo "$1: $n files identical; cand.p1 $(wc -l < "$3/cand.p1.tsv") rows, cand.p2 $(wc -l < "$3/cand.p2.tsv"), pairs $(wc -l < "$3/flank_groups.tsv")"
}

# toy 1
python3 "$SD/flankscan/tests/make_toy_flank.py" "$W/toy1" 1 > /dev/null
samtools faidx "$W/toy1/genome.fa"
for mode in awk new; do
    if [[ $mode == awk ]]; then E="FS9_AWK=1 FS9_THREADS=1"; else E="FS9_AWK=0 FS9_THREADS=3"; fi
    env $E SOFTMASK=both bash "$SD/flankscan/fs9_twins.sh" "$W/toy1/genome.fa" "$W/toy1/copies.bed" "$W/toy1/$mode" > "$W/toy1/$mode.log" 2>&1 ||
        { echo "FAIL toy1 $mode run"; tail -5 "$W/toy1/$mode.log"; exit 1; }
done
cmp_dirs "toy1 (make_toy_flank, SOFTMASK=both)" "$W/toy1/awk" "$W/toy1/new"

# toy 2
python3 - "$W/toy2" <<'PY'
import os, random, sys
out = sys.argv[1]; os.makedirs(out, exist_ok=True); r = random.Random(11)
B = "ACGT"
rs = lambda n: "".join(r.choice(B) for _ in range(n))
mut = lambda s, p: "".join(r.choice(B.replace(c, "")) if r.random() < p else c for c in s)
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
rc = lambda s: "".join(comp[c] for c in reversed(s))
blocks = [rs(160) for _ in range(6)]
body = rs(150)
ctg, bed, rep = {}, [], set()
for c in range(8):
    name = "c%d" % c if c % 2 else "chrom_%d_long" % c          # names of different lengths and shapes
    seq = [rs(r.randint(5, 60))]                                 # first copy close to the contig start: short 5' flank
    pos = len(seq[0])
    for k in range(50):
        b5 = blocks[r.randrange(6)] if r.random() < 0.7 else rs(160)
        b3 = blocks[r.randrange(6)] if r.random() < 0.7 else rs(160)
        sh5 = r.choice([0, 0, 0, 3, 7, 30]); sh3 = r.choice([0, 0, 0, 2, 5, 40])  # same block at several offsets
        f5 = mut(b5, r.choice([0.0, 0.02, 0.08, 0.2]))[sh5:sh5 + 120]
        f3 = mut(b3, r.choice([0.0, 0.02, 0.08, 0.2]))[sh3:sh3 + 120]
        el = mut(body, 0.1); strand = r.choice("+-")
        if strand == "-": f5, f3, el = rc(f3), rc(f5), rc(el)
        piece = f5 + el + f3
        st = pos + len(f5); en = st + len(el)
        bed.append((name, st, en, strand))
        seq.append(piece); pos += len(piece)
        if r.random() < 0.3: rep.add(piece[r.randrange(0, 100):][:20])
    if c == 7: seq.append(rs(10))                                    # last copy 10 bp from the contig end
    else: seq.append(rs(200))
    ctg[name] = "".join(seq)
with open(out + "/genome.fa", "w") as f:
    for n, s in ctg.items(): f.write(">%s\n%s\n" % (n, s))
with open(out + "/copies.bed", "w") as f:
    for n, s, e, st in sorted(bed): f.write("%s\t%d\t%d\t%s:%d-%d(%s)\t0\t%s\n" % (n, s, e, n, s, e, st, st))
with open(out + "/kcount.txt", "w") as f:                          # jellyfish dump -c style, space separated, some below GCAP
    for k in sorted(rep):
        if len(k) == 20: f.write("%s %d\n" % (min(k, rc(k)), r.choice([5, 25, 300])))
PY
samtools faidx "$W/toy2/genome.fa"
for mode in awk new; do
    if [[ $mode == awk ]]; then E="FS9_AWK=1 FS9_THREADS=1"; else E="FS9_AWK=0 FS9_THREADS=4"; fi
    env $E KCOUNT="$W/toy2/kcount.txt" CAP=6 bash "$SD/flankscan/fs9_twins.sh" "$W/toy2/genome.fa" "$W/toy2/copies.bed" "$W/toy2/$mode" > "$W/toy2/$mode.log" 2>&1 ||
        { echo "FAIL toy2 $mode run"; tail -5 "$W/toy2/$mode.log"; exit 1; }
done
cmp_dirs "toy2 (stress, CAP=6, KCOUNT)" "$W/toy2/awk" "$W/toy2/new"
# fs9_seeds.py itself with many small blocks (several workers, block boundaries inside the family)
python3 "$SD/flankscan/fs9_seeds.py" "$W/toy2/new/seqs.tsv" 111010010100110111 3 --cap 6 --threads 4 --block-hits 50 | sort -u > "$W/toy2/cand.p2.blocks"
if cmp -s "$W/toy2/cand.p2.blocks" "$W/toy2/awk/cand.p2.tsv"; then echo "toy2: fs9_seeds.py with 50-pair blocks identical"; else echo "FAIL toy2: fs9_seeds.py with small blocks differs"; FAILS=$((FAILS + 1)); fi
[[ $(wc -l < "$W/toy2/awk/cand.p2.tsv") -gt 0 ]] || { echo "FAIL toy2: no candidates at all, the stress toy tests nothing"; FAILS=$((FAILS + 1)); }
(( FAILS == 0 )) && echo "test_fs9_engines: passed" || { echo "test_fs9_engines: $FAILS failures"; exit 1; }
