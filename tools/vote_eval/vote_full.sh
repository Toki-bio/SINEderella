#!/bin/bash
# vote_full.sh CONS.fa COPIES.fa OUTDIR [MODE=base|pad] [PART=20000] [THREADS=32]
# The step-2 vote of SINEderella reproduced with the same ssearch36 call and the same library parts (seqkit split by PART), but
# every cycle keeps, per copy, the best alignment of each of its top 6 consensuses (not only the winner), with coordinates:
#   OUTDIR/cyc_<n>.tsv : copy  consensus  bits  cons_start  cons_end  copy_start  copy_end  part
# so the ten cycles can be re-tallied offline in any way (flat, family-then-subfamily, coverage rules). Read-only for the run.
set -euo pipefail
CONS=$1; COP=$2; OUT=$3; MODE=${4:-base}; PART=${5:-20000}; T=${6:-32}
mkdir -p "$OUT/parts"; rm -f "$OUT"/parts/* "$OUT"/cyc_*.tsv
seqkit split2 -s "$PART" -O "$OUT/parts" "$COP" 2> /dev/null
if [[ $MODE == pad ]]; then
  for f in "$OUT"/parts/*; do
    n=$(grep -c '>' "$f"); (( n < PART )) || continue
    python3 - "$COP" "$f" $(( PART - n )) <<'EOF'
import random, sys
src, part, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
seqs, buf = [], []
for l in open(src):
    l = l.rstrip("\n")
    if l.startswith(">"):
        if buf: seqs.append("".join(buf))
        buf = []
    else: buf.append(l)
if buf: seqs.append("".join(buf))
R = random.Random(11)
with open(part, "a") as o:
    for i in range(k):
        s = list(R.choice(seqs)); R.shuffle(s)
        o.write(">DECOY_%d\n%s\n" % (i, "".join(s)))
EOF
  done
fi
echo "vote_full $MODE: $(grep -c '>' "$COP") copies, parts: $(for f in "$OUT"/parts/*; do grep -c '>' $f; done | tr '\n' ' ')"
for c in $(seq 1 10); do
  for f in "$OUT"/parts/*; do
    ssearch36 -g -3 -T "$T" -Q -n -z 11 -E 2 -w 95 -W 70 -m 8 "$CONS" "$f" 2> /dev/null \
      | gawk -v P=$(basename $f) -v OFS='\t' '$2 !~ /^DECOY_/ { k = $2 SUBSEP $1; if (!(k in b) || $12 > b[k]) { b[k] = $12; r[k] = $7 OFS $8 OFS $9 OFS $10 } }
          END { for (k in b) { split(k, x, SUBSEP); print x[1], x[2], b[k], r[k], P } }' \
      | sort -t$'\t' -k1,1 -k3,3gr | gawk -F'\t' '$1 != p { p = $1; n = 0 } ++n <= 6' >> "$OUT/cyc_$c.tsv"
  done
done
echo "vote_full done: $(wc -l < "$OUT/cyc_1.tsv") rows in cycle 1"
