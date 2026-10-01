#!/usr/bin/env bash
# fs6b_ends.sh OUT_DIR RUN_DIR [THREADS=8]
#
# Stage 6b of flankscan: are the ENDS of each kept candidate closed? A candidate is built from a pair of
# units (stages 3-5 pair every copy with its nearest neighbour per side), so in a chain of three or more
# units (rsi r1 + r3 + r3) the pairs are found separately, and a pair built alone is cut short: its copies
# still continue into a known unit beyond the end, yet stage 6 accepts it because its copies read as full
# units of it. Found on rsi P12 (r3 + r3, 153 copies, "accept" 88 %): its 5' flank is r1 in 57 / 60 copies.
# The check itself is fs_endcheck in fs_lib.sh (flank vs bank at the junction, open when >= 50 % of copies);
# stage 6c extends the candidates whose end is open.
#
# In : OUT/candidates.tsv, OUT/cand/NAME.bed, OUT/cons.masked.fa ; RUN_DIR/genome.clean.fa(.fai)
# Out: OUT/ends.tsv   name  side5_share side5_unit side5_second side3_share side3_unit side3_second open
#                      (share = fraction of the candidate's copies open on that side; open = 5 / 3 / 53 / -)
#      OUT/cand/ends/NAME.h5, .h3   the per-copy junction hits (stage 6c reads them)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/fs_lib.sh"
OUT=${1:?OUT_DIR}; RUN=$(readlink -f "${2:?RUN_DIR}"); T=${3:-8}
cd "$OUT"
printf "name\tside5_share\tside5_unit\tside5_second\tside3_share\tside3_unit\tside3_second\topen\n" > ends.tsv
gawk -F'\t' 'NR > 1 && $14 == "kept" { print $1 }' candidates.tsv | while read -r N; do
    [[ -s "cand/$N.bed" ]] || continue
    fs_endcheck "$N" >> ends.tsv
done
echo "fs6b: $(($(wc -l < ends.tsv) - 1)) candidates checked, $(gawk -F'\t' 'NR > 1 && $8 != "-"' ends.tsv | wc -l) with an open end -> $OUT/ends.tsv" >&2
