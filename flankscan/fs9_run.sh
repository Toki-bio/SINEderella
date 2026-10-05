#!/usr/bin/env bash
# fs9_run.sh RUN_DIR [THREADS=8] [MIN_COPIES=20]
#
# Flank uniqueness over all firmly assigned copies of every family of a SINEderella run (flankscan stage 9, fs9_twins.sh),
# run by the orchestrator after assignment. Independent insertions have unrelated flanks; copies that were multiplied together
# with their neighbourhood (segmental duplications, array units the satellite screen did not catch, copies carried inside
# another element) share them. Found necessary on rsi MEG-RS (2026-10-05): of the copies left after the satellite screen,
# groups on four contigs had 96-100 % identical flanks, and the plates showed them as independent copies.
#
# In : RUN_DIR/genome.clean.fa (the run's genome; its headers are the sanitized ones used in results/assignment_full.tsv),
#      RUN_DIR/results/assignment_full.tsv (col 1 locus contig:start-end(strand), col 2 family, col 5 "assigned")
# Out: RUN_DIR/results/flank_twins/<family>/copy_status.tsv, flank_groups.tsv (fs9_twins.sh output)
#      RUN_DIR/results/flank_twins.tsv   family copies twin1 twin2 array masked untestable unique pct_twin
#      RUN_DIR/genome.clean.fa.k20gt20.tsv  canonical 20-mers seen more than 20 times (jellyfish; kept, it is the repeat mask of
#                                           every stage-9 run on this genome; genomes <= 300 Mb are counted by fs9_twins itself)
# The run's genome has no soft-masking (the header cleaning writes upper case), so the k-mer mask is the only repeat filter here.
# Env: FS9_SKIP=1 skips; MIN_COPIES; anything fs9_twins.sh reads (W, MINLEN, ID1, ...).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN=$(readlink -f "${1:?RUN_DIR}"); T=${2:-8}; MINC=${3:-${MIN_COPIES:-20}}
G="$RUN/genome.clean.fa"; AF="$RUN/results/assignment_full.tsv"; OUT="$RUN/results/flank_twins"
[[ -s "$G" && -s "$AF" ]] || { echo "fs9_run: need $G and $AF" >&2; exit 1; }
[[ -s "$G.fai" ]] || samtools faidx "$G"
mkdir -p "$OUT"
export LC_ALL=C

# repeat mask: canonical 20-mers above 20 copies; jellyfish for genomes over 300 Mb (fs9_twins counts smaller ones in awk)
SZ=$(gawk '{ s += $2 } END { print s + 0 }' "$G.fai")
KC="$G.k20gt20.tsv"
if (( SZ > 300000000 )) && [[ ! -s "$KC" ]]; then
    command -v jellyfish > /dev/null || { echo "fs9_run: genome > 300 Mb and no jellyfish on PATH: the 20-mer count cannot be made, stage 9 skipped" >&2; exit 0; }
    echo "fs9_run: counting canonical 20-mers of $(( SZ / 1000000 )) Mb with jellyfish ($T threads)"
    jellyfish count -m 20 -s 4G -t "$T" -C -o "$G.k20.jf" "$G"
    jellyfish dump -c -L 21 "$G.k20.jf" > "$KC"
    rm -f "$G.k20.jf"
    echo "fs9_run: $(wc -l < "$KC") 20-mers above 20 copies -> $KC"
fi

printf "family\tcopies\ttwin1\ttwin2\tarray\tmasked\tuntestable\tunique\tpct_twin\n" > "$OUT.tsv"
gawk -F'\t' 'NR > 1 && $5 == "assigned" { print $2 }' "$AF" | sort | uniq -c | gawk -v M="$MINC" '$1 >= M { print $2 }' | while read -r FAM; do
    D="$OUT/$FAM"; mkdir -p "$D"
    gawk -F'\t' -v OFS='\t' -v F="$FAM" 'NR > 1 && $2 == F && $5 == "assigned" {
        if (match($1, /^(.+):([0-9]+)-([0-9]+)\(([+-])\)$/, m)) print m[1], m[2], m[3], $1, 0, m[4] }' "$AF" | sort -k1,1 -k2,2n > "$D/copies.bed"
    if (( SZ > 300000000 )); then KCOUNT="$KC" bash "$HERE/fs9_twins.sh" "$G" "$D/copies.bed" "$D" > "$D/fs9.log" 2>&1 || { echo "fs9_run: $FAM failed, see $D/fs9.log" >&2; continue; }
    else bash "$HERE/fs9_twins.sh" "$G" "$D/copies.bed" "$D" > "$D/fs9.log" 2>&1 || { echo "fs9_run: $FAM failed, see $D/fs9.log" >&2; continue; }
    fi
    gawk -F'\t' -v OFS='\t' -v F="$FAM" '{ n[$2]++; t++ }
        END { printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\n", F, t, n["twin1"], n["twin2"], n["array"], n["masked"], n["untestable"], n["unique"],
                     t ? 100 * (n["twin1"] + n["twin2"]) / t : 0 }' "$D/copy_status.tsv" >> "$OUT.tsv"
    rm -f "$D"/keys.p?.tsv "$D"/hits.p?.tsv "$D"/f5.fa "$D"/f3.fa "$D"/seqs.tsv "$D"/rep20.tsv
done
echo "fs9_run: $(( $(wc -l < "$OUT.tsv") - 1 )) families -> $OUT.tsv"
column -t "$OUT.tsv" 2> /dev/null || cat "$OUT.tsv"
