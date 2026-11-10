#!/usr/bin/env bash
# fs9_run.sh RUN_DIR [THREADS=8] [MIN_COPIES=2]
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
#                                           every stage-9 run on this genome; made again when older than the genome)
# The run's genome has no soft-masking (the header cleaning writes upper case), so the k-mer mask is the only repeat filter here.
# Env: FS9_SKIP=1 skips; MIN_COPIES; anything fs9_twins.sh reads (W, MINLEN, ID1, ...).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# MIN_COPIES 2: any two copies can be compared; until 2026-10-21 it was 20, and smaller families never got a copy verdict
RUN=$(readlink -f "${1:?RUN_DIR}"); T=${2:-8}; MINC=${3:-${MIN_COPIES:-2}}
export BODYCHECK_THREADS=${BODYCHECK_THREADS:-$T}   # parallel ssearch36 calls in fs9_bodycheck.py
export FS9_THREADS=${FS9_THREADS:-$T}               # workers of fs9_seeds.py and of the confirmation awks
G="$RUN/genome.clean.fa"; AF="$RUN/results/assignment_full.tsv"; OUT="$RUN/results/flank_twins"
[[ -s "$G" && -s "$AF" ]] || { echo "fs9_run: need $G and $AF" >&2; exit 1; }
[[ -s "$G.fai" ]] || samtools faidx "$G"
# results of an earlier stage-9 run in this run dir go first: a family that fails now (or a stage skipped for want of
# jellyfish) kept its old copy_status.tsv / flank_groups.tsv, read by flank_verdict.py and step8a as this run's result; with
# none, its copies are grey (GLM audit 2026-11-08, finding D16)
rm -rf "$OUT" "$OUT.tsv"
mkdir -p "$OUT"
export LC_ALL=C

# repeat mask: canonical 20-mers above 20 copies, counted ONCE per genome with jellyfish and passed to every family (KCOUNT).
# Until 2026-11-10 only genomes over 300 Mb were counted here; for a smaller one each family's fs9_twins.sh counted the whole
# genome again in awk: 6.7 s and 0.31 GB of memory per Mb (2 and 4 Mb measured), i.e. about 31 GB per family at 100 Mb.
# Without jellyfish: genomes up to 300 Mb fall back to that awk count (a warning says so), larger ones skip stage 9 as before.
SZ=$(gawk '{ s += $2 } END { print s + 0 }' "$G.fai")
KC="$G.k20gt20.tsv"; USEKC=1
if [[ ! -s "$KC" || "$KC" -ot "$G" ]]; then
    if command -v jellyfish > /dev/null; then
        echo "fs9_run: counting canonical 20-mers of $(( SZ / 1000000 )) Mb with jellyfish ($T threads)"
        JS=4G; (( SZ > 300000000 )) || JS=$(( SZ > 1000000 ? SZ : 1000000 ))      # hash size: as before for large genomes
        jellyfish count -m 20 -s "$JS" -t "$T" -C -o "$G.k20.jf" "$G"
        jellyfish dump -c -L 21 "$G.k20.jf" > "$KC.tmp"
        # a genome with no 20-mer above 20 copies gives an empty dump; fs9_twins.sh needs a KCOUNT with a line in it
        [[ -s "$KC.tmp" ]] || echo "NNNNNNNNNNNNNNNNNNNN 0" > "$KC.tmp"
        mv "$KC.tmp" "$KC"; rm -f "$G.k20.jf"
        echo "fs9_run: $(wc -l < "$KC") 20-mers above 20 copies -> $KC"
    elif (( SZ > 300000000 )); then
        echo "fs9_run: genome > 300 Mb and no jellyfish on PATH: the 20-mer count cannot be made, stage 9 skipped" >&2; exit 0
    else
        echo "fs9_run: no jellyfish on PATH: every family counts the genome's 20-mers in awk (about 0.3 GB of memory and 7 s per Mb)" >&2
        USEKC=0
    fi
fi

printf "family\tcopies\ttwin1\ttwin2\tarray\tmasked\tuntestable\tunique\tpct_twin\n" > "$OUT.tsv"
gawk -F'\t' 'NR > 1 && $5 == "assigned" { print $2 }' "$AF" | sort | uniq -c | gawk -v M="$MINC" '$1 >= M { print $2 }' | while read -r FAM; do
    D="$OUT/$FAM"; mkdir -p "$D"
    gawk -F'\t' -v OFS='\t' -v F="$FAM" 'NR > 1 && $2 == F && $5 == "assigned" {
        if (match($1, /^(.+):([0-9]+)-([0-9]+)\(([+-])\)$/, m)) print m[1], m[2], m[3], $1, 0, m[4] }' "$AF" | sort -k1,1 -k2,2n > "$D/copies.bed"
    if (( USEKC )); then KCOUNT="$KC" bash "$HERE/fs9_twins.sh" "$G" "$D/copies.bed" "$D" > "$D/fs9.log" 2>&1 || { echo "fs9_run: $FAM failed, see $D/fs9.log" >&2; continue; }
    else bash "$HERE/fs9_twins.sh" "$G" "$D/copies.bed" "$D" > "$D/fs9.log" 2>&1 || { echo "fs9_run: $FAM failed, see $D/fs9.log" >&2; continue; }
    fi
    gawk -F'\t' -v OFS='\t' -v F="$FAM" '{ n[$2]++; t++ }
        END { printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\n", F, t, n["twin1"], n["twin2"], n["array"], n["masked"], n["untestable"], n["unique"],
                     t ? 100 * (n["twin1"] + n["twin2"]) / t : 0 }' "$D/copy_status.tsv" >> "$OUT.tsv"
    rm -f "$D"/keys.p?.tsv "$D"/hits.p?.tsv "$D"/f5.fa "$D"/f3.fa "$D"/seqs.tsv "$D"/rep20.tsv "$D"/flank20.tsv
    # intermediates no other step reads (2026-10-24; cand.p1/p2 were 0.6 GB on tbr VES); FS9_KEEP=1 keeps them for debugging
    [[ ${FS9_KEEP:-0} == 1 ]] || rm -f "$D"/cand.p?.tsv "$D"/conf.p?.tsv "$D"/cand.nomask.tsv "$D"/conf.nomask.tsv "$D"/pairs*.body "$D"/flank_groups.body "$D"/seqs.raw.tsv
done
echo "fs9_run: $(( $(wc -l < "$OUT.tsv") - 1 )) families -> $OUT.tsv"
column -t "$OUT.tsv" 2> /dev/null || cat "$OUT.tsv"
