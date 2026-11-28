#!/usr/bin/env bash
# tools/piece_only_drop.sh on a planted run (bash tests/test_piece_only_drop.sh; needs bedtools). Prints PASS/FAIL per case.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); TOOL=$HERE/../tools/piece_only_drop.sh
T=$(mktemp -d "${TMPDIR:-/tmp}/t_pod.XXXXXX"); trap 'rm -rf "$T"' EXIT
fail=0; chk(){ if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fail=1; fi; }
mk(){   # RUN with_family_view
  R=$T/$1; mkdir -p "$R/step2/step2_output" "$R/genome.clean_step1/searches" "$R/family_view"
  printf 'copy\tsubfamily\tassignment_status\tverdict\treason\n' > "$R/step2/step2_output/loci_verdict.tsv"
  printf '%s\n' \
    "c1:100-200(+)	A	rejected_low_bitscore	not_SINE	not assigned" \
    "c1:1000-1069(+)	A	assigned	SINE	ok" \
    "c1:2000-2069(-)	A	rejected_low_bitscore	not_SINE	not assigned" \
    "c1:3000-3069(+)	A	no_unanimous	not_SINE	not assigned" \
    "c1:4000-4069(+)	A	rejected_low_bitscore	not_SINE	not assigned" \
    "ch:r:5000-5069(+)	A	rejected_low_bitscore	not_SINE	not assigned" >> "$R/step2/step2_output/loci_verdict.tsv"
  printf 'c1\t120\t180\tA\n' > "$R/genome.clean_step1/searches/gen-consensuses.clean.part_A.bed"
  printf 'c1\t1000\t1069\nc1\t2000\t2069\nc1\t3000\t3069\nc1\t4000\t4069\nch:r\t5000\t5069\n' > "$R/genome.clean_step1/searches/gen-pieces.part_A__piece_1-69.bed"
  if [[ $2 == 1 ]]; then
    printf 'Sequence\tFamily\tBest\tStretch\tCoverage_of_best\tFamily_score\tStretch_threshold\tStatus\tStretch_id\n' > "$R/family_view/partial_copies.tsv"
    printf 'c1:2000-2069(-)\tF\tA\tA:1-69\t0.2\t1\t1\tpartial\tP1\nc1:3000-3069(+)\tF\tA\tA:1-69\t0.2\t1\t1\tbelow_stretch_threshold\tP1\n' >> "$R/family_view/partial_copies.tsv"
  fi
}
ids(){ awk -F'\t' 'NR > 1 { print $1 }' "$1" | tr '\n' ' '; }

mk fam 1; bash "$TOOL" "$T/fam" 2> "$T/fam.log"; L=$T/fam/step2/step2_output
chk "whole hit, unassigned: kept"            "[[ \"\$(ids $L/loci_verdict.tsv)\" == *'c1:100-200(+)'* ]]"
chk "piece-only assigned: kept"              "[[ \"\$(ids $L/loci_verdict.tsv)\" == *'c1:1000-1069(+)'* ]]"
chk "piece-only partial: kept"               "[[ \"\$(ids $L/loci_verdict.tsv)\" == *'c1:2000-2069(-)'* ]]"
chk "piece-only below stretch bar: dropped"  "[[ \"\$(ids $L/piece_only_dropped.tsv)\" == *'c1:3000-3069(+)'* ]]"
chk "piece-only rejected: dropped"           "[[ \"\$(ids $L/piece_only_dropped.tsv)\" == *'c1:4000-4069(+)'* ]]"
chk "contig name with a colon: dropped"      "[[ \"\$(ids $L/piece_only_dropped.tsv)\" == *'ch:r:5000-5069(+)'* ]]"
chk "kept + dropped = all, header kept"      "[[ \$(wc -l < $L/loci_verdict.tsv) == 4 && \$(wc -l < $L/piece_only_dropped.tsv) == 4 && \$(head -1 $L/loci_verdict.tsv | cut -f1) == copy ]]"
chk "log line"                               "awk 'index(\$0, \"PIECE_ONLY piece_only=5 dropped=3 kept_assigned=1 kept_partial=1\") { f = 1 } END { exit !f }' $T/fam.log"
cp $L/loci_verdict.tsv $T/lv1; cp $L/piece_only_dropped.tsv $T/d1
bash "$TOOL" "$T/fam" 2> /dev/null
chk "second call changes nothing"            "cmp -s $T/lv1 $L/loci_verdict.tsv && cmp -s $T/d1 $L/piece_only_dropped.tsv"

mk nofam 0; bash "$TOOL" "$T/nofam" 2> /dev/null; L=$T/nofam/step2/step2_output
chk "no family_view: partial copy dropped"   "[[ \"\$(ids $L/piece_only_dropped.tsv)\" == *'c1:2000-2069(-)'* ]]"

mk off 1; L=$T/off/step2/step2_output; cp $L/loci_verdict.tsv $T/lv0; PIECE_ONLY_DROP=0 bash "$TOOL" "$T/off" 2> /dev/null
chk "PIECE_ONLY_DROP=0: unchanged"           "cmp -s $T/lv0 $L/loci_verdict.tsv && [[ ! -e $L/piece_only_dropped.tsv ]]"

mk nopiece 1; rm "$T/nopiece/genome.clean_step1/searches/gen-pieces.part_A__piece_1-69.bed"; L=$T/nopiece/step2/step2_output
cp $L/loci_verdict.tsv $T/lv2; bash "$TOOL" "$T/nopiece" 2> /dev/null
chk "no pieces searched: unchanged"          "cmp -s $T/lv2 $L/loci_verdict.tsv && [[ \$(wc -l < $L/piece_only_dropped.tsv) == 1 ]]"
exit $fail
