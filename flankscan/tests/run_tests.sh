#!/usr/bin/env bash
# run_tests.sh [STAGE ...] - build the toy (make_toy.sh), run flankscan stages on it, check each
# stage against the planted truth. PASS/FAIL per check; exit 1 if any FAIL.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; FS="$(dirname "$HERE")"
T=${TOY:-$HOME/tmp/fs_toy}; O=$T/fs
rm -rf "$T"; bash "$HERE/make_toy.sh" "$T" >&2
fails=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails+1)); fi; }

# ---------- stage 1: extraction
bash "$FS/fs1_extract.sh" "$T" "$O" 1000
N=$(($(wc -l < "$T/truth.tsv") - 1))
chk "fs1: one window per planted copy ($N)" "[ \$(grep -c '^>' $O/windows.fa) -eq $N ] && [ \$((\$(wc -l < $O/loci.tsv)-1)) -eq $N ]"
# the core inside every window must be exactly the planted core (extracted independently)
gawk -F'\t' 'NR>1 { print $4 "\t" $5 "\t" $6 "\t" $1 "\t0\t" $7 }' "$O/loci.tsv" > "$O/core.bed"
bedtools getfasta -fi "$T/genome.clean.fa" -bed "$O/core.bed" -s -nameOnly | sed '/^>/s/([+-])$//' \
  | gawk '/^>/{k=substr($0,2); next} {print k "\t" $0}' | sort > "$O/core.tsv"
gawk -F'\t' 'FNR==NR { if (NR>1) { cs[$1]=$8; ce[$1]=$9 }; next }
             /^>/ { k=substr($0,2); next }
             { print k "\t" substr($0, cs[k], ce[k]-cs[k]+1) }' "$O/loci.tsv" <(seqkit seq -w 0 "$O/windows.fa") \
  | sort > "$O/core_from_window.tsv"
chk "fs1: core inside each window = the planted core (all $N)" "cmp -s $O/core.tsv $O/core_from_window.tsv"
chk "fs1: contig-end copy has clamp5 = 800" \
    "gawk -F'\t' 'FNR==NR{if(\$2==\"contigend\")L=\$1;next} \$2==L{exit (\$11==800 && \$12==0)?0:1}' $T/truth.tsv $O/loci.tsv"
# dimer: TA head and TB part are two loci exactly 39 bp apart (the linker)
chk "fs1: every dimer_left has its TB neighbour 39 bp downstream (40)" \
    "[ \$(gawk -F'\t' 'FNR==NR{if(\$2==\"dimer_left\")d[\$1];next} (\$2 in d) && \$14==39 && \$16==\"TB\"' $T/truth.tsv $O/loci.tsv | wc -l) -eq 40 ]"

for s in "$@"; do [[ -f "$HERE/checks_$s.sh" ]] && source "$HERE/checks_$s.sh"; done
echo "fails=$fails"; exit $((fails > 0))
