#!/usr/bin/env bash
# test_fs9_parallel.sh [WORKDIR] - stage 9 run side by side (2026-11-26) gives the same files as one step after the other, and a
# failing tool in any parallel child stops it.
#  fs9_twins.sh: the unmasked pass runs beside the masked passes, pass 1 beside pass 2. Same files at FS9_THREADS 1 and 5 (SOFTMASK
#    0 and both). Stubs: gawk failing in the pass-2 / unmasked confirmation, python3 failing in pass-2 seeds / in the repeat mask:
#    fs9_twins.sh must exit non-zero and leave no copy_status.tsv.
#  fs9_run.sh: families at once (FS9_JOBS). FS9_JOBS=1 and 3 give byte-identical files and flank_twins.tsv in the family order;
#    a family whose fs9_twins.sh fails gets no row and a message and the stage exits 0 (as before, test_audit7 D16); a failing
#    sort of one family's copies stops the stage non-zero.
# Needs bedtools, samtools, gawk, python3 + numpy (body check off: BODYCHECK=0; jellyfish used by fs9_run.sh if on PATH).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD="$(dirname "$HERE")"; FS="$SD/flankscan"
W=${1:-${TMPDIR:-/tmp}/fs9_parallel}; rm -rf "$W"; mkdir -p "$W"; cd "$W" || exit 1
export LC_ALL=C BODYCHECK=0
FAILS=0; PASS=0
chk() { if eval "$2"; then echo "PASS $1"; PASS=$((PASS + 1)); else echo "FAIL $1"; FAILS=$((FAILS + 1)); fi; }
same_dirs() {   # DIR1 DIR2: the same file names and every file byte-identical (fs9.log aside)
    diff <(cd "$1" && find . -type f ! -name fs9.log | sort) <(cd "$2" && find . -type f ! -name fs9.log | sort) > /dev/null || return 1
    local f; while IFS= read -r f; do cmp -s "$1/$f" "$2/$f" || return 1; done < <(cd "$1" && find . -type f ! -name fs9.log)
}

python3 "$FS/tests/make_toy_flank.py" toy 1 > /dev/null
samtools faidx toy/genome.fa
REALPY=$(command -v python3); REALCONF=$(command -v mawk || command -v gawk)   # the confirmation awk of fs9_twins.sh

# ---------- fs9_twins.sh: thread counts
for SM in 0 both; do
    for T in 1 5; do
        FS9_THREADS=$T SOFTMASK=$SM FS9_KEEP=1 bash "$FS/fs9_twins.sh" toy/genome.fa toy/copies.bed "tw_${SM}_$T" > "tw_${SM}_$T.log" 2>&1 ||
            { echo "FAIL fs9_twins SOFTMASK=$SM FS9_THREADS=$T run"; tail -5 "tw_${SM}_$T.log"; FAILS=$((FAILS + 1)); }
    done
    chk "fs9_twins SOFTMASK=$SM: FS9_THREADS 1 and 5 give the same files ($(ls "tw_${SM}_1" | wc -l) files, $(( $(wc -l < "tw_${SM}_1/flank_groups.tsv") - 1 )) pairs)" \
        "same_dirs tw_${SM}_1 tw_${SM}_5"
    chk "fs9_twins SOFTMASK=$SM: no work copy of the unmasked flanks left" "[ ! -e tw_${SM}_1/seqs.nomask.tsv ]"
done
chk "fs9_twins: the toy has twin pairs (the comparison tests something)" "[ \$(wc -l < tw_0_1/flank_groups.tsv) -gt 5 ]"
# the confirmation in gawk and in mawk (the default when on PATH): the same files
if command -v mawk > /dev/null; then
    for A in gawk mawk; do
        FS9_CONF_AWK=$(command -v $A) FS9_THREADS=3 SOFTMASK=both FS9_KEEP=1 bash "$FS/fs9_twins.sh" toy/genome.fa toy/copies.bed "aw_$A" > "aw_$A.log" 2>&1 ||
            { echo "FAIL fs9_twins FS9_CONF_AWK=$A run"; FAILS=$((FAILS + 1)); }
    done
    chk "fs9_twins: confirmation in gawk and in mawk give the same files ($(wc -l < aw_gawk/conf.p1.tsv) + $(wc -l < aw_gawk/conf.p2.tsv) confirmed rows)" "same_dirs aw_gawk aw_mawk"
else echo "SKIP confirmation gawk vs mawk (no mawk on PATH)"; fi

# ---------- fs9_twins.sh: a failing tool in a parallel child stops it
mkstub() {   # DIR NAME REAL PATTERN: NAME exits 3 when its arguments match PATTERN, else runs REAL
    mkdir -p "$1"; printf '#!/usr/bin/env bash\ncase "$*" in %s) echo "stub %s: failing on purpose" >&2; exit 3 ;; esac\nexec "%s" "$@"\n' "$4" "$2" "$3" > "$1/$2"
    chmod +x "$1/$2"
}
mkstub stub_conf_p2 confawk "$REALCONF" '*TAG=p2*'
mkstub stub_conf_nomask confawk "$REALCONF" '*TAG=nomask*'
mkstub stub_seeds_p2 python3 "$REALPY" '*fs9_seeds.py*111010010100110111*'
mkstub stub_repmask python3 "$REALPY" '*fs9_repmask.py\ mask*'
for S in stub_conf_p2 stub_conf_nomask stub_seeds_p2 stub_repmask; do
    CA=; [[ -e $W/$S/confawk ]] && CA=$W/$S/confawk
    PATH="$W/$S:$PATH" FS9_CONF_AWK=${CA:-$REALCONF} FS9_THREADS=4 bash "$FS/fs9_twins.sh" toy/genome.fa toy/copies.bed "out_$S" > "$S.log" 2>&1; RC=$?
    chk "fs9_twins: $S makes it fail (exit $RC, stub ran: $(gawk '/failing on purpose/ { n++ } END { print n + 0 }' "$S.log"))" \
        "[ $RC -ne 0 ] && [ ! -e out_$S/copy_status.tsv ] && gawk '/failing on purpose/ { f = 1 } END { exit !f }' $S.log"
    # stopped where the branches are joined, not only later by a missing file
    chk "fs9_twins: $S is caught where the parallel branches end" \
        "gawk '/fs9: (pass 1 \\/ pass 2|unmasked pass \\/ masked passes) failed/ { f = 1 } END { exit !f }' $S.log"
done

# ---------- fs9_run.sh: families at once
mkdir -p run/results; cp toy/genome.fa run/genome.clean.fa; samtools faidx run/genome.clean.fa
{ printf "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo\n"
  gawk -F'\t' '{ f = (NR == 7) ? "F4" : "F" (NR % 3 + 1); st = (NR % 11 == 0) ? "unassigned" : "assigned"
                 print $1 ":" $2 "-" $3 "(" $6 ")\t" f "\t500\t10\t" st "\t100\tNA" }' toy/copies.bed; } > run/results/assignment_full.tsv
for J in 1 3; do
    FS9_JOBS=$J bash "$FS/fs9_run.sh" run 4 > "run_$J.log" 2>&1; RC=$?
    chk "fs9_run FS9_JOBS=$J: exit 0" "[ $RC -eq 0 ]"
    rm -rf "res_$J"; cp -r run/results "res_$J"
done
chk "fs9_run: FS9_JOBS 1 and 3 give the same files ($(cd res_1 && find . -type f | wc -l) files)" "same_dirs res_1 res_3"
chk "fs9_run: rows in family order F1 F2 F3, F4 (1 copy) left out" \
    "[ \"\$(gawk -F'\t' 'NR > 1 { printf \"%s \", \$1 }' res_3/flank_twins.tsv)\" = 'F1 F2 F3 ' ]"
chk "fs9_run: no work directory left" "[ ! -e res_3/flank_twins.work ]"
# a family whose fs9_twins.sh fails: no row, a message, exit 0 (as before)
mkdir -p fsfail; cp "$FS"/*.sh "$FS"/*.py fsfail/
printf '#!/usr/bin/env bash\ncase "$2" in */F2/*) echo "fs9_twins stub: boom" >&2; exit 1 ;; esac\nexec bash "%s/fs9_twins.sh" "$@"\n' "$FS" > fsfail/fs9_twins.sh
FS9_JOBS=3 bash fsfail/fs9_run.sh run 4 > runfail.log 2>&1; RC=$?
chk "fs9_run: a failing family is reported, gets no row, the stage exits 0 (exit $RC)" \
    "[ $RC -eq 0 ] && gawk '/F2 failed/ { f = 1 } END { exit !f }' runfail.log && [ \"\$(gawk -F'\t' 'NR > 1 { printf \"%s \", \$1 }' run/results/flank_twins.tsv)\" = 'F1 F3 ' ]"
# a failing sort of one family's copies stops the stage
REALSORT=$(command -v sort)
mkstub stub_sort sort "$REALSORT" '*F3/copies.bed.unsorted*'
PATH="$W/stub_sort:$PATH" FS9_JOBS=3 bash "$FS/fs9_run.sh" run 4 > runsort.log 2>&1; RC=$?
chk "fs9_run: a failing sort in one family's job stops the stage non-zero (exit $RC)" \
    "[ $RC -ne 0 ] && gawk '/stub sort: failing on purpose/ { f = 1 } END { exit !f }' runsort.log"

echo "test_fs9_parallel: $PASS PASS, $FAILS FAIL"
(( FAILS == 0 ))
