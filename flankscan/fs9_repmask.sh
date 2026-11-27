#!/usr/bin/env bash
# fs9_repmask.sh rep  SEQS.tsv KCOUNT GCAP [THREADS] > rep20.tsv        repeated 20-mers that occur in the flanks
# fs9_repmask.sh mask SEQS.tsv REP20.tsv  [THREADS] > seqs.masked.tsv   flank bases covered by a repeated 20-mer become N
#
# The repeat filter of fs9_twins.sh in awk, the same rows as fs9_repmask.py (and as the FS9_AWK=1 reference):
#  rep : canonical 20-mers (the smaller of the 20-mer and its reverse complement, string order = A < C < G < T) of the flanks,
#        windows with an N skipped, that KCOUNT ("KMER COUNT", whitespace separated, jellyfish dump -c) lists with COUNT > GCAP;
#        the KCOUNT k-mer is upper-cased. Written one per KCOUNT line, in KCOUNT order. Exit 1 if KCOUNT has no line with two fields.
#  mask: SEQS rows (name TAB side TAB sequence) in input order, every base covered by a 20-mer of REP20 as N.
#  A 5' flank (side 5) is stored reversed (outward from the junction); its 20-mers are taken in genome orientation.
# How: the reversed and the complemented sequence of every row come from rev and tr (C tools; an awk per-base loop was the slow
# part), so the awk takes each window and its reverse complement with two substr. The rows are cut into THREADS consecutive
# chunks, one awk each; mask outputs are joined in chunk order. rep: each awk holds the KCOUNT k-mers above GCAP in a hash and
# prints the flank k-mers found there; a last awk prints the KCOUNT lines whose k-mer was found, in KCOUNT order.
# Env: FS9_REP_AWK (default mawk if on PATH, else gawk; both give the same rows), TMPDIR (chunk files, about the size of SEQS).
set -euo pipefail
export LC_ALL=C
MODE=${1:-}; SQ=${2:-}
case $MODE in
    rep)  KC=${3:?KCOUNT}; C=${4:?GCAP}; P=${5:-1} ;;
    mask) REP=${3:?REP20.tsv}; P=${4:-1} ;;
    *) echo "usage: fs9_repmask.sh rep SEQS KCOUNT GCAP [THREADS] | mask SEQS REP20 [THREADS]" >&2; exit 2 ;;
esac
[[ -r $SQ ]] || { echo "fs9_repmask.sh: cannot read $SQ" >&2; exit 1; }
AWK=${FS9_REP_AWK:-$(command -v mawk || command -v gawk)}
[[ $P =~ ^[1-9][0-9]*$ ]] || P=1
TD=$(mktemp -d "${TMPDIR:-/tmp}/fs9_repmask.XXXXXX"); trap 'rm -rf "$TD"' EXIT
# chunks: consecutive rows c000, c001, ... (the order of SEQS)
NR_=$(wc -l < "$SQ"); (( NR_ > 0 )) || P=1
"$AWK" -v P="$P" -v N="$NR_" -v D="$TD" '{ print > sprintf("%s/c%03d", D, int((FNR - 1) * P / N)) }' "$SQ"
[[ -e $TD/c000 ]] || : > "$TD/c000"
# per row: side, stored sequence, its reverse, its complement, the complement of its reverse; name kept for mask
cols() {   # CHUNK
    cut -f3 "$1" > "$1.s"; rev "$1.s" > "$1.r"
    paste <(cut -f1,2 "$1") "$1.s" "$1.r" <(tr ACGT TGCA < "$1.s") <(tr ACGT TGCA < "$1.r")
    rm -f "$1.s" "$1.r"
}
# the genome-oriented sequence g and its reverse complement r: side 3 g = s, r = comp(rev s); side 5 g = rev s, r = comp(s).
# window i of g (1-based) and its reverse complement substr(r, n - i - 18, 20)
WIN='{ if ($2 == 5) { g = $4; r = $5 } else { g = $3; r = $6 }; n = length(g); hasn = index(g, "N") }'
if [[ $MODE == rep ]]; then
    # 1) per chunk: flank k-mers that are KCOUNT k-mers above GCAP, each once per chunk
    run_rep() {   # CHUNK
        cols "$1" | "$AWK" -v C="$C" '
            FILENAME == ARGV[1] { if (NF >= 2 && $2 + 0 > C) kc[toupper($1)] = 1; next }   # KCOUNT: blank-separated
            !tabs { FS = "\t"; tabs = 1; $0 = $0 }                                         # the flank rows: tab-separated
            '"$WIN"'
            { for (i = 1; i + 19 <= n; i++) { m = substr(g, i, 20); if (hasn && index(m, "N")) continue
                  x = substr(r, n - i - 18, 20); if (x < m) m = x
                  if ((m in kc) && !(m in hit)) { hit[m] = 1; print m } } }' "$KC" - > "$1.hit"
    }
    pids=(); for c in "$TD"/c[0-9][0-9][0-9]; do run_rep "$c" & pids+=($!); done
    rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
    (( rc == 0 )) || { echo "fs9_repmask.sh: rep chunk failed" >&2; exit 1; }
    # 2) KCOUNT lines whose upper-cased k-mer was found, in KCOUNT order; no line with two fields: exit 1 (wrong file)
    cat "$TD"/c[0-9][0-9][0-9].hit > "$TD/hits"
    "$AWK" -v C="$C" 'FILENAME == ARGV[1] { h[$1]; next }
        NF >= 2 { any = 1; m = toupper($1); if ($2 + 0 > C && (m in h)) print m }
        END { if (!any) exit 1 }' "$TD/hits" "$KC"
else
    run_mask() {   # CHUNK: rows in order; covered runs built in genome orientation, written in stored orientation
        cols "$1" | "$AWK" -F'\t' -v REP="$REP" '
            BEGIN { while ((getline l < REP) > 0) if (length(l) == 20) rep[l] = 1; close(REP)
                    NN = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" }
            function ns(k) { while (length(NN) < k) NN = NN NN; return substr(NN, 1, k) }
            '"$WIN"'
            { nr = 0; e = 0                                     # merged covered runs [rs[k], re[k]] in g coordinates
              for (i = 1; i + 19 <= n; i++) { m = substr(g, i, 20); if (hasn && index(m, "N")) continue
                  x = substr(r, n - i - 18, 20); if (x < m) m = x
                  if (m in rep) { if (nr && i <= e + 1) { e = i + 19; re[nr] = e } else { rs[++nr] = i; e = i + 19; re[nr] = e } } }
              if (!nr) { print $1 "\t" $2 "\t" $3; next }
              o = ""; p = 1                                     # p: next stored position not yet written
              if ($2 == 5) { for (k = nr; k >= 1; k--) { a = n - re[k] + 1; b = n - rs[k] + 1; o = o substr($3, p, a - p) ns(b - a + 1); p = b + 1 } }
              else         { for (k = 1; k <= nr; k++) { a = rs[k]; b = re[k]; o = o substr($3, p, a - p) ns(b - a + 1); p = b + 1 } }
              print $1 "\t" $2 "\t" o substr($3, p) }' > "$1.out"
    }
    pids=(); for c in "$TD"/c[0-9][0-9][0-9]; do run_mask "$c" & pids+=($!); done
    rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
    (( rc == 0 )) || { echo "fs9_repmask.sh: mask chunk failed" >&2; exit 1; }
    for c in "$TD"/c[0-9][0-9][0-9]; do cat "$c.out"; done
fi
