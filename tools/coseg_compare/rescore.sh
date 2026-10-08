#!/usr/bin/env bash
# Rebuild scores.txt of a finished run directory from the files already there (no re-running of any tool).
#   rescore.sh WORKDIR [WORKDIR ...]
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PY=${PYTHON:-python3.12}
for W in "$@"; do
    ( cd "$W" || exit 1
      SPECS=()
      for a in cs_drop_m*/cs_drop.seqs.assign; do [ -f "$a" ] && SPECS+=("$a@cs_drop.names:COSEG_drop_$(basename "$(dirname "$a")" | sed 's/.*_//')"); done
      for a in cs_keep_m*/cs_keep.seqs.assign; do [ -f "$a" ] && SPECS+=("$a@cs_keep.names:COSEG_keep_$(basename "$(dirname "$a")" | sed 's/.*_//')"); done
      for n in $(ls -d sf_n* 2>/dev/null | sed 's/sf_n//' | sort -n -r); do [ -f "sf_n$n/copies.chunks.tsv" ] && SPECS+=("sf_n$n/copies.chunks.tsv:SubFam_n$n"); done
      for n in $(ls -d peel_n* 2>/dev/null | sed 's/peel_n//' | sort -n -r); do [ -s "peel_n$n/copy_groups.tsv" ] && SPECS+=("peel_n$n/copy_groups.tsv:SubFam+peel_n$n"); done
      $PY "$HERE/score.py" copies.labels "${SPECS[@]}" > scores.txt ) && echo "rescored $W" || echo "FAILED $W"
done
