#!/usr/bin/env bash
# Print the commit of the SINEderella tree that is running, for manifest.txt and verdict_chain.done.
#
#   tools/pipeline_commit.sh [TREE]        TREE defaults to this script's parent directory
#
# A run must say which code made it, in every way the code is installed. In order:
#   1. $PIPELINE_COMMIT from the environment, if set and not empty. It wins, so a caller that knows
#      better (a build system, a wrapper) can say so.
#   2. git, but only when TREE itself is a checkout (TREE/.git exists). Without that test git walks up
#      and an unpacked tree sitting inside somebody else's repository would be stamped with that
#      repository's commit.
#   3. TREE/PIPELINE_COMMIT, a one-line file. The container build writes it (there is no git in the
#      image, containers/sinederella.def), and tools/pack_tree.sh writes it into the tarball.
#   4. unknown. The caller says so in its log; see docs/REPRODUCIBILITY.md section 12.
set -uo pipefail

tree="${1:-}"
if [[ -z $tree ]]; then
  tree="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fi

if [[ -n ${PIPELINE_COMMIT:-} ]]; then
  printf '%s\n' "$PIPELINE_COMMIT"
  exit 0
fi

if [[ -e "$tree/.git" ]]; then
  if c="$(git -C "$tree" describe --always --dirty 2>/dev/null)" && [[ -n $c ]]; then
    printf '%s\n' "$c"
    exit 0
  fi
fi

if [[ -s "$tree/PIPELINE_COMMIT" ]]; then
  # first non-empty line, CR stripped, no surrounding blanks
  c="$(awk '{ sub(/\r$/, ""); gsub(/^[ \t]+|[ \t]+$/, ""); if ($0 != "") { print; exit } }' "$tree/PIPELINE_COMMIT")"
  if [[ -n $c ]]; then
    printf '%s\n' "$c"
    exit 0
  fi
fi

printf 'unknown\n'
