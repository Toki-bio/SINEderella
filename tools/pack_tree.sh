#!/usr/bin/env bash
# Pack the SINEderella tree for a machine that has no git, with the commit inside it.
#
#   tools/pack_tree.sh [OUT.tar.gz] [COMMITTISH]
#
# Default output: ../SINEderella_<commit>.tar.gz next to the checkout; default committish: HEAD.
# The tarball holds the tree of that commit (git archive, so no .git and no untracked files) plus a
# PIPELINE_COMMIT file. tools/pipeline_commit.sh reads that file when git is not available, so a run
# from the unpacked tree records the real commit instead of "unknown"
# (docs/REPRODUCIBILITY.md section 12).
#
# Unpack and run:
#   mkdir SINEderella && tar -xzf SINEderella_<commit>.tar.gz -C SINEderella
#   SINEderella/SINEderella genome.fa consensuses.fa
set -euo pipefail

REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO" || exit 1

[[ -e "$REPO/.git" ]] || { echo "pack_tree: $REPO is not a git checkout" >&2; exit 1; }

REF="${2:-HEAD}"
COMMIT="$(git -C "$REPO" describe --always --dirty "$REF" 2>/dev/null || git -C "$REPO" rev-parse --short "$REF")"
# --dirty only means something for HEAD; for an explicit ref git describe ignores it
if [[ $REF != "HEAD" ]]; then COMMIT="${COMMIT%-dirty}"; fi
OUT="${1:-$REPO/../SINEderella_${COMMIT}.tar.gz}"

if [[ $COMMIT == *-dirty ]]; then
  echo "pack_tree: the checkout has uncommitted changes; git archive packs the committed tree, so the" >&2
  echo "           tarball will NOT contain them. Commit first, or pack a clean checkout." >&2
  exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
git -C "$REPO" archive --format=tar "$REF" -o "$TMP/tree.tar"
printf '%s\n' "$COMMIT" > "$TMP/PIPELINE_COMMIT"
tar -rf "$TMP/tree.tar" -C "$TMP" PIPELINE_COMMIT
gzip -c "$TMP/tree.tar" > "$OUT"

echo "pack_tree: $OUT" >&2
echo "pack_tree: PIPELINE_COMMIT $COMMIT" >&2
