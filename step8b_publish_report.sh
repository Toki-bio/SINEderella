#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# step8b_publish_report.sh
#
# Wire already-generated per-subfamily alignment files into an existing
# SINEderella HTML report so it becomes a complete, self-contained,
# publishable result. This script does NOT extract or align anything — it
# only reads results/alignments/manifest.tsv and edits results/report.html.
#
# Usage:
#   step8b_publish_report.sh <RUN_ROOT> <SPECIES_CODE> \
#       [--raw-base-url URL] [--msa-viewer-url URL]
###############################################################################

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ printf '[%s] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2; exit 1; }

usage(){
  cat >&2 <<'EOF'
Usage: step8b_publish_report.sh <RUN_ROOT> <SPECIES_CODE> \
       [--raw-base-url URL] [--msa-viewer-url URL]

  RUN_ROOT        SINEderella run directory containing results/report.html
                  and results/alignments/manifest.tsv
  SPECIES_CODE    Species prefix used in alignment filenames
  --raw-base-url URL
                  If given, alignment links open in an MSA viewer pointing
                  at this raw file base URL. If omitted, links point
                  directly at the local relative alignment file path.
  --msa-viewer-url URL
                  Base URL of the MSA viewer app (only used together with
                  --raw-base-url). Default:
                  https://toki-bio.github.io/MSA-viewer/
EOF
  exit 2
}

###############################################################################
# Parse arguments
###############################################################################
[[ $# -ge 2 ]] || usage

RUN_ROOT="$1"; shift
SPECIES_CODE="$1"; shift

RAW_BASE_URL=""
MSA_VIEWER_URL="https://toki-bio.github.io/MSA-viewer/"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw-base-url)
      shift
      [[ $# -gt 0 ]] || die "--raw-base-url requires a URL argument"
      RAW_BASE_URL="$1"
      shift
      ;;
    --msa-viewer-url)
      shift
      [[ $# -gt 0 ]] || die "--msa-viewer-url requires a URL argument"
      MSA_VIEWER_URL="$1"
      shift
      ;;
    -h|--help)
      usage
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

command -v python3 >/dev/null 2>&1 || die "python3 is required but not found in PATH"

REPORT="$RUN_ROOT/results/report.html"
MANIFEST="$RUN_ROOT/results/alignments/manifest.tsv"

[[ -f "$REPORT"  ]] || die "Missing required input: $REPORT"
[[ -f "$MANIFEST" ]] || die "Missing required input: $MANIFEST"

log "Publishing subfamily alignments into: $REPORT"

###############################################################################
# Build the alignments table fragment and splice it into report.html.
# Python is used for the in-place HTML edit so the multi-line insertion is
# robust (no fragile sed). newline="" preserves line endings byte-for-byte
# outside the changed region.
###############################################################################
python3 - "$REPORT" "$MANIFEST" "$SPECIES_CODE" "$RAW_BASE_URL" "$MSA_VIEWER_URL" <<'PYEOF'
import sys
import re
import html
import urllib.parse

report_path, manifest_path, species, raw_base, msa_url = sys.argv[1:6]

# --- Read manifest.tsv (do NOT scan the alignments/ directory) ---
rows = []
with open(manifest_path, encoding="utf-8", newline="") as f:
    for line in f:
        line = line.rstrip("\n").rstrip("\r")
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        # Skip a header row: the has_top100 column must be 0 or 1
        if cols[1] not in ("0", "1"):
            continue
        subfam = cols[0]
        has_top = cols[1] == "1"
        has_rand = cols[2] == "1"
        has_sub = cols[3] == "1"
        n_members = cols[4] if len(cols) > 4 else ""
        rows.append((subfam, has_top, has_rand, has_sub, n_members))

# --- Build a single link element for one alignment variant ---
def make_link(subfam, variant, filename, title):
    if raw_base:
        file_url = raw_base.rstrip("/") + "/" + filename
        href = (
            msa_url
            + "?url=" + urllib.parse.quote(file_url, safe="")
            + "&title=" + urllib.parse.quote(title, safe="")
        )
    else:
        href = "alignments/" + urllib.parse.quote(filename, safe="")
    return '<a href="{}">{}</a>'.format(
        html.escape(href, quote=True),
        html.escape(variant),
    )

def cell(link_html):
    return "<td>{}</td>".format(link_html) if link_html else "<td>&mdash;</td>"

# --- Build the HTML fragment: one table row per subfamily, in manifest order ---
body_rows = []
for subfam, has_top, has_rand, has_sub, n_members in rows:
    cells = ["<td>{}</td>".format(html.escape(subfam))]

    top_link = rand_link = sub_link = ""
    if has_top:
        top_link = make_link(
            subfam, "top100",
            "{}_{}_top100.aln.fa".format(species, subfam),
            "{} {} top100".format(species, subfam),
        )
    if has_rand:
        rand_link = make_link(
            subfam, "rand100",
            "{}_{}_rand100.aln.fa".format(species, subfam),
            "{} {} rand100".format(species, subfam),
        )
    if has_sub:
        sub_link = make_link(
            subfam, "subfam",
            "{}_{}_subfam.aln.fa".format(species, subfam),
            "{} {} subfam".format(species, subfam),
        )

    cells.append(cell(top_link))
    cells.append(cell(rand_link))
    cells.append(cell(sub_link))
    cells.append("<td>{}</td>".format(html.escape(n_members)))
    body_rows.append("<tr>{}</tr>".format("".join(cells)))

fragment = (
    '<section id="alignments">\n'
    '<h2>Subfamily Alignments</h2>\n'
    '<table>\n'
    '<thead><tr><th>Subfamily</th><th>top100</th>'
    '<th>rand100</th><th>subfam</th><th>n_members</th></tr></thead>\n'
    '<tbody>\n'
    + "\n".join(body_rows)
    + '\n</tbody>\n'
    '</table>\n'
    '</section>\n'
)

# --- Read the existing report (preserve bytes outside the changed region) ---
with open(report_path, encoding="utf-8", newline="") as f:
    content = f.read()

# Case-insensitive search via regex (avoids str.lower() length-shift issues
# with certain Unicode code points, which would misalign byte indices).
needle = "alignments not yet available"
m = re.search(re.escape(needle), content, re.IGNORECASE)

if m:
    # Replace the placeholder text with the fragment.
    content = content[:m.start()] + fragment + content[m.end():]
else:
    # No placeholder found: insert right before the closing </body> tag.
    body_matches = list(re.finditer(r"</body>", content, re.IGNORECASE))
    if body_matches:
        body_idx = body_matches[-1].start()
        content = content[:body_idx] + fragment + content[body_idx:]
    else:
        content = content + fragment

with open(report_path, "w", encoding="utf-8", newline="") as f:
    f.write(content)
PYEOF

log "Done: alignment links published into $REPORT"
