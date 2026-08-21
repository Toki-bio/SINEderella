#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# step7_boundary_refine.sh
#
# Usage:
#   step7_boundary_refine.sh <RUN_ROOT> [FLANK_BASE=50] [STEP_BP=50] [MAX_EXT_BP=1000]
#
# This script refines SINE candidate boundaries by extending flanks until
# they become background-level similarity.
###############################################################################

RUN_ROOT="${1:-}"
FLANK_BASE="${2:-50}"
STEP_BP="${3:-50}"
MAX_EXT_BP="${4:-1000}"

if [[ -z "$RUN_ROOT" ]]; then
  echo "Usage: $0 <RUN_ROOT> [FLANK_BASE=50] [STEP_BP=50] [MAX_EXT_BP=1000]" >&2
  exit 1
fi

GENOME="$RUN_ROOT/genome.clean.fa"
FAI="$GENOME.fai"
ASSIGNED="$RUN_ROOT/step2/step2_output/assigned.fasta"
OUT_DIR="$RUN_ROOT/step2/step2_output"
OUT_TSV="$OUT_DIR/boundary_refinement.tsv"

for f in "$GENOME" "$FAI" "$ASSIGNED"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing required file: $f" >&2; exit 1; }
done

mkdir -p "$OUT_DIR"
TMP_DIR="$(mktemp -d "$RUN_ROOT/.boundary_refine_XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

log "Starting boundary refinement for $RUN_ROOT"

###############################################################################
# Step A: Background identity
###############################################################################
log "Step A: Measuring background identity..."
awk -F'\t' '$2 >= 200 {print $1"\t"$2}' "$FAI" > "$TMP_DIR/contigs.tsv"

cat > "$TMP_DIR/bg_regions.py" <<'EOF'
import random, sys
random.seed(42)
contigs = []
with open(sys.argv[1]) as f:
    for line in f:
        name, length = line.strip().split("\t")
        length = int(length)
        if length >= 200:
            contigs.append((name, length))
if not contigs:
    sys.exit(1)
regions = []
for _ in range(300):
    name, length = random.choice(contigs)
    start = random.randint(0, length - 70)
    end = start + 70
    regions.append(f"{name}:{start+1}-{end}")
with open(sys.argv[2], "w") as f:
    f.write("\n".join(regions) + "\n")
EOF

python3 "$TMP_DIR/bg_regions.py" "$TMP_DIR/contigs.tsv" "$TMP_DIR/bg_regions.txt"
samtools faidx "$GENOME" -r "$TMP_DIR/bg_regions.txt" > "$TMP_DIR/bg_seqs.fasta"

cat > "$TMP_DIR/calc_identity.py" <<'EOF'
import random, sys
random.seed(42)
seqs = []
with open(sys.argv[1]) as f:
    seq = ""
    for line in f:
        if line.startswith(">"):
            if seq:
                seqs.append(seq)
            seq = ""
        else:
            seq += line.strip()
    if seq:
        seqs.append(seq)
if len(seqs) < 2:
    print(0.0)
    sys.exit(0)
n_pairs = int(sys.argv[2]) if len(sys.argv) > 2 else 150
n_pairs = min(n_pairs, len(seqs) * (len(seqs)-1) // 2)
diffs = []
for _ in range(n_pairs):
    a, b = random.sample(seqs, 2)
    min_len = min(len(a), len(b))
    if min_len == 0:
        continue
    diff = sum(1 for i in range(min_len) if a[i] != b[i]) / min_len
    diffs.append(diff)
if not diffs:
    print(0.0)
else:
    identity = (1 - sum(diffs)/len(diffs)) * 100
    print(identity)
EOF

BG_IDENTITY=$(python3 "$TMP_DIR/calc_identity.py" "$TMP_DIR/bg_seqs.fasta" 2000)
IDENTITY_THRESHOLD=$(python3 -c "print($BG_IDENTITY + 10)")
log "Background identity: $BG_IDENTITY%, threshold: $IDENTITY_THRESHOLD%"

###############################################################################
# Step B: Per-subfamily boundary refinement
###############################################################################
log "Step B: Parsing assigned.fasta..."
awk '
/^>/ {
  hdr=$0
  sub(/^>/, "", hdr)
  n=split(hdr, a, "|")
  subfam=a[2]
  loc=a[1]
  strand="+"
  if (match(loc, /\([^\)]*\)/)) {
    strand=substr(loc, RSTART+1, RLENGTH-2)
    sub(/\([^\)]*\)/, "", loc)
  }
  colon_pos = 0
  for (i=1; i<=length(loc); i++) {
    if (substr(loc, i, 1) == ":") colon_pos = i
  }
  ctg = substr(loc, 1, colon_pos-1)
  coords = substr(loc, colon_pos+1)
  split(coords, c, "-")
  start=c[1]
  end=c[2]
  print subfam "\t" ctg "\t" start "\t" end "\t" strand
}
' "$ASSIGNED" > "$TMP_DIR/assigned_parsed.tsv"

cut -f1 "$TMP_DIR/assigned_parsed.tsv" | sort -u > "$TMP_DIR/subfams.txt"

printf "subfamily\tside\tboundary_bp\tstatus\tfinal_identity_pct\tbackground_identity_pct\n" > "$OUT_TSV"

total_subfams=0
confirmed_count=0
undetermined_count=0
insufficient_count=0

while IFS= read -r subfam; do
  total_subfams=$((total_subfams + 1))
  awk -v sub="$subfam" '$1 == sub' "$TMP_DIR/assigned_parsed.tsv" > "$TMP_DIR/members.tsv"
  count=$(wc -l < "$TMP_DIR/members.tsv")
  if (( count < 20 )); then
    printf "%s\tupstream\t0\tinsufficient_data\tNA\t%s\n" "$subfam" "$BG_IDENTITY" >> "$OUT_TSV"
    printf "%s\tdownstream\t0\tinsufficient_data\tNA\t%s\n" "$subfam" "$BG_IDENTITY" >> "$OUT_TSV"
    insufficient_count=$((insufficient_count + 2))
    continue
  fi
  if (( count > 2000 )); then
    python3 -c "
import random
random.seed(42)
with open('$TMP_DIR/members.tsv') as f:
    lines = f.readlines()
if len(lines) > 2000:
    lines = random.sample(lines, 2000)
with open('$TMP_DIR/members.tsv', 'w') as f:
    f.writelines(lines)
"
  fi

  for side in upstream downstream; do
    ext=$FLANK_BASE
    last_identity="NA"
    status="undetermined"
    boundary_bp=$MAX_EXT_BP
    while (( ext <= MAX_EXT_BP )); do
      awk -v fai="$FAI" -v side="$side" -v ext="$ext" -v STEP_BP="$STEP_BP" '
      BEGIN {
        while ((getline line < fai) > 0) {
          split(line, f, "\t")
          contig_len[f[1]] = f[2]
        }
        FS="\t"
      }
      {
        ctg=$2; start=$3; end=$4; strand=$5
        len = contig_len[ctg]
        if (len == "") len = end
        if (side == "upstream") {
          if (strand ~ /^-/) {
            w_start = end + ext
            w_end = w_start + STEP_BP
          } else {
            w_start = start - ext - STEP_BP
            w_end = start - ext
          }
        } else {
          if (strand ~ /^-/) {
            w_start = start - ext - STEP_BP
            w_end = start - ext
          } else {
            w_start = end + ext
            w_end = w_start + STEP_BP
          }
        }
        if (w_start < 0) w_start = 0
        if (w_end > len) w_end = len
        if (w_start < w_end) {
          print ctg ":" (w_start+1) "-" w_end
        }
      }
      ' "$TMP_DIR/members.tsv" > "$TMP_DIR/regions.txt"

      if [[ ! -s "$TMP_DIR/regions.txt" ]]; then
        last_identity="0"
        status="confirmed"
        boundary_bp=$ext
        break
      fi

      samtools faidx "$GENOME" -r "$TMP_DIR/regions.txt" > "$TMP_DIR/seqs.fasta"
      identity=$(python3 "$TMP_DIR/calc_identity.py" "$TMP_DIR/seqs.fasta" 150)
      last_identity="$identity"
      if (( $(python3 -c "print(1 if $identity <= $IDENTITY_THRESHOLD else 0)") )); then
        status="confirmed"
        boundary_bp=$ext
        break
      fi
      ext=$((ext + STEP_BP))
    done

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$subfam" "$side" "$boundary_bp" "$status" "$last_identity" "$BG_IDENTITY" >> "$OUT_TSV"
    if [[ "$status" == "confirmed" ]]; then
      confirmed_count=$((confirmed_count + 1))
    else
      undetermined_count=$((undetermined_count + 1))
    fi
  done
done < "$TMP_DIR/subfams.txt"

log "Done."
echo "Boundary refinement summary:"
echo "  Total subfamilies processed: $total_subfams"
echo "  Sides confirmed: $confirmed_count"
echo "  Sides undetermined: $undetermined_count"
echo "  Sides insufficient_data: $insufficient_count"
echo "  Background identity: $BG_IDENTITY%"
echo "  Output: $OUT_TSV"
