#!/bin/bash
set -euo pipefail

# --- CONFIGURATION ---
CONSENSUS="$1"
GENOME="$2"
SAMPLE_SIZE="${3:-30000}"
BINSIZE="${4:-50}"

# Validate SAMPLE_SIZE as a positive integer
if ! [[ "$SAMPLE_SIZE" =~ ^[0-9]+$ ]] || [[ "$SAMPLE_SIZE" -le 0 ]]; then
    echo "ERROR: SAMPLE_SIZE must be a positive integer (got: '$SAMPLE_SIZE')"
    exit 1
fi

# --- VALIDATION ---
[[ -z "${CONSENSUS:-}" || -z "${GENOME:-}" ]] && {
    echo "Usage: $0 <consensus.fa> <genome.fa> [sample_size] [bin_size]"
    exit 1
}

# Resolve paths for tools, but use original names for OUTDIR
GENOME_ORIGINAL="$GENOME"  # Save original name before resolution
CONSENSUS="$(readlink -f "$CONSENSUS")"
GENOME="$(readlink -f "$GENOME")"

[[ ! -f "$CONSENSUS" ]] && { echo "ERROR: Consensus not found: $CONSENSUS"; exit 1; }
[[ ! -f "$GENOME" ]] && { echo "ERROR: Genome not found: $GENOME"; exit 1; }

# Masked search (docs/MASKING.md): SEARCH_GENOME is a copy of GENOME with known families as N,
# same headers and lengths. sear searches it; extraction below stays on GENOME, so copies and
# their flanks are real bases and coordinates need no mapping.
SEARCH_GENOME="${SEARCH_GENOME:-$GENOME}"
SEARCH_GENOME="$(readlink -f "$SEARCH_GENOME")"
[[ ! -f "$SEARCH_GENOME" ]] && { echo "ERROR: SEARCH_GENOME not found: $SEARCH_GENOME"; exit 1; }
if [[ "$SEARCH_GENOME" != "$GENOME" ]]; then
    echo "[$(date)] Masked search: searching $SEARCH_GENOME, extracting from $GENOME"
fi

# Check required tools
for tool in samtools seqkit sear bedtools SubFam mafft; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: Required tool '$tool' not found in PATH"; exit 1;
    }
done

# --- SETUP ---
shopt -s nullglob

if [[ ! -f "${GENOME}.fai" ]]; then
    echo "[$(date)] Indexing genome..."
    samtools faidx "$GENOME"
fi

# Use original name for output directory (not resolved symlink)
GENOME_BASENAME=$(basename "$GENOME_ORIGINAL")
OUTDIR="${GENOME_BASENAME%.*}_step1"
mkdir -p "$OUTDIR/searches"
echo "[$(date)] Output directory: $OUTDIR"

# Save metadata for reproducibility
{
    echo "# SINEderella Step1 Run"
    echo "Date: $(date)"
    echo "Consensus: $CONSENSUS ($(grep -c "^>" "$CONSENSUS" 2>/dev/null || echo "?") sequences)"
    echo "Genome: $GENOME_ORIGINAL"
    echo "Sample size: $SAMPLE_SIZE"
    echo "Bin size: $BINSIZE"
    echo "Command: $0 $*"
    echo "Hostname: $(hostname)"
    echo "PID: $$"
} > "$OUTDIR/run_info.txt"

# Save paths for step2 (use original names to avoid confusion)
echo "$GENOME_ORIGINAL" > "$OUTDIR/genome.path"
echo "$CONSENSUS" > "$OUTDIR/consensus.path"

cd "$OUTDIR"

# --- SPLIT CONSENSUS ---
echo "[$(date)] Splitting consensus library..."
seqkit split -i "$CONSENSUS" -O searches
# seqkit < 0.15 names the per-sequence files <name>.id_<ID>.<ext>; newer versions use .part_<ID>. Normalise to .part_
for _f in searches/*.id_*; do
    [[ -e "$_f" ]] || continue
    _d=$(dirname "$_f"); _b=$(basename "$_f")
    mv "$_f" "$_d/${_b/.id_/.part_}"
done

CONSENSUS_BASENAME=$(basename "$CONSENSUS")
CONSENSUS_NAME="${CONSENSUS_BASENAME%.*}"
CONSENSUS_EXT="${CONSENSUS_BASENAME##*.}"

# --- SEARCH ---
cd searches
echo "[$(date)] Starting searches..."

QUERY_FILES=("${CONSENSUS_NAME}".part_*."${CONSENSUS_EXT}")

if [[ ${#QUERY_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No split files found matching pattern: ${CONSENSUS_NAME}.part_*.${CONSENSUS_EXT}"
    exit 1
fi

echo "[$(date)] Found ${#QUERY_FILES[@]} consensus sequences to search"

# KEY CHANGE: keep genome split parts during this entire Step1 run (reuse across queries)
# sear must implement -k/--keep as discussed
for query in "${QUERY_FILES[@]}"; do
    echo "[$(date)] Searching with $(basename "$query")"
    sear -k "$query" "$SEARCH_GENOME" 0.8 65 50

    # Report sear results to stderr (visible in SINEderella_multi terminal)
    _qbase="${query%.*}"
    _sfname="${_qbase##*.part_}"
    _bed="gen-${_qbase}.bed"
    if [[ -f "$_bed" ]]; then
        _nhits=$(wc -l < "$_bed")
        echo "[$(date)]   -> ${_sfname}: ${_nhits} hits" >&2
    else
        echo "[$(date)]   -> ${_sfname}: 0 hits (no BED output)" >&2
    fi
done

# KEY CHANGE: delete genome split parts ONCE, after last query
echo "[$(date)] Cleaning genome split parts (end of Step1 searches)..."
rm -f *.2k.bnk *.2k.bnk.fai *.2k.part_* *.2k.part_*.bnk.fai *_s.fai 2>/dev/null || true

cd ..

# --- SATELLITE SCREEN (docs/SATELLITES.md): SINE-derived satellites and SINE-containing tandem arrays, found from the hits
#     and TRF on the hit windows, before extraction / SubFam / assignment see them. SKIP_SATELLITES=1 skips it;
#     SATELLITE_EXCLUDE=0 writes the tables but removes nothing; SATELLITE_EXCLUDE_B=verified|flagged|long|all|none for kind-B
#     arrays (default verified: runs whose units are near-identical, docs/SATELLITES.md 5h; the tool's own default is the same).
_TOOLS="${SINEDERELLA_TOOLS:-$(dirname "$(readlink -f "$0")")/tools}"
if [[ "${SKIP_SATELLITES:-0}" == "1" ]]; then
    echo "[$(date)] SKIP_SATELLITES=1 - satellite screen skipped"
elif [[ -f "$_TOOLS/satellite_stage.py" ]] && command -v trf >/dev/null 2>&1; then
    echo "[$(date)] Satellite screen (tools/satellite_stage.py) on searches/"
    _SAT_ARGS=(--searches searches --genome "$GENOME" --cons "$CONSENSUS" --out satellites --threads "${THREADS:-$(nproc 2>/dev/null || echo 8)}")
    [[ "${SATELLITE_EXCLUDE:-1}" == "0" ]] && _SAT_ARGS+=(--no-exclude)
    _SAT_ARGS+=(--exclude-b "${SATELLITE_EXCLUDE_B:-verified}")
    python3 "$_TOOLS/satellite_stage.py" "${_SAT_ARGS[@]}" 2>&1 | tee satellites.log >&2         || echo "[$(date)] WARNING: satellite screen failed (see satellites.log); continuing with the unfiltered hits" >&2
else
    echo "[$(date)] Satellite screen skipped (tools/satellite_stage.py or trf not found)"
fi

# --- REQUIRE sear output ---
BED_FILES=(searches/*.bed)
if [[ ${#BED_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BED files generated by sear"
    exit 1
fi

# --- BUILD FILES REQUIRED BY postprocess.sh ---
# Creates:
#   searches/all_hits.labeled.bed
#   searches/regions.by_subfam.bed
echo "[$(date)] Building searches/all_hits.labeled.bed and searches/regions.by_subfam.bed..."

(
  cd searches || exit 1

  : > all_hits.labeled.bed

  for bed in *.bed; do
    [[ -s "$bed" ]] || continue

    # bed basename -> query basename
    # sear produces gen-<querybase>.bed in your runs, so strip that
    base="${bed%.bed}"
    base="${base#gen-}"

    # find the corresponding split FASTA file deterministically
    # prefer exact extension match, otherwise first match
    qfile=""
    if [[ -f "${base}.${CONSENSUS_EXT}" ]]; then
      qfile="${base}.${CONSENSUS_EXT}"
    else
      qfile="$(ls -1 "${base}".* 2>/dev/null | head -n1 || true)"
    fi

    # subfamily label: from FASTA header (first token) if possible, else fallback to base
    sf="$base"
    if [[ -n "$qfile" && -f "$qfile" ]]; then
      sf="$(awk '/^>/{h=$1; sub(/^>/,"",h); print h; exit}' "$qfile")"
      [[ -n "$sf" ]] || sf="$base"
    fi

    # Write as 6-col BED: chr start end subfam score strand
    # Be robust to different .bed formats
    awk -v sf="$sf" 'BEGIN{OFS="\t"}
      {
        score = (NF>=5 ? $5 : 1)
        strand = (NF>=6 ? $6 : ".")
        print $1,$2,$3,sf,score,strand
      }' "$bed" >> all_hits.labeled.bed
  done

  [[ -s all_hits.labeled.bed ]] || { echo "ERROR: all_hits.labeled.bed ended up empty"; exit 1; }

  bedtools sort -i all_hits.labeled.bed > all_hits.labeled.bed.tmp
  mv all_hits.labeled.bed.tmp all_hits.labeled.bed

  # Conflict regions: merge all labeled hits and collect distinct subfam labels
  # Output is 4-col: chr start end subfam_or_comma_list
  bedtools merge -i all_hits.labeled.bed -c 4 -o distinct > regions.by_subfam.bed

  [[ -s regions.by_subfam.bed ]] || { echo "ERROR: regions.by_subfam.bed ended up empty"; exit 1; }

  reg_n=$(wc -l < regions.by_subfam.bed)
  conf_n=$(awk -F'\t' '$4 ~ /,/' regions.by_subfam.bed | wc -l)
  echo "[$(date)] regions.by_subfam.bed: $reg_n regions ($conf_n with conflicts)"
)

# --- MERGE AND EXTRACT (your original behavior) ---
echo "[$(date)] Found ${#BED_FILES[@]} BED files. Merging..."
cat "${BED_FILES[@]}" | bedtools sort -i - | \
    bedtools merge -i - -c 4,5,6 -o max,max,distinct > merged_hits.bed

BED_COUNT=$(wc -l < merged_hits.bed)
echo "[$(date)] Merged into $BED_COUNT unique intervals"

# --- SATELLITE REGIONS ON THE MERGED LOCI: the screen above filtered each consensus' own hit list, but a locus inside a verified
#     array can come back through another consensus' hit at the same place (rsi MEG-RS: 6 of 129 remaining copies sat inside
#     excluded spans, 2026-10-05). Remove every merged interval that overlaps an excluded region (kind-A loci, verified arrays).
if [[ "${SATELLITE_EXCLUDE:-1}" != "0" && -s satellites/exclude_regions.bed ]]; then
    cp merged_hits.bed merged_hits.bed.before_satellites
    bedtools intersect -v -a merged_hits.bed -b satellites/exclude_regions.bed > merged_hits.bed.tmp
    mv merged_hits.bed.tmp merged_hits.bed
    _after=$(wc -l < merged_hits.bed)
    echo "[$(date)] Satellite regions: $((BED_COUNT - _after)) of $BED_COUNT merged intervals lie in excluded regions (satellites/exclude_regions.bed) and were removed (kept in merged_hits.bed.before_satellites)"
    BED_COUNT=$_after
fi

echo "[$(date)] Extracting sequences..."
bedtools getfasta -s -fi "$GENOME" -bed merged_hits.bed > extracted.fasta

EXTRACTED=$(grep -c "^>" extracted.fasta)
echo "[$(date)] Extracted $EXTRACTED sequences"
# every merged interval must come back as one sequence (getfasta skips intervals it cannot read)
if [[ "$EXTRACTED" -ne "$BED_COUNT" ]]; then
    echo "ERROR: extracted $EXTRACTED sequences from $BED_COUNT merged intervals" >&2; exit 1
fi

# --- SAMPLING ---
if [[ $EXTRACTED -le $SAMPLE_SIZE ]]; then
    echo "[$(date)] Fewer sequences ($EXTRACTED) than sample size, using all"
    cp extracted.fasta sampled_"${SAMPLE_SIZE}".fasta
else
    echo "[$(date)] Sampling $SAMPLE_SIZE sequences..."
    [[ -f "extracted.fasta" ]] || { echo "ERROR: extracted.fasta not found"; exit 1; }
    seqkit sample -n "$SAMPLE_SIZE" extracted.fasta > sampled_"${SAMPLE_SIZE}".fasta
fi

# --- CLUSTERING ---
echo "[$(date)] Running SubFam..."
mkdir -p subfam_input
cp sampled_"${SAMPLE_SIZE}".fasta subfam_input/input.fasta
cd subfam_input
SubFam input.fasta "$BINSIZE" 2>&1 | tee subfam_log.txt

# --- MAFFT ALIGNMENT ---
CONSENSUS_PATH=$(cat ../consensus.path)
[[ -f "input.clw" ]] || { echo "ERROR: SubFam failed to produce input.clw"; exit 1; }

FORMAT=$(head -n 1 input.clw)
if echo "$FORMAT" | grep -q "^CLUSTAL"; then
    echo "[$(date)] Converting input.clw from Clustal to FASTA..."
    awk '
    /^CLUSTAL/ {next}
    /^$/ {next}
    {
      if (NF >= 2) {
        name = $1
        sequence = $2
        gsub(/[-.]/, "", sequence)
        if (!(name in seq)) {
          order[++count] = name
        }
        seq[name] = seq[name] sequence
      }
    }
    END {
      for (i=1; i<=count; i++) {
        print ">" order[i]
        print seq[order[i]]
      }
    }' input.clw > input_reps.fasta
else
    echo "[$(date)] input.clw appears to be in FASTA format, copying..."
    cp input.clw input_reps.fasta
fi

CLUSTERS=$(grep -c "^>" input_reps.fasta || echo 0)
echo "[$(date)] Aligning $CLUSTERS clusters with consensus using MAFFT..."

cat input_reps.fasta "$CONSENSUS_PATH" > combined_input.fasta

mafft --thread "${THREADS:-$(nproc)}" --threadit 0 \
      --localpair \
      --maxiterate 1000 \
      --ep 0.123 \
      --nuc \
      --reorder \
      --preservecase \
      --quiet \
      combined_input.fasta > input.clw.al

[[ -s "input.clw.al" ]] || { echo "ERROR: MAFFT alignment failed or produced empty output"; exit 1; }

echo "[$(date)] MAFFT completed successfully"
echo "[$(date)] Final output: $(pwd)/input.clw.al"
cd ..

echo "Number of hits per query:"
for bed in searches/gen-*.bed; do
    if [[ -f "$bed" ]]; then
        count=$(wc -l < "$bed")
        # Extract subfamily name from gen-consensuses.clean.part_<NAME>.bed
        _bn="$(basename "$bed" .bed)"
        _sfname="${_bn#gen-consensuses.clean.part_}"
        [[ "$_sfname" == "$_bn" ]] && _sfname="$_bn"
        echo "  $_sfname: $count hits"
    fi
done
echo "Total merged hits: $BED_COUNT"
