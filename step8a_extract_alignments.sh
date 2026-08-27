#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# step8a_extract_alignments.sh
#
# Build per-subfamily alignment files from a completed SINEderella run.
#
# Usage: step8a_extract_alignments.sh <RUN_ROOT> <SPECIES_CODE>
###############################################################################

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <RUN_ROOT> <SPECIES_CODE>" >&2
    exit 1
fi

RUN_ROOT="$(cd -- "$1" && pwd -P)"
SPECIES_CODE="$2"

# Validate inputs
[[ -d "$RUN_ROOT" ]] || { echo "ERROR: RUN_ROOT not a directory: $RUN_ROOT" >&2; exit 1; }

GENOME="$RUN_ROOT/genome.clean.fa"
CONSENSUS="$RUN_ROOT/consensuses.clean.fa"
[[ -s "$GENOME" ]]    || { echo "ERROR: Missing genome.clean.fa" >&2; exit 1; }
[[ -s "$CONSENSUS" ]] || { echo "ERROR: Missing consensuses.clean.fa" >&2; exit 1; }

# Find step2 output directory
STEP2_OUT="$(ls -dt "$RUN_ROOT"/step2/step2_output* 2>/dev/null | head -n1 || true)"
[[ -n "$STEP2_OUT" && -d "$STEP2_OUT" ]] || { echo "ERROR: Cannot find step2 output directory" >&2; exit 1; }
ASSIGNED="$STEP2_OUT/assigned.fasta"
[[ -s "$ASSIGNED" ]] || { echo "ERROR: Missing assigned.fasta in $STEP2_OUT" >&2; exit 1; }

# Boundary refinement file (optional)
BOUNDARY_TSV="$STEP2_OUT/boundary_refinement.tsv"

# Output directory
OUTDIR="$RUN_ROOT/results/alignments"
mkdir -p "$OUTDIR"

# Temp working directory under RUN_ROOT
TMPDIR="$(mktemp -d "${RUN_ROOT}/.step8a_XXXXXX")"
trap 'rm -rf "$TMPDIR"' EXIT

# Ensure genome is indexed for bedtools
if [[ ! -s "$GENOME.fai" ]]; then
    samtools faidx "$GENOME"
fi

# MAFFT parameters (shared across all variants)
MAFFT_ARGS=(--localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase)

# Base flank sizes
BASE_UP=50
BASE_DOWN=70

###############################################################################
# --reorder above reorders MAFFT's OUTPUT by guide tree similarity, so putting
# the consensus first in the INPUT does not guarantee it stays first in the
# OUTPUT. Reorder the aligned output afterward so the consensus record is
# always first, regardless of where mafft placed it.
###############################################################################
cat > "$TMPDIR/reorder_consensus_first.py" <<'PYEOF'
import sys
cons_header = None
with open(sys.argv[2]) as f:
    for line in f:
        if line.startswith('>'):
            cons_header = line[1:].strip()
            break
recs = []
name = None
seq = []
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('>'):
            if name is not None:
                recs.append((name, ''.join(seq)))
            name = line[1:].strip()
            seq = []
        else:
            seq.append(line.rstrip('\n'))
    if name is not None:
        recs.append((name, ''.join(seq)))
cons_idx = None
for i, (n, s) in enumerate(recs):
    if n == cons_header:
        cons_idx = i
        break
if cons_idx is not None and cons_idx != 0:
    recs.insert(0, recs.pop(cons_idx))
for n, s in recs:
    sys.stdout.write('>' + n + '\n' + s + '\n')
PYEOF

reorder_consensus_first() {
    # args: aligned_fasta cons_fasta -> writes reordered fasta to stdout
    python3 "$TMPDIR/reorder_consensus_first.py" "$1" "$2"
}

###############################################################################
# Parse boundary_refinement.tsv (optional)
# Columns: subfamily, population, side, boundary_bp, status, ...
#
# A boundary measured on a random/general sample is NOT valid for the
# top100 variant: top100's bitscore-ranked members are the least-diverged
# copies, not a random draw, and can stay non-independent far past where
# the general population reaches background (found live against scorpion
# g01 -- see step7_boundary_refine.sh's comment on this same issue). So
# top100 uses the "top100"-population rows; rand100 and subfam use the
# "general"-population rows.
# Output per population file: subfamily \t upstream_ext \t downstream_ext
###############################################################################
build_boundaries_file(){
    local population="$1" outfile="$2"
    if [[ -s "$BOUNDARY_TSV" ]]; then
        # Only use boundary_bp when status=="confirmed". "undetermined" means
        # the test hit MAX_EXT_BP without ever confirming the flank reached
        # background -- it is NOT a validated extension, just wherever the
        # loop ran out of room. Applying it anyway (found live: g17
        # downstream got a 1070bp flank -- base 70 + the full 1000bp cap --
        # for a side that was explicitly never confirmed unique, and made
        # mafft --localpair --maxiterate 1000 extremely slow on ~100
        # sequences that long) asserts confidence the test never earned.
        # Fall back to ext=0 (base flank only) on anything not confirmed.
        awk -F'\t' -v pop="$population" 'NR>1 && $2==pop && $5=="confirmed" {
            if ($3 == "upstream") up[$1] = $4 + 0
            else if ($3 == "downstream") down[$1] = $4 + 0
        }
        END {
            for (sf in up) {
                d = (sf in down) ? down[sf] : 0
                print sf "\t" up[sf] "\t" d
            }
            for (sf in down) {
                if (!(sf in up)) print sf "\t" 0 "\t" down[sf]
            }
        }' "$BOUNDARY_TSV" > "$outfile"
    else
        : > "$outfile"
    fi
}
build_boundaries_file "general" "$TMPDIR/boundaries_general.tsv"
build_boundaries_file "top100" "$TMPDIR/boundaries_top100.tsv"

###############################################################################
# Parse assigned.fasta -> loci.tsv
# Header format: >ctg:start-end(strand)|subfamily|bitscore
# Output TSV: subfamily, bitscore, ctg, start, end, strand
###############################################################################
awk 'BEGIN{OFS="\t"}
/^>/ {
    line = substr($0, 2)
    gsub(/\r/, "", line)
    n = split(line, parts, "|")
    loc = parts[1]
    subfam = parts[2]
    bitscore = parts[3] + 0

    # Extract strand from parentheses
    strand = loc
    sub(/.*\(/, "", strand)
    sub(/\)/, "", strand)

    # Remove (strand) suffix
    sub(/\([^)]*\)$/, "", loc)

    # Split ctg:start-end
    ctg = loc
    sub(/:[0-9]+-[0-9]+$/, "", ctg)
    coords = loc
    sub(/.*:/, "", coords)
    split(coords, se, "-")
    start = se[1]
    end = se[2]

    # Sanitize strand for BED (handle "+,-" etc.)
    if (strand != "+" && strand != "-") strand = "+"

    print subfam, bitscore, ctg, start, end, strand
}' "$ASSIGNED" > "$TMPDIR/loci.tsv"

###############################################################################
# Subfamily list and member counts
###############################################################################
cut -f1 "$TMPDIR/loci.tsv" | sort | uniq -c | awk '{print $2 "\t" $1}' > "$TMPDIR/counts.tsv"

###############################################################################
# Helper: extract consensus for a subfamily from consensuses.clean.fa
###############################################################################
extract_consensus() {
    local name="$1" outfile="$2"
    awk -v name="$name" '
        BEGIN{want=0}
        /^>/{
            hdr=substr($0,2)
            sub(/[ \t].*$/,"",hdr)
            want=(hdr==name)
        }
        want{print}
    ' "$CONSENSUS" > "$outfile"
}

###############################################################################
# Helper: create BED (0-based) from loci TSV, extract with flanks, align
# Args: loci_tsv, cons_fa, up_flank, down_flank, outfile
###############################################################################
extract_flank_align() {
    local loci_tsv="$1" cons_fa="$2" up_flank="$3" down_flank="$4" outfile="$5"

    # Create BED: chrom, start-1, end, name, score, strand
    awk -F'\t' 'BEGIN{OFS="\t"} {s=$4-1; if (s<0) s=0; print $3, s, $5, NR, $2, $6}' \
        "$loci_tsv" > "$TMPDIR/cur.bed" || return 1

    # Strand-aware slop
    bedtools slop -i "$TMPDIR/cur.bed" -g "$GENOME.fai" \
        -l "$up_flank" -r "$down_flank" -s > "$TMPDIR/cur_slop.bed" || return 1

    # Extract sequences (strand-aware)
    bedtools getfasta -fi "$GENOME" -bed "$TMPDIR/cur_slop.bed" -s \
        > "$TMPDIR/cur_extracted.fa" 2>/dev/null || return 1

    # Append consensus
    cat "$cons_fa" "$TMPDIR/cur_extracted.fa" > "$TMPDIR/cur_combined.fa" || return 1

    # Align
    mafft "${MAFFT_ARGS[@]}" "$TMPDIR/cur_combined.fa" \
        > "$TMPDIR/cur_aligned.fa" 2>"$TMPDIR/cur_mafft.err" || return 1

    # --reorder can move the consensus away from the top; force it back first
    reorder_consensus_first "$TMPDIR/cur_aligned.fa" "$cons_fa" > "$TMPDIR/cur_reordered.fa" || return 1

    # Restore @U@ -> _ in headers and write output
    sed '/^>/s/@U@/_/g' "$TMPDIR/cur_reordered.fa" > "$outfile" || return 1
}

###############################################################################
# Process each subfamily
###############################################################################
MANIFEST="$TMPDIR/manifest.tsv"
printf "subfamily\thas_top100\thas_rand100\thas_subfam\tn_members\n" > "$MANIFEST"

idx=0
while IFS=$'\t' read -r subfam count; do
    [[ -n "$subfam" ]] || continue
    idx=$((idx + 1))

    echo "[$(date '+%F %T')] Processing subfamily: $subfam ($count members)" >&2

    # Extract this subfamily's loci
    awk -F'\t' -v sf="$subfam" '$1==sf' "$TMPDIR/loci.tsv" > "$TMPDIR/loci_${idx}.tsv"

    # Determine flank sizes from boundary_refinement.tsv -- top100 and
    # rand100/subfam use DIFFERENT boundary populations (see the comment
    # above build_boundaries_file: a general/random-sample boundary is not
    # valid for the bitscore-ranked top100 set).
    general_up_ext=0; general_down_ext=0
    boundary_line="$(awk -F'\t' -v sf="$subfam" '$1==sf{print $2"\t"$3; exit}' "$TMPDIR/boundaries_general.tsv" || true)"
    if [[ -n "$boundary_line" ]]; then
        general_up_ext="$(printf '%s' "$boundary_line" | cut -f1)"
        general_down_ext="$(printf '%s' "$boundary_line" | cut -f2)"
    fi
    general_up_flank=$((BASE_UP + general_up_ext))
    general_down_flank=$((BASE_DOWN + general_down_ext))

    top100_up_ext=0; top100_down_ext=0
    boundary_line="$(awk -F'\t' -v sf="$subfam" '$1==sf{print $2"\t"$3; exit}' "$TMPDIR/boundaries_top100.tsv" || true)"
    if [[ -n "$boundary_line" ]]; then
        top100_up_ext="$(printf '%s' "$boundary_line" | cut -f1)"
        top100_down_ext="$(printf '%s' "$boundary_line" | cut -f2)"
    fi
    top100_up_flank=$((BASE_UP + top100_up_ext))
    top100_down_flank=$((BASE_DOWN + top100_down_ext))

    # Extract consensus for this subfamily
    extract_consensus "$subfam" "$TMPDIR/cons_${idx}.fa"
    if [[ ! -s "$TMPDIR/cons_${idx}.fa" ]]; then
        echo "WARNING: No consensus found for $subfam -- skipping" >&2
        printf "%s\t0\t0\t0\t%s\n" "$subfam" "$count" >> "$MANIFEST"
        continue
    fi

    has_top=0
    has_rand=0
    has_sub=0

    # -- top100: 100 highest-bitscore members --
    # Sort by bitscore descending -> temp file, then head (avoids SIGPIPE)
    sort -t$'\t' -k2,2nr "$TMPDIR/loci_${idx}.tsv" > "$TMPDIR/sorted_${idx}.tsv"
    head -100 "$TMPDIR/sorted_${idx}.tsv" > "$TMPDIR/top100_${idx}.tsv"

    if [[ -s "$TMPDIR/top100_${idx}.tsv" ]]; then
        set +e
        extract_flank_align "$TMPDIR/top100_${idx}.tsv" "$TMPDIR/cons_${idx}.fa" \
            "$top100_up_flank" "$top100_down_flank" \
            "$OUTDIR/${SPECIES_CODE}_${subfam}_top100.aln.fa"
        rc=$?
        set -e
        if (( rc == 0 )); then
            has_top=1
        else
            echo "WARNING: top100 alignment failed for $subfam (rc=$rc)" >&2
            tail -n 20 "$TMPDIR/cur_mafft.err" >&2 || true
        fi
    fi

    # -- rand100: 100 randomly sampled members --
    shuf "$TMPDIR/loci_${idx}.tsv" > "$TMPDIR/shuffled_${idx}.tsv"
    head -100 "$TMPDIR/shuffled_${idx}.tsv" > "$TMPDIR/rand100_${idx}.tsv"

    if [[ -s "$TMPDIR/rand100_${idx}.tsv" ]]; then
        set +e
        extract_flank_align "$TMPDIR/rand100_${idx}.tsv" "$TMPDIR/cons_${idx}.fa" \
            "$general_up_flank" "$general_down_flank" \
            "$OUTDIR/${SPECIES_CODE}_${subfam}_rand100.aln.fa"
        rc=$?
        set -e
        if (( rc == 0 )); then
            has_rand=1
        else
            echo "WARNING: rand100 alignment failed for $subfam (rc=$rc)" >&2
            tail -n 20 "$TMPDIR/cur_mafft.err" >&2 || true
        fi
    fi

    # -- subfam: only for subfamilies with >= 400 members --
    if (( count >= 400 )); then
        # Fixed-seed shuffle (seed=42) -> take up to 10000
        # Write to temp file first to avoid SIGPIPE with head
        awk -F'\t' 'BEGIN{srand(42); OFS="\t"} {print rand(), $0}' \
            "$TMPDIR/loci_${idx}.tsv" | \
            sort -t$'\t' -k1,1g | cut -f2- > "$TMPDIR/seeded_${idx}.tsv"
        head -10000 "$TMPDIR/seeded_${idx}.tsv" > "$TMPDIR/subfam_sample_${idx}.tsv"

        if [[ -s "$TMPDIR/subfam_sample_${idx}.tsv" ]]; then
            # Extract element sequences (no flanks)
            awk -F'\t' 'BEGIN{OFS="\t"} {s=$4-1; if (s<0) s=0; print $3, s, $5, NR, $2, $6}' \
                "$TMPDIR/subfam_sample_${idx}.tsv" > "$TMPDIR/subfam_${idx}.bed"

            # NOT wrapping this in set +e/-e caused a real bug (found live
            # against scorpion g02, 30284 members): bedtools getfasta failing
            # here killed the whole script silently under set -e, with the
            # real error swallowed by 2>/dev/null -- reproduced deterministically
            # with bash -x tracing, which is the only reason this was caught.
            set +e
            bedtools getfasta -fi "$GENOME" -bed "$TMPDIR/subfam_${idx}.bed" -s \
                > "$TMPDIR/subfam_elems_${idx}.fa" 2>"$TMPDIR/subfam_getfasta_${idx}.err"
            getfasta_rc=$?
            set -e

            subfam_ran=0
            if (( getfasta_rc != 0 )) || [[ ! -s "$TMPDIR/subfam_elems_${idx}.fa" ]]; then
                echo "WARNING: bedtools getfasta failed for subfam variant of $subfam (rc=$getfasta_rc) -- skipping subfam variant" >&2
                tail -n 20 "$TMPDIR/subfam_getfasta_${idx}.err" >&2 || true
                subfam_rc=1
            else
            subfam_ran=1

            # Run SubFam in scratch directory
            scratch="$TMPDIR/subfam_scratch_${idx}"
            mkdir -p "$scratch"
            cp "$TMPDIR/subfam_elems_${idx}.fa" "$scratch/input.fasta"

            set +e
            (cd "$scratch" && SubFam input.fasta 50)
            subfam_rc=$?
            set -e
            fi

            if (( subfam_rc == 0 )) && [[ -s "$scratch/input.clw" ]]; then
                # Degap the .clw file (FASTA-like, gapped per-chunk consensuses)
                awk '!/^>/{gsub(/-/, "")}1' "$scratch/input.clw" \
                    > "$TMPDIR/subfam_degapped_${idx}.fa"

                # Align degapped subfam consensuses with the subfamily's own consensus
                # (consensus first, so it ends up on top of the alignment)
                cat "$TMPDIR/cons_${idx}.fa" "$TMPDIR/subfam_degapped_${idx}.fa" \
                    > "$TMPDIR/subfam_combined_${idx}.fa"

                set +e
                mafft "${MAFFT_ARGS[@]}" "$TMPDIR/subfam_combined_${idx}.fa" \
                    > "$TMPDIR/subfam_aligned_${idx}.fa" 2>"$TMPDIR/subfam_mafft_${idx}.err"
                mafft_rc=$?
                set -e

                if (( mafft_rc == 0 )); then
                    reorder_consensus_first "$TMPDIR/subfam_aligned_${idx}.fa" "$TMPDIR/cons_${idx}.fa" \
                        > "$TMPDIR/subfam_reordered_${idx}.fa"
                    sed '/^>/s/@U@/_/g' "$TMPDIR/subfam_reordered_${idx}.fa" \
                        > "$OUTDIR/${SPECIES_CODE}_${subfam}_subfam.aln.fa"
                    has_sub=1
                else
                    echo "WARNING: subfam mafft failed for $subfam (rc=$mafft_rc)" >&2
                    tail -n 20 "$TMPDIR/subfam_mafft_${idx}.err" >&2 || true
                fi
            elif (( subfam_ran == 1 )); then
                echo "WARNING: SubFam failed for $subfam (rc=$subfam_rc) -- skipping subfam variant" >&2
            fi
        fi
    fi

    # Manifest entry
    printf "%s\t%s\t%s\t%s\t%s\n" "$subfam" "$has_top" "$has_rand" "$has_sub" "$count" >> "$MANIFEST"

done < "$TMPDIR/counts.tsv"

# Copy manifest to output
cp "$MANIFEST" "$OUTDIR/manifest.tsv"

echo "[$(date '+%F %T')] Done. Alignments in $OUTDIR, manifest in $OUTDIR/manifest.tsv" >&2
