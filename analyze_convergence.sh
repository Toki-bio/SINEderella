#!/bin/bash
# analyze_convergence.sh
# Analyzes sine_consensus.sh logs to identify slow convergers and noisy families

set -euo pipefail

LOGDIR="${1:-.}"
OUTFILE="/tmp/convergence_analysis.txt"

echo "=== SINE Convergence Quality Analysis ===" > "$OUTFILE"
echo "Analyzing logs in: $LOGDIR" >> "$OUTFILE"
echo "" >> "$OUTFILE"

# Parse each *_consensus.log file
echo "NAME,SEQS,ITERS,BP,STATUS,AVG_CHANGE" > /tmp/convergence_metrics.csv

for logfile in "$LOGDIR"/*_consensus.log; do
    [ -f "$logfile" ] || continue
    
    name=$(basename "$logfile" _consensus.log)
    
    # Extract metrics
    seqs=$(grep "Input:" "$logfile" | grep -oP '\(\K[0-9]+(?= sequences)' || echo "?")
    iters=$(grep "Converged at iter" "$logfile" | grep -oP 'iter \K[0-9]+' || echo "50+")
    bp=$(grep "Output:" "$logfile" | grep -oP '\(\K[0-9]+(?= bp)' || echo "?")
    
    # Check if it converged or hit max
    if grep -q "Converged" "$logfile"; then
        status="CONVERGED"
    elif grep -q "did not converge" "$logfile"; then
        status="NOT_CONVERGED"
    else
        status="UNKNOWN"
    fi
    
    # Calculate average per-iteration change (rough proxy for variance)
    avg_change=$(awk '/Iter [0-9]+: change=/{sum+=$NF; count++} END {if(count>0) print sum/count; else print "N/A"}' "$logfile" || echo "N/A")
    
    echo "$name,$seqs,$iters,$bp,$status,$avg_change" >> /tmp/convergence_metrics.csv
    
done

# Print analysis
echo "=== Convergence Summary ===" >> "$OUTFILE"
cat /tmp/convergence_metrics.csv >> "$OUTFILE"

echo "" >> "$OUTFILE"
echo "=== Interpretation ===" >> "$OUTFILE"
echo "ITERS > 30: Slow convergence (likely mixed/FP data)" >> "$OUTFILE"
echo "NOT_CONVERGED: Hit max 50 iters without convergence (definite problem)" >> "$OUTFILE"
echo "AVG_CHANGE > 0.01: High iteration-to-iteration variance (noisy)" >> "$OUTFILE"

echo "" >> "$OUTFILE"
echo "=== PROBLEMATIC SINES ===" >> "$OUTFILE"
echo "HIGH_ITERS (potential garbage):" >> "$OUTFILE"
awk -F, 'NR>1 && $3+0 > 30 {print "  " $1 " (" $3 " iters, " $2 " seqs)"}' /tmp/convergence_metrics.csv >> "$OUTFILE"

echo "" >> "$OUTFILE"
echo "NOT_CONVERGED (definite problems):" >> "$OUTFILE"
awk -F, 'NR>1 && $5 == "NOT_CONVERGED" {print "  " $1 " (" $2 " seqs)"}' /tmp/convergence_metrics.csv >> "$OUTFILE"

cat "$OUTFILE"
