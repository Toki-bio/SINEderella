#!/usr/bin/env python3
"""pair_aln.py GENOME COPIES.bed CONSENSUS.fa BODY_TSV N MINBODY MAXBODY OUTDIR PREFIX [--mintier T]
Pick N pairs from BODY_TSV (pair_body_id.py output) with MINBODY <= body_identity <= MAXBODY (seeded order), extract
each copy with 150 bp flanks in element orientation, align with the consensus (mafft --localpair), write one FASTA per pair.
The flanks are lower case where the genome is soft-masked; row names carry copy, locus and the two identities."""
import random, subprocess, sys, os
G, BED, CONS, BT, N, LO, HI, OUT, PRE = sys.argv[1:10]; N, LO, HI = int(N), float(LO), float(HI)
loc = {l.split("\t")[3]: l.rstrip("\n").split("\t") for l in open(BED)}
rows = [l.rstrip("\n").split("\t") for i, l in enumerate(open(BT)) if i]
seen, pick = set(), []
for r in sorted(rows, key=lambda r: (r[0], r[1])):
    k = (r[0], r[1])
    if k in seen: continue
    seen.add(k)
    if LO <= float(r[6]) <= HI: pick.append(r)
random.Random(1).shuffle(pick); pick = pick[:N]
os.makedirs(OUT, exist_ok=True)
cons = open(CONS).read()
for i, r in enumerate(pick, 1):
    fa = cons
    for c in (r[0], r[1]):
        ch, s, e, _, _, st = loc[c]
        reg = f"{ch}:{max(1, int(s) - 150 + 1)}-{int(e) + 150}"
        seq = subprocess.run(["samtools", "faidx", G, reg] + (["-i"] if st == "-" else []), capture_output=True, text=True).stdout
        body = "".join(seq.splitlines()[1:])
        fa += f">{c}_{ch}_{s}_{st}_flank{r[3]}_body{r[6]}\n{body}\n"
    tmp = os.path.join(OUT, "tmp.fa"); open(tmp, "w").write(fa)
    aln = subprocess.run(["mafft", "--quiet", "--localpair", "--maxiterate", "1000", "--preservecase", tmp], capture_output=True, text=True).stdout
    open(os.path.join(OUT, f"{PRE}_{i}.aln.fasta"), "w").write(aln); os.remove(tmp)
    print(f"{PRE}_{i}\t{r[0]}\t{r[1]}\tflank {r[3]}\tbody {r[6]}")
