# Real-data test of the fs9 twin rules: Squam1 and Squam3C (2026-10-14)

Two squamate SINEs that multiply into arrays and satellites, in soft-masked assemblies:

| | genome | soft-masked | full-length copies |
|---|---|---|---|
| Squam1 | *Anolis carolinensis*, UCSC anoCar2 | 36 % | 79,225 |
| Squam3C | *Python bivittatus*, NCBI RefSeq GCF_000186305.1 | 24 % | 6,713 |

Consensuses: `Squam1.fa` and `Squam3C.fa`, from Vassetzky, Kosushkin & Ryskov 2023, *Mobile DNA* 14:21, Additional files 9 and 2.
Copies: blastn hits covering at least 80 % of the consensus at 65 % identity or more, overlapping hits of one strand merged.
Run: `run_real_twins.sh WORKDIR [THREADS]` (downloads both genomes; blastn, bedtools, samtools, jellyfish, ssearch36, mafft).
The body test of a pair aligns the two SINE bodies (ssearch36, best HSP); the background is 2,000 random copy pairs of the family (seed 1).

## Results

Share of pairs whose SINE bodies are more alike than the 95th percentile of random copy pairs:

| pairs | Anolis Squam1 | Python Squam3C |
|---|---|---|
| strict flank test (masked bases as N, SOFTMASK=1) | 87 % of 1,875 | 79 % of 57 |
| added when masked bases are kept (SOFTMASK=0) | 18 % of 151,502 (sample of 2,000) | 4 % of 870 |
| added by the masked 95 % pass of SOFTMASK=both | not measured | 8 % of 38 |

Twin copies after the body check (fs9_bodycheck.py, partner-corrected threshold):

| | Anolis Squam1 | Python Squam3C |
|---|---|---|
| SOFTMASK=1 + body check | 1,952 | 37 |
| SOFTMASK=0 + body check (default) | 5,358 | 60 |
| SOFTMASK=0, flank test alone | 17,068 | 174 |

The extra twins of SOFTMASK=0 + body check have body identity that tracks flank identity (median difference -0.1 in Anolis,
+0.7 in Python), as copies duplicated together with their neighbourhood must; hot-spot pairs have flanks of 85-98 % and bodies at
the family background (83-85 %).

## Alignments to look at

Each file: the consensus, then the two copies with 150 bp of flank on each side in element orientation, lower case where the genome is
soft-masked. Row names end in the flank and the body identity of the pair. Look at whether the similarity of the two copies runs
through the element at the same level as through the flanks (duplication) or drops to background inside the element (hot spot).

| what | Anolis Squam1 | Python Squam3C |
|---|---|---|
| twin, kept (flanks and bodies alike) | [aca_Squam1_twin_1](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/aca_Squam1_twin_1.aln.fasta&title=aca_Squam1_twin_1), [aca_Squam1_twin_2](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/aca_Squam1_twin_2.aln.fasta&title=aca_Squam1_twin_2) | [pyt_Squam3C_twin_1](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/pyt_Squam3C_twin_1.aln.fasta&title=pyt_Squam3C_twin_1), [pyt_Squam3C_twin_2](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/pyt_Squam3C_twin_2.aln.fasta&title=pyt_Squam3C_twin_2) |
| hot spot, removed (flanks alike, bodies at background) | [aca_Squam1_hotspot_1](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/aca_Squam1_hotspot_1.aln.fasta&title=aca_Squam1_hotspot_1), [aca_Squam1_hotspot_2](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/aca_Squam1_hotspot_2.aln.fasta&title=aca_Squam1_hotspot_2) | [pyt_Squam3C_hotspot_1](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/pyt_Squam3C_hotspot_1.aln.fasta&title=pyt_Squam3C_hotspot_1), [pyt_Squam3C_hotspot_2](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/pyt_Squam3C_hotspot_2.aln.fasta&title=pyt_Squam3C_hotspot_2) |
| added by the masked 95 % pass, removed by the body check | | [pyt_Squam3C_masked95_1](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/pyt_Squam3C_masked95_1.aln.fasta&title=pyt_Squam3C_masked95_1), [pyt_Squam3C_masked95_2](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/pyt_Squam3C_masked95_2.aln.fasta&title=pyt_Squam3C_masked95_2) |
| tandem array unit (1.8 kb apart) | [aca_Squam1_array_1](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINEderella/8c7835829cf760edd52f6cb6367f428245186fa9/flankscan/tests/real_twins/aca_Squam1_array_1.aln.fasta&title=aca_Squam1_array_1) | |
