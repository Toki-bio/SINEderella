# What SINEderella is for (owner's words, 2026-10-17)

This file is the reference for every change. A change that does not serve this aim does not belong in the workflow.

## Subject
We study **interspersed repeats**: elements whose copies have **unique flanks**. Copies whose flanks are shared between
members (segmental duplications, array units, insertion hot spots, anything multiplied with its neighbourhood) are not
the subject. They are never shown in alignments and never asked about; at most they are counted.

## Outputs
SINEderella must bring:
1. **The technical result:** the list of SINE loci in a genome.
2. **A user-readable output:** an HTML page with
   - links to alignments,
   - conclusions on alignment parts,
   - calculations such as copy numbers, leak, and divergence from consensus.

## Requirements
The workflow must be logical, effective, reproducible and reliable. Everything above is produced by the workflow
itself, on every run, without manual checking.
