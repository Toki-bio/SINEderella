#!/usr/bin/env python3
"""random_pairs.py COPIES.bed N SEED: N random copy pairs (seeded) in flank_groups.tsv layout, for the body background"""
import random, sys
cop = [l.split("\t")[3] for l in open(sys.argv[1])]; rnd = random.Random(int(sys.argv[3]))
print("copy_a\tcopy_b\tside\tidentity\tlength\ttier\tclass")
for _ in range(int(sys.argv[2])):
    a, b = rnd.sample(cop, 2); print(f"{a}\t{b}\t5\t0\t0\tNA\trandom")
