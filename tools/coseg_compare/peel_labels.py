#!/usr/bin/env python3
"""Glue between SubFam chunks and the peel (SINE-discriminator peel_features.py).

usage: peel_labels.py truth CHUNKS.tsv LABELS.tsv OUT_TRUTH.json
         writes the majority label of each chunk (peel_features.py takes it to report purity)
       peel_labels.py groups CHUNKS.tsv PEEL.json OUT.tsv
         writes  copy id <TAB> group  for every copy whose chunk was peeled into a group;
         copies in unpeeled chunks (the residue) are left out, so coverage shows what was not placed.
CHUNKS.tsv = SubFam's  copy id, chunk, strand.
"""
import collections
import json
import sys

mode = sys.argv[1]
if mode == "truth":
    chunks_tsv, labels_tsv, out = sys.argv[2:5]
    lab = dict(l.rstrip("\n").split("\t")[:2] for l in open(labels_tsv))
    by = collections.defaultdict(collections.Counter)
    for l in open(chunks_tsv):
        f = l.rstrip("\n").split("\t")
        if f[0] in lab:
            by[f[1]][lab[f[0]]] += 1
    json.dump({c: cnt.most_common(1)[0][0] for c, cnt in by.items()}, open(out, "w"), indent=1)
elif mode == "groups":
    chunks_tsv, peel_json, out = sys.argv[2:5]
    pj = json.load(open(peel_json))
    chunk_group = {}
    for k, g in enumerate(pj["peeled"], 1):
        for m in g["members"]:
            chunk_group.setdefault(m, "pg%d" % k)
    n_in = n_out = 0
    with open(out, "w") as fo:
        for l in open(chunks_tsv):
            f = l.rstrip("\n").split("\t")
            n_in += 1
            if f[1] in chunk_group:
                fo.write("%s\t%s\n" % (f[0], chunk_group[f[1]])); n_out += 1
    print("%d groups from the peel; %d of %d copies placed" % (len(pj["peeled"]), n_out, n_in), file=sys.stderr)
else:
    sys.exit(__doc__)
