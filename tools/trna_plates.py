"""trna_plates.py (run in <step2_output>/trna by tools/trna_heads.sh, after trna_heads.tsv): one plate per consensus of
the bank - the consensus against its best-scoring tRNA gene (the hit behind best_tRNA in trna_heads.tsv), aligned by
the same ssearch36 that scored it, so the plate shows exactly the alignment behind the numbers. When the best gene of a
different amino acid scores within TIE_BITS (2) bits of the best, the isotype is not resolved and that gene is a third
row (both alignments to the consensus, merged on the consensus coordinates).

Rows: the consensus, then the tRNA gene(s); aligned stretches in upper case, the unaligned rest in lower case.
Composite stretches (query names with ':') get no plate of their own: their consensus' plate shows its best head.
Out: plates/<consensus>.trna.aln.fa and plates/index.tsv (consensus, tRNA gene, ident%, aln_bp, q_span, bits,
null_max_bits, call: head = bits above every shuffle, below_null otherwise; tie_gene, tie_ident%, tie_aln_bp, tie_bits
or '-')."""
import collections
import os
import re
import subprocess

TIE_BITS = float(os.environ.get("TIE_BITS", "2"))

seq, order, n = {}, [], None
for l in open("queries.fa"):
    l = l.strip()
    if l.startswith(">"):
        n = l[1:]; order.append(n); seq[n] = ""
    else:
        seq[n] += l
ref, rn = {}, None
for l in open("ref.fa"):
    l = l.rstrip("\n")
    if l.startswith(">"):
        rn = l[1:].split()[0]; ref[rn] = ""
    else:
        ref[rn] += l.strip().upper()


def aa(t):
    m = re.search(r"tRNA-([A-Za-z]+)-", t)
    return m.group(1) if m else t


hits = collections.defaultdict(list)
for l in open("ss.m8"):
    f = l.split("\t")
    if "tRNA" in f[1]:
        hits[f[0]].append((float(f[11]), f[1], float(f[2]), int(f[3]), "%s-%s" % (f[6], f[7])))
null = collections.defaultdict(float)
for l in open("shuf.m8"):
    f = l.split("\t")
    q = re.sub(r"-shuffled-\d+$", "", f[0])
    null[q] = max(null[q], float(f[11]))


def markx10(qname, tname):
    """ssearch36 -m 10 of one query against one gene -> [(al_start, al_stop, display_start, display)] for query, gene"""
    with open("plate_q.fa", "w") as fo:
        fo.write(">%s\n%s\n" % (qname, seq[qname]))
    with open("plate_t.fa", "w") as fo:
        fo.write(">%s\n%s\n" % (tname, ref[tname]))
    out = subprocess.run(["ssearch36", "-3", "-m", "10", "-b", "1", "-d", "1", "plate_q.fa", "plate_t.fa"],
                         capture_output=True, text=True, check=True).stdout
    part = out.split("\n>>", 1)[1].split("\n>>>", 1)[0]       # the one hit
    recs = []
    for blk in re.split(r"\n>(?!>)", part)[1:3]:               # query block, library block
        blk = blk.split("\n; al_cons")[0]                     # the match-marker lines follow it
        info = dict(re.findall(r"^; (\w+): (\S+)", blk, re.M))
        disp = "".join(x.strip() for x in blk.split("\n")[1:] if x.strip() and not x.startswith(";"))
        recs.append((int(info["al_start"]), int(info["al_stop"]), int(info["al_display_start"]), disp))
    return recs


def on_consensus(q, t):
    """the gene laid on the consensus: col[i] = gene letter opposite consensus residue i (or '-'), ins[i] = gene letters
    inserted before consensus residue i (ins[L] = after the end); aligned gene letters upper case, the rest lower case.
    Returns col, ins, consensus aligned span (1-based)."""
    (qa, qb, qd, qdisp), (ta, tb, td, tdisp) = markx10(q, t)
    g, L = ref[t], len(seq[q])
    col, ins = ["-"] * L, [""] * (L + 1)
    c = lambda k: g[k].upper() if ta <= k + 1 <= tb else g[k].lower()
    qk, tk = qd - 1, td - 1
    ins[qk] += "".join(c(k) for k in range(tk))              # gene 5' part not displayed
    for j in range(max(len(qdisp), len(tdisp))):
        x = qdisp[j] if j < len(qdisp) else "-"
        y = tdisp[j] if j < len(tdisp) else "-"
        if x != "-" and qk >= L:
            x = "-"
        if x != "-" and y != "-":
            col[qk] = c(tk); qk += 1; tk += 1
        elif x != "-":
            qk += 1
        elif y != "-":
            ins[qk] += c(tk); tk += 1
    ins[qk] += "".join(c(k) for k in range(tk, len(g)))     # gene 3' part not displayed
    return col, ins, (qa, qb)


def plate(q, genes):
    """rows [(name, aligned string)] of the consensus and every gene, merged on the consensus coordinates"""
    lays = [on_consensus(q, t) for t in genes]
    qa, qb = lays[0][2]
    s, L = seq[q], len(seq[q])
    first = min(min(i for i in range(L + 1) if ins[i] or (i < L and col[i] != "-")) for col, ins, _ in lays)
    rows = [[] for _ in range(len(genes) + 1)]
    for i in range(L + 1):
        w = max(len(ins[i]) for _, ins, _ in lays)
        if w:
            rows[0].append("-" * w)
            for r, (_, ins, _) in zip(rows[1:], lays):
                r.append(ins[i].rjust(w, "-") if i <= first else ins[i].ljust(w, "-"))
        if i < L:
            rows[0].append(s[i].upper() if qa <= i + 1 <= qb else s[i].lower())
            for r, (col, _, _) in zip(rows[1:], lays):
                r.append(col[i])
    rows = ["".join(r) for r in rows]
    assert rows[0].replace("-", "").upper() == s.upper(), q
    for t, r in zip(genes, rows[1:]):
        assert r.replace("-", "").upper() == ref[t], t
    assert len(set(map(len, rows))) == 1, q
    return rows


os.makedirs("plates", exist_ok=True)
with open("plates/index.tsv", "w") as ix:
    ix.write("consensus\ttRNA_gene\tident%\taln_bp\tq_span\tbits\tnull_max_bits\tcall\ttie_gene\ttie_ident%\ttie_aln_bp\t"
             "tie_bits\tplate\n")
    for q in order:
        if ":" in q or not hits.get(q):
            continue
        hs = sorted(hits[q], key=lambda h: -h[0])           # equal scores: first hit, as ssearch36 lists them
        best = hs[0]
        tie = next((h for h in hs if aa(h[1]) != aa(best[1])), None)
        if tie and best[0] - tie[0] > TIE_BITS:
            tie = None
        genes = [best[1]] + ([tie[1]] if tie else [])
        rows = plate(q, genes)
        fn = "plates/%s.trna.aln.fa" % re.sub(r"[^A-Za-z0-9_.-]", "_", q)
        with open(fn, "w") as fo:
            fo.write(">%s\n%s\n" % (q, rows[0]))
            for h, r in zip([best] + ([tie] if tie else []), rows[1:]):
                fo.write(">%s %.1f%%_over_%dnt_bits_%.1f_shuffle_max_%.1f\n%s\n" % (h[1], h[2], h[3], h[0], null.get(q, 0), r))
        call = "head" if best[0] > null.get(q, 0) else "below_null"
        tcols = ("%s\t%.1f\t%d\t%.1f" % (tie[1], tie[2], tie[3], tie[0])) if tie else "-\t-\t-\t-"
        ix.write("%s\t%s\t%.1f\t%d\t%s\t%.1f\t%.1f\t%s\t%s\t%s\n" % (q, best[1], best[2], best[3], best[4], best[0],
                                                                 null.get(q, 0), call, tcols, os.path.basename(fn)))
for f in ("plate_q.fa", "plate_t.fa"):
    os.remove(f)
