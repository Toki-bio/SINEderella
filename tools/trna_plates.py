"""trna_plates.py (run in <step2_output>/trna by tools/trna_heads.sh, after trna_heads.tsv): one plate per consensus of
the bank - the consensus against its best-scoring tRNA gene (the hit behind best_tRNA in trna_heads.tsv), aligned by
the same ssearch36 that scored it, so the plate shows exactly the alignment behind the numbers.

Rows: the consensus, then the tRNA gene; the aligned stretch in upper case, the unaligned rest in lower case.
Composite stretches (query names with ':') get no plate of their own: their consensus' plate shows its best head.
Out: plates/<consensus>.trna.aln.fa and plates/index.tsv (consensus, tRNA gene, ident%, aln_bp, q_span, bits,
null_max_bits, call: head = bits above every shuffle, below_null otherwise)."""
import collections
import os
import re
import subprocess

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
best = {}
for l in open("ss.m8"):
    f = l.split("\t")
    if "tRNA" in f[1] and (f[0] not in best or float(f[11]) > best[f[0]][0]):
        best[f[0]] = (float(f[11]), f[1], float(f[2]), int(f[3]), "%s-%s" % (f[6], f[7]))
null = collections.defaultdict(float)
for l in open("shuf.m8"):
    f = l.split("\t")
    q = re.sub(r"-shuffled-\d+$", "", f[0])
    null[q] = max(null[q], float(f[11]))


def markx10(qname, tname):
    """ssearch36 -m 10 of one query against one gene -> [(al_start, al_stop, display_start, display), ...] for both"""
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


def row(s, a, b, ds, disp):
    """prefix / aligned display / suffix of a sequence, aligned part upper case"""
    k = ds - 1
    mid = []
    for c in disp:
        if c == "-":
            mid.append("-")
        else:
            k += 1
            mid.append(c.upper() if a <= k <= b else c.lower())
    return s[:ds - 1].lower(), "".join(mid), s[k:].lower()


os.makedirs("plates", exist_ok=True)
with open("plates/index.tsv", "w") as ix:
    ix.write("consensus\ttRNA_gene\tident%\taln_bp\tq_span\tbits\tnull_max_bits\tcall\tplate\n")
    for q in order:
        if ":" in q or q not in best:
            continue
        bits, t, ident, alen, span = best[q]
        (qa, qb, qd, qdisp), (ta, tb, td, tdisp) = markx10(q, t)
        qp, qm, qs = row(seq[q], qa, qb, qd, qdisp)
        tp, tm, ts = row(ref[t], ta, tb, td, tdisp)
        w = max(len(qm), len(tm)); qm, tm = qm.ljust(w, "-"), tm.ljust(w, "-")
        w = max(len(qp), len(tp)); qp, tp = qp.rjust(w, "-"), tp.rjust(w, "-")
        w = max(len(qs), len(ts)); qs, ts = qs.ljust(w, "-"), ts.ljust(w, "-")
        assert (qp + qm + qs).replace("-", "").upper() == seq[q].upper(), q
        assert (tp + tm + ts).replace("-", "").upper() == ref[t], t
        assert len(qp + qm + qs) == len(tp + tm + ts), q
        call = "head" if bits > null.get(q, 0) else "below_null"
        fn = "plates/%s.trna.aln.fa" % re.sub(r"[^A-Za-z0-9_.-]", "_", q)
        with open(fn, "w") as fo:
            fo.write(">%s\n%s\n" % (q, qp + qm + qs))
            fo.write(">%s %.1f%%_over_%dnt_bits_%.1f_shuffle_max_%.1f\n%s\n" % (t, ident, alen, bits, null.get(q, 0), tp + tm + ts))
        ix.write("%s\t%s\t%.1f\t%d\t%s\t%.1f\t%.1f\t%s\t%s\n" % (q, t, ident, alen, span, bits, null.get(q, 0), call,
                                                               os.path.basename(fn)))
for f in ("plate_q.fa", "plate_t.fa"):
    os.remove(f)
