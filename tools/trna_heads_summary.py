"""trna_heads_summary.py (run in <step2_output>/trna by tools/trna_heads.sh): one row per query -> stdout (trna_heads.tsv).

Columns: query, len, best_tRNA (amino acid-anticodon of the best tRNA gene hit), ident%, aln_bp, q_span, bits,
null_max_bits (best of the dinucleotide shuffles), SINEBase_60/60 (>= 60 % over >= 60 nt), 2nd_isotype (best hit of a
different amino acid, bits), Rfam (best two cmscan hits), tRNAscan (type-anticodon span score ... isotype model), A_box, B_box."""
import collections
import re

order, seq, n = [], {}, None
for l in open("queries.fa"):
    l = l.strip()
    if l.startswith(">"):
        n = l[1:]; order.append(n); seq[n] = ""
    else:
        seq[n] += l


def iso(t):
    m = re.search(r"tRNA-([A-Za-z]+)-([A-Z]{3})", t)
    return "%s-%s" % (m.group(1), m.group(2)) if m else t.split("_")[-1][:12]


hits = collections.defaultdict(list)
for l in open("ss.m8"):
    f = l.split("\t")
    hits[f[0]].append((float(f[11]), float(f[2]), int(f[3]), int(f[6]), int(f[7]), f[1]))
null = collections.defaultdict(float)
for l in open("shuf.m8"):
    f = l.split("\t")
    q = re.sub(r"-shuffled-\d+$", "", f[0])
    null[q] = max(null[q], float(f[11]))
cm = collections.defaultdict(list)
for l in (open("cm.tbl") if __import__("os").path.isfile("cm.tbl") else []):
    if l.startswith("#"):
        continue
    f = l.split()
    cm[f[2]].append((float(f[15]), f[0], int(f[7]), int(f[8]), float(f[14])))
ts = collections.defaultdict(list)
for l in (open("trnascan.out") if __import__("os").path.isfile("trnascan.out") else []):
    f = l.split("\t") if "\t" in l else l.split()
    if len(f) < 9 or not f[1].strip().isdigit():
        continue
    name = f[0].strip()
    a, b = int(f[2]), int(f[3])
    typ, ac, score = f[4].strip(), f[5].strip(), f[8].strip()
    note = " ".join(x.strip() for x in f[9:] if x.strip())
    ts[name].append("%s-%s %s-%s %s%s" % (typ, ac, min(a, b), max(a, b), score, (" " + note) if note else ""))
BOXA = re.compile(r"T[AG]G[CT]..A[AG][CT]GG")      # TRGYnnARYGG (A box, Tatosyan et al. 2020 style)
BOXB = re.compile(r"GGTTC[AG]A.[CT]C[CT]")         # GGTTCRAnYCY (B box)

print("query\tlen\tbest_tRNA\tident%\taln_bp\tq_span\tbits\tnull_max_bits\tSINEBase_60/60\t2nd_isotype\tRfam\ttRNAscan\tA_box\tB_box")
for q in order:
    hs = sorted(hits.get(q, []), reverse=True)
    trna = [h for h in hs if "tRNA" in h[5]]
    if trna:
        b = trna[0]
        best = iso(b[5])
        sec = next((h for h in trna if iso(h[5]).split("-")[0] != best.split("-")[0]), None)
        ok = "yes" if (b[1] >= 60 and b[2] >= 60) else "no"
        row = [best, "%.1f" % b[1], str(b[2]), "%d-%d" % (b[3], b[4]), "%.1f" % b[0], "%.1f" % null.get(q, 0), ok,
               ("%s %.1f" % (iso(sec[5]), sec[0])) if sec else "-"]
    else:
            row = ["-", "-", "-", "-", "-", "%.1f" % null.get(q, 0), "no", "-"]
    c = sorted(cm.get(q, []))
    rf = "; ".join("%s %d-%d E=%.0e" % (x[1], x[2], x[3], x[0]) for x in c[:2]) or "-"
    tr = " | ".join(ts.get(q.split()[0], [])) or "-"
    s = seq[q].upper()
    ab = ",".join(str(m.start() + 1) for m in BOXA.finditer(s)) or "-"
    bb = ",".join(str(m.start() + 1) for m in BOXB.finditer(s)) or "-"
    print("\t".join([q, str(len(s))] + row + [rf, tr, ab, bb]))
