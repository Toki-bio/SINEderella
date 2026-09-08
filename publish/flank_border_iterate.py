#!/usr/bin/env python3
"""Border-loop iteration 1 (5 prime) for oma_SINE5 top100.

Proper geometry:
  genome re-extract -> MAFFT (localpair) -> postprocess_flanks (unaligned flanks)
  validation search up to MAX_FLANK bp per side; display 50L/70R
  all consensus methods in ONE alignment (no per-method row-0 splice files)
"""
from __future__ import annotations

import html
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import quote

import numpy as np

OCC, CONS = 0.50, 0.45
MISS_EL = 8
GAPS = set("-.")


def column_stats(rows, j):
    b = [r[j] for r in rows if j < len(r) and r[j] not in GAPS]
    if not b:
        return 0.0, 0.0
    occ = len(b) / float(len(rows))
    u = [x.upper() for x in b]
    top = max(u.count(c) for c in set(u))
    return occ, top / float(len(u))


def element_window(cons, rows):
    L = len(cons)
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if not nz or not rows:
        return 0, L - 1, {}
    old_lo, old_hi = nz[0], nz[-1]
    ok = {}

    def good(j):
        if j not in ok:
            o, c = column_stats(rows, j)
            ok[j] = (o >= OCC and c >= CONS)
        return ok[j]

    core = [j for j in nz if good(j)]
    if not core:
        return old_lo, old_hi, {"note": "no supported consensus column"}
    lo, hi = core[0], core[-1]

    miss, j = 0, lo - 1
    while j >= 0 and miss < MISS_EL:
        if good(j):
            lo, miss = j, 0
        else:
            miss += 1
        j -= 1

    miss, j = 0, hi + 1
    while j < L and miss < MISS_EL:
        if good(j):
            hi, miss = j, 0
        else:
            miss += 1
        j += 1

    d = {
        "old_lo": old_lo, "old_hi": old_hi, "new_lo": lo, "new_hi": hi,
        "extended_left": max(0, old_lo - lo), "extended_right": max(0, hi - old_hi),
        "old_span": old_hi - old_lo + 1, "new_span": hi - lo + 1,
    }
    return lo, hi, d

ENV_PATH = "/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin"
MAX_FLANK = 1000
DISP_L, DISP_R = 50, 70
MIN_RUN = 6
Z_CUT = 8.0
MIN_COPIES = 8
MAFFT_OPTS = (
    "--localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase --quiet"
)
# Validation alignments (1000 bp slop × 100 copies): fast tree, no iterative refinement
MAFFT_VAL = "--retree 1 --maxiterate 0 --nuc --reorder --preservecase --quiet"
MAFFT_DISP = "--auto --quiet --preservecase"


def read_fa(path):
    names, seqs, cur, buf = [], [], None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if cur is not None:
                    seqs.append("".join(buf))
                cur = line[1:]
                names.append(cur)
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


def write_fa(path, pairs):
    with open(path, "w") as fh:
        for n, s in pairs:
            fh.write(">%s\n%s\n" % (n, s))


def sh(cmd, cwd=None, timeout=7200):
    full = "export PATH=%s:$PATH; %s" % (ENV_PATH, cmd)
    r = subprocess.run(full, shell=True, cwd=cwd, capture_output=True, text=True,
                       timeout=timeout)
    return r.returncode, r.stdout, r.stderr


def ungap(s):
    return re.sub(r"[^ACGTacgt]", "", s).upper()


def parse_header(hdr):
    """Return dict ctg, start, end, strand, sf, bits or None."""
    parts = hdr.split("|")
    if len(parts) < 2:
        return None
    loc = parts[0]
    strand = "+"
    if re.search(r"\(-\)$", loc):
        strand = "-"
    loc = re.sub(r"\([^)]*\)$", "", loc)
    m = re.match(r"^(.+):(\d+)-(\d+)$", loc)
    if not m:
        return None
    bits = 0
    if len(parts) > 2:
        try:
            bits = int(float(parts[2]))
        except ValueError:
            bits = 0
    return {
        "ctg": m.group(1), "start": int(m.group(2)), "end": int(m.group(3)),
        "strand": strand, "sf": parts[1], "bits": bits,
        "hdr": hdr,
    }


def load_top100(top100_aln, assigned_fa, sf_name, n=100):
    """Copy list from top100 alignment; fall back to assigned bitscore sort."""
    copies = []
    if top100_aln and os.path.isfile(top100_aln):
        names, _ = read_fa(top100_aln)
        copies = [nm for nm in names if "CONSENSUS" not in nm.upper()]
    if len(copies) >= n:
        return copies[:n]
    recs = []
    name, seq, buf = None, None, []
    with open(assigned_fa) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    p = parse_header(name)
                    if p and p["sf"] == sf_name:
                        recs.append((p["bits"], name))
                name = line[1:]
            elif name:
                pass
    if name:
        p = parse_header(name)
        if p and p["sf"] == sf_name:
            recs.append((p["bits"], name))
    recs.sort(key=lambda t: t[0], reverse=True)
    return [h for _, h in recs[:n]]


def headers_to_bed(headers, work):
    bed = os.path.join(work, "loci.bed")
    with open(bed, "w") as fh:
        for h in headers:
            p = parse_header(h)
            if not p:
                continue
            fh.write("%s\t%d\t%d\t%s\t0\t%s\n" % (
                p["ctg"], p["start"], p["end"], h, p["strand"]))
    return bed


def extract_flanked(bed, genome, flank_l, flank_r, out_fa, work):
    """Genomic extract with optional asymmetric flanks (0 = none on that side)."""
    sizes = os.path.join(work, "genome.sizes")
    sh("cut -f1,2 %s.fai > %s" % (genome, sizes))
    slop = os.path.join(work, "slop.bed")
    code, _, err = sh(
        "bedtools slop -s -l %d -r %d -g %s -i %s > %s 2>/dev/null"
        % (flank_l, flank_r, sizes, bed, slop))
    if code != 0 or not os.path.getsize(slop):
        return False, err
    code, _, err = sh(
        "bedtools getfasta -s -nameOnly -fi %s -bed %s > %s 2>/dev/null"
        % (genome, slop, out_fa))
    if code != 0 or not os.path.getsize(out_fa):
        return False, err
    return True, ""


def extract_for_side(bed, genome, side, flank_bp, out_fa, work):
    """Extract element plus ONE flank only — never both flanks in one validation MAFFT.

    5prime: left flank + element (no right flank in the alignment).
    3prime: element + right flank (no left flank).
    """
    if side in ("5prime", "L", "left"):
        return extract_flanked(bed, genome, flank_bp, 0, out_fa, work)
    if side in ("3prime", "R", "right"):
        return extract_flanked(bed, genome, 0, flank_bp, out_fa, work)
    return extract_flanked(bed, genome, flank_bp, flank_bp, out_fa, work)


def mafft(in_fa, out_fa, threads=8, opts=None):
    o = opts or MAFFT_OPTS
    code, _, err = sh(
        "mafft --thread %d %s %s > %s 2>/dev/null"
        % (threads, o, in_fa, out_fa))
    return code == 0 and os.path.getsize(out_fa) > 0, err


def extract_elements(bed, genome, out_fa, work):
    code, _, err = sh(
        "bedtools getfasta -s -nameOnly -fi %s -bed %s > %s 2>/dev/null"
        % (genome, bed, out_fa))
    return code == 0 and os.path.getsize(out_fa) > 0, err


def postprocess_flanks(aln_file, cons_name):
    """Python port of extract_alignments.sh postprocess_flanks."""
    names, seqs = read_fa(aln_file)
    cons_idx = None
    for i, n in enumerate(names):
        h = n.split()[0]
        if h == cons_name or cons_name in h:
            cons_idx = i
            break
    if cons_idx is None:
        for i, n in enumerate(names):
            if "CONSENSUS" in n.upper():
                cons_idx = i
                break
    if cons_idx is None:
        return False

    cs = seqs[cons_idx]
    alen = len(cs)
    lbound = rbound = 0
    for j, c in enumerate(cs):
        if c not in GAPS:
            if not lbound:
                lbound = j + 1
            rbound = j + 1
    if not lbound:
        lbound, rbound = 1, alen

    ord_idx = [cons_idx] + [i for i in range(len(names)) if i != cons_idx]
    out_lines = []
    for i in ord_idx:
        s = seqs[i]
        out_lines.append(">" + names[i])
        if i == cons_idx:
            for j in range(0, alen, 80):
                out_lines.append(s[j:j + 80])
            continue
        lf_bases = "".join(c.lower() for c in s[:lbound - 1] if c not in GAPS)
        lf_len = lbound - 1
        lf_out = "-" * (lf_len - len(lf_bases)) + lf_bases
        body = s[lbound - 1:rbound]
        rf_bases = "".join(c.lower() for c in s[rbound:] if c not in GAPS)
        rf_len = alen - rbound
        rf_out = rf_bases + "-" * (rf_len - len(rf_bases))
        full = lf_out + body + rf_out
        for j in range(0, len(full), 80):
            out_lines.append(full[j:j + 80])
    with open(aln_file, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")
    return True


def consensus_index(names):
    for i, n in enumerate(names):
        if "CONSENSUS" in n.upper():
            return i
    return 0


def flank_blocks(path, side=None):
    names, seqs = read_fa(path)
    k = consensus_index(names)
    cons = seqs[k]
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if len(nz) < 20:
        return None
    lo, hi = nz[0], nz[-1]
    rows = [s.upper() for i, s in enumerate(seqs) if i != k]
    L = np.array([list(r[:lo]) for r in rows]) if lo > 0 else None
    R = np.array([list(r[hi + 1:]) for r in rows]) if hi + 1 < len(rows[0]) else None
    if side == "L":
        return L
    if side == "R":
        return R
    return L, R


def island_scan_block(block):
    if block is None or not block.size:
        return None
    counts = {b: int((block == b).sum()) for b in "ACGT"}
    tot = sum(counts.values())
    if tot < 500:
        return None
    p0 = sum((c / float(tot)) ** 2 for c in counts.values())
    zs, n_col, n_isl, run = [], 0, 0, 0
    for j in range(block.shape[1]):
        col = block[:, j]
        col = col[(col != "-") & np.isin(col, list("ACGT"))]
        n = len(col)
        if n < MIN_COPIES:
            zs.append(0.0)
            continue
        pairs = n * (n - 1) / 2.0
        _, cnt = np.unique(col, return_counts=True)
        match = float((cnt * (cnt - 1) / 2.0).sum())
        p = match / pairs
        sd = np.sqrt(p0 * (1 - p0) / pairs)
        zs.append((p - p0) / sd if sd > 0 else 0.0)
    hot = np.asarray(zs) > Z_CUT
    for h in hot:
        if h:
            run += 1
        else:
            if run >= MIN_RUN:
                n_isl += 1
                n_col += run
            run = 0
    if run >= MIN_RUN:
        n_isl += 1
        n_col += run
    n_flank = block.shape[1]
    return {
        "islands": n_isl, "island_cols": n_col, "flank_cols": n_flank,
        "island_fraction": round(n_col / float(n_flank), 4) if n_flank else None,
        "max_z": round(float(max(zs)) if zs else 0.0, 1),
    }


def island_scan_side(path, side):
    b = flank_blocks(path, side)
    return island_scan_block(b)


def extend_cols_left(cons, rows, anchor_lo):
    """How many columns left of anchor element look copy-supported (outward move)."""
    extend = 0
    j = anchor_lo - 1
    miss = 0
    while j >= 0 and miss < 8:
        o, c = column_stats(rows, j)
        if o >= OCC and c >= CONS:
            extend += 1
            miss = 0
        else:
            miss += 1
        j -= 1
    return extend


def update_bed_5prime(bed_in, extend_bp, work):
    bed_out = os.path.join(work, "loci_ext.bed")
    sizes = {}
    with open(os.path.join(work, "genome.sizes")) as fh:
        for line in fh:
            ctg, ln = line.rstrip().split("\t")
            sizes[ctg] = int(ln)
    with open(bed_in) as fin, open(bed_out, "w") as fout:
        for line in fin:
            ctg, s, e, name, sc, strand = line.rstrip().split("\t")
            s, e = int(s), int(e)
            if extend_bp <= 0:
                fout.write(line)
                continue
            if strand == "+":
                s = max(0, s - extend_bp)
            else:
                e = min(sizes.get(ctg, e + extend_bp), e + extend_bp)
            fout.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (ctg, s, e, name, sc, strand))
    return bed_out


def run_conse(aln_fa, pct, work):
    os.makedirs(work, exist_ok=True)
    out_cons = os.path.join(work, os.path.basename(aln_fa) + ".cons")
    aln_abs = os.path.abspath(aln_fa)
    code, _, err = sh("conse %s %d" % (aln_abs, pct), cwd=work)
    if not os.path.isfile(out_cons):
        out_cons = aln_abs + ".cons"
    if code != 0 or not os.path.isfile(out_cons):
        return None, err[-300:]
    _, seqs = read_fa(out_cons)
    return ungap(seqs[0]) if seqs else None, None


def run_sine(copies_fa, work, script):
    scw = os.path.join(work, "sine")
    os.makedirs(scw, exist_ok=True)
    sh("cp %s %s/copies.fa" % (copies_fa, scw))
    code, _, err = sh("bash %s copies.fa 100 50 0.01" % (script, ), cwd=scw, timeout=7200)
    hits = list(Path(scw).glob("*_consensus.fasta")) + list(Path(scw).glob("final_consensus.fasta"))
    if not hits:
        return None, err[-300:]
    _, seqs = read_fa(str(hits[0]))
    return ungap(seqs[0]) if seqs else None, None


def element_copies_from_aln(path, lo, hi):
    names, seqs = read_fa(path)
    k = consensus_index(names)
    out = []
    for i, (n, s) in enumerate(zip(names, seqs)):
        if i == k:
            continue
        out.append((re.sub(r"\s+", "_", n)[:90], ungap(s[lo:hi + 1])))
    return out


def verify_unaligned_flanks(path):
    """Return checks on copy rows (consensus row is element-only, no flank)."""
    names, seqs = read_fa(path)
    k = consensus_index(names)
    cons = seqs[k]
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    lo, hi = nz[0], nz[-1]
    ulens = []
    copy_lower = 0
    for i, s in enumerate(seqs):
        if i == k:
            continue
        ul = "".join(c for c in s[:lo] if c not in GAPS)
        ulens.append(len(ul))
        copy_lower += sum(1 for c in s if c.islower())
    spread = max(ulens) - min(ulens) if ulens else 0
    # column test: in MAFFT-aligned flanks, many columns have identical non-gap across copies
    aligned_cols = 0
    checked = 0
    for j in range(max(0, lo - 20), lo):
        col = [s[j] for s in seqs[1:21] if j < len(s) and s[j] not in GAPS]
        if len(col) >= 8:
            checked += 1
            if len(set(col)) == 1:
                aligned_cols += 1
    return {
        "copy_lowercase_bases": copy_lower,
        "left_ungap_len_min": min(ulens) if ulens else 0,
        "left_ungap_len_max": max(ulens) if ulens else 0,
        "left_ungap_spread": spread,
        "identical_col_frac_last20": round(aligned_cols / float(checked or 1), 3),
        "elem_cols": hi - lo + 1,
        "n_consensus_rows": sum(1 for n in names if "CONSENSUS" in n.upper()),
        "n_copies": len(ulens),
    }


def main():
    if len(sys.argv) < 4:
        sys.exit(
            "usage: flank_border_iterate.py TOP100_ALN GENOME OUT_DIR "
            "[PEEL_ALN] [ASSIGNED] [SINE_SCRIPT]")
    top100_aln = sys.argv[1]
    genome = sys.argv[2]
    out_dir = sys.argv[3]
    peel_aln = sys.argv[4] if len(sys.argv) > 4 else ""
    assigned = sys.argv[5] if len(sys.argv) > 5 else ""
    sine_script = sys.argv[6] if len(sys.argv) > 6 else \
        "/staging/tmp/SINEderella/sine_consensus.sh"
    tag = Path(top100_aln).stem.replace(".aln", "")

    os.makedirs(out_dir, exist_ok=True)
    work = os.path.join(out_dir, "work")
    os.makedirs(work, exist_ok=True)
    aln_pub = os.path.join(out_dir, "alignments")
    os.makedirs(aln_pub, exist_ok=True)

    sf = "SINE5"
    if "SINE5" in tag:
        sf = "SINE5"
    headers = load_top100(top100_aln, assigned, sf, 100)
    bed = headers_to_bed(headers, work)

    # --- peel + anno seed consensuses ---
    consensuses = {}
    peel_cons = peel_aln + ".cons" if peel_aln else ""
    peel_work = os.path.join(work, "peel")
    seq = None
    if peel_cons and os.path.isfile(peel_cons) and os.path.getsize(peel_cons) > 0:
        _, ps = read_fa(peel_cons)
        if ps:
            seq = ungap(ps[0])
    if not seq and peel_aln and os.path.isfile(peel_aln):
        seq, _ = run_conse(peel_aln, 50, peel_work)
    if seq:
        consensuses["peel_loci_conse50"] = seq
    _, anno_seqs = read_fa(top100_aln)
    consensuses["anno_seed"] = ungap(anno_seqs[0])

    log = {"side": "5prime", "iteration": 1, "steps": []}

    # --- validation extract: element + ONE flank (never both in one MAFFT) ---
    vf = os.path.join(work, "val_pre.flank.fa")
    extract_for_side(bed, genome, "5prime", MAX_FLANK, vf, work)
    # start from peel or anno for anchor
    anchor_key = "peel_loci_conse50" if "peel_loci_conse50" in consensuses else "anno_seed"
    anchor_name = "CONSENSUS_" + anchor_key
    anchor_seq = consensuses[anchor_key]
    write_fa(os.path.join(work, "anchor.fa"), [(anchor_name, anchor_seq)])
    sh("cat %s/work/anchor.fa %s > %s/work/val_pre.combined.fa"
       % (out_dir, vf, out_dir))
    val_pre_aln = os.path.join(work, "val_pre.aln.fa")
    mafft(os.path.join(work, "val_pre.combined.fa"), val_pre_aln)
    postprocess_flanks(val_pre_aln, anchor_name)
    names, seqs = read_fa(val_pre_aln)
    k = consensus_index(names)
    rows = [s.upper() for i, s in enumerate(seqs) if i != k]
    cons = seqs[k].upper()
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    anchor_lo, anchor_hi = nz[0], nz[-1]
    lo, hi, bdiag = element_window(cons, rows)
    isl_L_pre = island_scan_side(val_pre_aln, "L")
    isl_R_pre = island_scan_side(val_pre_aln, "R")
    extend_bp = bdiag.get("extended_left", 0)
    log["steps"].append({
        "phase": "pre_5prime", "flank_slop": MAX_FLANK,
        "anchor": anchor_key, "anchor_elem_cols": anchor_hi - anchor_lo + 1,
        "copy_supported": bdiag, "island_left": isl_L_pre, "island_right": isl_R_pre,
        "extend_cols_proposed": extend_bp,
    })

    # --- outward border move (5 prime): extend BED starts ---
    bed_ext = update_bed_5prime(bed, extend_bp, work) if extend_bp > 0 else bed

    # --- rebuild consensus from extended element window ---
    vf2 = os.path.join(work, "val_post.flank.fa")
    extract_for_side(bed_ext, genome, "5prime", MAX_FLANK, vf2, work)
    sh("cat %s/work/anchor.fa %s > %s/work/val_post.combined.fa"
       % (out_dir, vf2, out_dir))
    val_post_aln = os.path.join(work, "val_post.aln.fa")
    mafft(os.path.join(work, "val_post.combined.fa"), val_post_aln)
    postprocess_flanks(val_post_aln, anchor_name)
    names2, seqs2 = read_fa(val_post_aln)
    k2 = consensus_index(names2)
    rows2 = [s.upper() for i, s in enumerate(seqs2) if i != k2]
    cons2 = seqs2[k2].upper()
    lo2, hi2, bdiag2 = element_window(cons2, rows2)
    elem_fa = os.path.join(work, "copies_element.fa")
    extract_elements(bed_ext, genome, elem_fa, work)
    elem_aln = os.path.join(work, "copies_element.aln.fa")
    mafft(elem_fa, elem_aln, opts="--retree 2 --maxiterate 0 --quiet")
    for pct in (35, 50):
        seq, err = run_conse(elem_aln, pct, work)
        if seq:
            consensuses["conse_%d" % pct] = seq
    seq, _ = run_sine(elem_fa, work, sine_script)
    if seq:
        consensuses["sine_consensus"] = seq
    seq50, _ = run_conse(elem_aln, 50, work)
    if seq50:
        consensuses["loop_conse50"] = seq50

    isl_L_post = island_scan_side(val_post_aln, "L")
    isl_R_post = island_scan_side(val_post_aln, "R")
    log["steps"].append({
        "phase": "post_5prime", "extend_bp": extend_bp,
        "copy_supported": bdiag2, "island_left": isl_L_post, "island_right": isl_R_post,
        "rebuilt_consensuses": {k: len(v) for k, v in consensuses.items()},
    })

    # --- display alignment: 50L/70R, ALL methods in one file ---
    disp_fa = os.path.join(work, "disp.flank.fa")
    extract_flanked(bed_ext, genome, DISP_L, DISP_R, disp_fa, work)
    cons_order = [
        "loop_conse50", "conse_35", "conse_50", "sine_consensus",
        "peel_loci_conse50", "anno_seed",
    ]
    cons_pairs = []
    for key in cons_order:
        if key in consensuses:
            cons_pairs.append(("CONSENSUS_" + key, consensuses[key]))
    combined = os.path.join(work, "disp.combined.fa")
    cnames, cseqs = read_fa(disp_fa)
    write_fa(combined, cons_pairs + list(zip(cnames, cseqs)))
    disp_aln = os.path.join(aln_pub, tag + "_border_iter1_all_methods.aln.fa")
    ok, err = mafft(combined, disp_aln, opts=MAFFT_DISP)
    if not ok:
        sys.exit("display mafft failed: %s" % err[-500:])
    anchor_disp = cons_pairs[0][0]
    if not postprocess_flanks(disp_aln, anchor_disp):
        postprocess_flanks(disp_aln, cons_pairs[-1][0])
    verify = verify_unaligned_flanks(disp_aln)
    log["display_alignment"] = os.path.basename(disp_aln)
    log["verify"] = verify

    json_path = os.path.join(out_dir, tag + "_border_iter1.json")
    with open(json_path, "w") as fh:
        json.dump(log, fh, indent=2)

    # --- HTML ---
    gh_base = sys.argv[7] if len(sys.argv) > 7 else \
        "https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/oma_border_iter1"
    msa_url = ("https://toki-bio.github.io/MSA-viewer/?url="
               + quote(gh_base.rstrip("/") + "/" + log["display_alignment"], safe="")
               + "&title=" + quote("oma SINE5 border iter1 all methods"))

    pre = log["steps"][0]
    post = log["steps"][1]

    def isl_cells(d):
        if not d:
            return "<td>—</td><td>—</td><td>—</td><td>—</td>"
        return ("<td class='n'>%d</td><td class='n'>%d</td><td class='n'>%d</td>"
                "<td class='n'>%s</td>" % (
                    d["islands"], d["island_cols"], d["flank_cols"], d["island_fraction"]))

    cons_list = ", ".join(html.escape(k) for k in cons_order if k in consensuses)
    page = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>%s — border loop iter 1 (5 prime)</title>
<style>
body{font-family:system-ui,sans-serif;max-width:960px;margin:24px auto;padding:0 16px}
table{border-collapse:collapse;width:100%%;font-size:13px;margin:12px 0}
th,td{border:1px solid #ccc;padding:6px 8px} td.n{text-align:right;font-family:monospace}
.note{background:#f5f5f0;border-left:3px solid #1f6f5c;padding:10px 14px}
code{font-family:monospace;background:#eee;padding:1px 4px}
</style></head><body>
<h1>%s — border loop iteration 1 (5 prime)</h1>
<p class="note">Genome re-extract → MAFFT localpair → postprocess_flanks. Validation
slop %d bp/side; display %dL/%dR. Outward 5 prime extension: <strong>%d bp</strong>
(from %d copy-supported columns). All consensus methods in one alignment.</p>
<h2>Island scan (per side, not pooled)</h2>
<table><thead><tr><th>phase</th><th>side</th><th>islands</th><th>island cols</th>
<th>flank cols</th><th>fraction</th></tr></thead><tbody>
<tr><td>pre</td><td>L</td>%s</tr>
<tr><td>pre</td><td>R</td>%s</tr>
<tr><td>post</td><td>L</td>%s</tr>
<tr><td>post</td><td>R</td>%s</tr>
</tbody></table>
<h2>Element span (copy-supported boundary walk)</h2>
<p>Pre: %d bp → post: %d bp</p>
<h2>Verification</h2>
<ul>
<li>copy lowercase bases (flanks postprocessed): %d</li>
<li>copy left-flank ungapped lengths: %d–%d bp</li>
<li>identical cols in last 20 flank columns: %s (want low)</li>
<li>consensus rows: %d; copies: %d</li>
</ul>
<h2>One alignment, all methods</h2>
<p><a href="%s">MSA-viewer</a></p>
<p>Consensus rows: %s</p>
</body></html>""" % (
        tag, tag, MAX_FLANK, DISP_L, DISP_R, extend_bp, extend_bp,
        isl_cells(pre.get("island_left")), isl_cells(pre.get("island_right")),
        isl_cells(post.get("island_left")), isl_cells(post.get("island_right")),
        bdiag.get("new_span", 0), bdiag2.get("new_span", 0),
        verify.get("copy_lowercase_bases", 0),
        verify.get("left_ungap_len_min"), verify.get("left_ungap_len_max"),
        verify.get("identical_col_frac_last20"),
        verify.get("n_consensus_rows", 0), verify.get("n_copies", 0),
        msa_url, cons_list,
    )
    html_path = os.path.join(out_dir, tag + "_border_iter1.html")
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(page)
    print("Wrote", html_path)
    print("Verify:", json.dumps(verify))
    print("Extend bp:", extend_bp)


if __name__ == "__main__":
    main()
