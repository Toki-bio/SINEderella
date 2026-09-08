#!/usr/bin/env python3
"""Per-subfamily border loop: 5 prime then 3 prime (HANDOFF section 68).

Generalizes flank_border_iterate.py (SINE5 prototype). Outward moves only.
Validation slop 1000 bp/side; peel-rebuilt or loop conse consensus preferred.
"""
import argparse
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

import flank_border_iterate as F  # noqa: E402


def extend_3prime_cols(cons, rows, lo, hi):
    """Columns supported 3 prime of consensus core within alignment."""
    nz = [i for i, c in enumerate(cons) if c not in F.GAPS]
    if not nz:
        return 0
    ok = {}

    def good(j):
        if j not in ok:
            o, c = F.column_stats(rows, j)
            ok[j] = (o >= F.OCC and c >= F.CONS)
        return ok[j]

    extend = 0
    miss = 0
    j = hi + 1
    while j < len(cons) and miss < F.MISS_EL:
        if good(j):
            extend = j - hi
            miss = 0
        else:
            miss += 1
        j += 1
    return extend


def extend_5prime_cols(cons, rows, lo):
    nz = [i for i, c in enumerate(cons) if c not in F.GAPS]
    if not nz:
        return 0
    ok = {}

    def good(j):
        if j not in ok:
            o, c = F.column_stats(rows, j)
            ok[j] = (o >= F.OCC and c >= F.CONS)
        return ok[j]

    extend = 0
    miss = 0
    j = lo - 1
    while j >= 0 and miss < F.MISS_EL:
        if good(j):
            extend = lo - j
            miss = 0
        else:
            miss += 1
        j -= 1
    return extend


def update_bed_3prime(bed_in, extend_bp, work):
    bed_out = os.path.join(work, "loci_ext3.bed")
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
                e = min(sizes.get(ctg, e + extend_bp), e + extend_bp)
            else:
                s = max(0, s - extend_bp)
            fout.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (ctg, s, e, name, sc, strand))
    return bed_out


def rebuild_consensus(bed, genome, work, sine_script):
    elem_fa = os.path.join(work, "copies_element.fa")
    if not F.extract_elements(bed, genome, elem_fa, work)[0]:
        return None, "extract_elements failed"
    elem_aln = os.path.join(work, "copies_element.aln.fa")
    F.mafft(elem_fa, elem_aln, opts="--retree 2 --maxiterate 0 --quiet")
    seq, _ = F.run_conse(elem_aln, 50, work)
    if not seq:
        seq, _ = F.run_sine(elem_fa, work, sine_script)
    return seq, None


def run_side(side, bed, genome, anchor_name, anchor_seq, work, out_log):
    """Validate one flank: extract element + that flank only, then score."""
    vf = os.path.join(work, "val_%s.flank.fa" % side)
    F.extract_for_side(bed, genome, side, F.MAX_FLANK, vf, work)
    F.write_fa(os.path.join(work, "anchor.fa"), [(anchor_name, anchor_seq)])
    combined = os.path.join(work, "val_%s.combined.fa" % side)
    F.sh("cat %s %s > %s" % (os.path.join(work, "anchor.fa"), vf, combined))
    aln = os.path.join(work, "val_%s.aln.fa" % side)
    F.mafft(combined, aln, opts=F.MAFFT_VAL)
    F.postprocess_flanks(aln, anchor_name)
    names, seqs = F.read_fa(aln)
    k = F.consensus_index(names)
    rows = [s.upper() for i, s in enumerate(seqs) if i != k]
    cons = seqs[k].upper()
    nz = [i for i, c in enumerate(cons) if c not in F.GAPS]
    lo, hi = nz[0], nz[-1]
    lo2, hi2, bdiag = F.element_window(cons, rows)
    scan_side = "L" if side == "5prime" else "R"
    if side == "5prime":
        ext_bp = bdiag.get("extended_left", 0)
        bed_new = F.update_bed_5prime(bed, ext_bp, work) if ext_bp else bed
    else:
        ext_bp = bdiag.get("extended_right", 0)
        bed_new = update_bed_3prime(bed, ext_bp, work) if ext_bp else bed
    out_log["steps"].append({
        "side": side, "extend_cols": ext_bp, "extend_bp": ext_bp,
        "copy_supported": bdiag,
        "island": F.island_scan_side(aln, scan_side),
        "extract": "element+%s_flank_only" % ("left" if side == "5prime" else "right"),
    })
    return bed_new, ext_bp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_root")
    ap.add_argument("subfam")
    ap.add_argument("--peel-aln", default="", help="peel loci alignment for conse")
    ap.add_argument("--sine-script",
                    default="/staging/tmp/SINEderella/sine_consensus.sh")
    args = ap.parse_args()
    run_root = args.run_root.rstrip("/")
    sf = args.subfam
    genome = os.path.join(run_root, "genome.clean.fa")
    assigned = os.path.join(run_root, "step2/step2_output/assigned.fasta")
    out_dir = os.path.join(run_root, "border_loop", sf)
    work = os.path.join(out_dir, "work")
    os.makedirs(work, exist_ok=True)

    headers = F.load_top100("", assigned, sf, 100)
    if not headers:
        sys.exit("no copies for %s" % sf)
    bed = F.headers_to_bed(headers, work)

    anchor_seq = None
    peel_cons = args.peel_aln + ".cons" if args.peel_aln else ""
    if peel_cons and os.path.isfile(peel_cons) and os.path.getsize(peel_cons) > 0:
        _, ps = F.read_fa(peel_cons)
        if ps:
            anchor_seq = F.ungap(ps[0])
    if not anchor_seq and args.peel_aln and os.path.isfile(args.peel_aln):
        anchor_seq, _ = F.run_conse(args.peel_aln, 50, work)
    if not anchor_seq:
        cons_fa = os.path.join(run_root, "consensuses.clean.fa")
        if os.path.isfile(cons_fa):
            for nm, sq in zip(*F.read_fa(cons_fa)):
                if nm.split()[0] == sf:
                    anchor_seq = F.ungap(sq)
                    break
    if not anchor_seq:
        sys.exit("no anchor consensus for %s" % sf)

    anchor_name = "CONSENSUS_peel_or_seed"
    log = {"subfam": sf, "steps": []}

    bed, ext5 = run_side("5prime", bed, genome, anchor_name, anchor_seq, work, log)
    seq5, err = rebuild_consensus(bed, genome, work, args.sine_script)
    anchor2 = seq5 or anchor_seq
    bed, ext3 = run_side("3prime", bed, genome, anchor_name, anchor2, work, log)
    final_seq, err2 = rebuild_consensus(bed, genome, work, args.sine_script)

    reb_dir = os.path.join(run_root, "rebuilt_consensus")
    os.makedirs(reb_dir, exist_ok=True)
    reb_path = os.path.join(reb_dir, "%s.fa" % sf)
    if final_seq:
        F.write_fa(reb_path, [(sf, final_seq)])
    adj_bed = os.path.join(out_dir, "loci.adjusted.bed")
    F.sh("cp %s %s" % (bed, adj_bed))
    log["rebuilt_consensus"] = reb_path if final_seq else None
    log["adjusted_bed"] = adj_bed
    log["extend_5prime_bp"] = ext5
    log["extend_3prime_bp"] = ext3
    with open(os.path.join(out_dir, "border_loop.json"), "w") as fh:
        json.dump(log, fh, indent=2)
    print("OK %s extend5=%s extend3=%s consensus=%s"
          % (sf, ext5, ext3, reb_path if final_seq else "unchanged"))


if __name__ == "__main__":
    main()
