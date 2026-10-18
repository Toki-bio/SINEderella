#!/usr/bin/env python3
"""Per-subfamily border loop: 5 prime then 3 prime (HANDOFF section 68).

Generalizes flank_border_iterate.py (SINE5 prototype). Outward moves only.
Validation slop 1000 bp/side; peel-rebuilt or loop conse consensus preferred.
"""
import argparse
import json
import os
import re
import shlex
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


ISLAND_FRAC_STOP = 0.10   # below this, treat the flank as independent (matches
                          # verdict.py's ISLAND_NOTE threshold -- same bar for
                          # "worth flagging" and "counts as resolved")
MAX_ROUNDS = 5            # rounds per side; each round is one extend+rescan


def _flank_seqs(bed, genome, side, n, work, tag):
    """the n bp next to each element on one side, in the element's orientation (bedtools flank -s; no alignment)"""
    sizes = os.path.join(work, "genome.sizes")
    F.sh("cut -f1,2 %s > %s" % (shlex.quote(genome + ".fai"), shlex.quote(sizes)))
    fb = os.path.join(work, "flank_%s.bed" % tag)
    ff = os.path.join(work, "flank_%s.fa" % tag)
    lr = ("-l %d -r 0" % n) if side == "5prime" else ("-l 0 -r %d" % n)
    F.sh("bedtools flank -s %s -g %s -i %s > %s 2>/dev/null" % (lr, shlex.quote(sizes), shlex.quote(bed), shlex.quote(fb)))
    F.sh("bedtools getfasta -s -nameOnly -fi %s -bed %s > %s 2>/dev/null" % (shlex.quote(genome), shlex.quote(fb), shlex.quote(ff)))
    out = {}
    if os.path.isfile(ff) and os.path.getsize(ff):
        for nm, sq in zip(*F.read_fa(ff)):
            out[re.sub(r"\([+-]\)$", "", nm.split()[0])] = sq.upper()
    return out


def run_side(side, bed, genome, anchor_name, anchor_seq, work, out_log,
             max_rounds=MAX_ROUNDS, island_frac_stop=ISLAND_FRAC_STOP):
    """Iteratively extend one flank until it's independent, or extension runs
    out, or MAX_FLANK/max_rounds is hit -- and says honestly which happened.

    The core (the copies' elements with the anchor consensus) is aligned ONCE per side; each round only takes the flank
    from the genome (bedtools, seconds) and joins it to the core, every copy's flank packed against its own junction
    (column d = d bp from the junction). Owner's design (2026-10-18): the core stays aligned, the expanded flanks are added
    separately and junctioned. Until then every round re-extracted element + the whole 1 kb flank of all copies and ran
    MAFFT on it - and postprocess_flanks then discarded the flank alignment by packing the flank bases against the
    junction anyway, so the tests below already read junction-packed flanks; rsi MEG-RS took 15 min.

    Each round: element_window on core + packed flank (copy-supported extension, measured from the anchor's edge, so a
    round's extension is the growth over the rounds before) -> island scan on the MAX_FLANK bp beyond the current
    element edge (as before). Stops as soon as ANY of:
      - island_fraction < island_frac_stop           -> genuinely independent
      - no further copy-supported extension is found -> extension exhausted
      - total extension for this side hits MAX_FLANK  -> capped, "STILL BAD"
      - max_rounds reached                            -> capped, "STILL BAD"
      - nothing left to scan                          -> untested, "STILL BAD"
    The capped and untested cases are reported as such, not silently accepted as clean.
    """
    bed0 = bed
    ef = os.path.join(work, "core_%s.elements.fa" % side)
    F.extract_elements(bed0, genome, ef, work)
    F.write_fa(os.path.join(work, "anchor.fa"), [(anchor_name, anchor_seq)])
    combined = os.path.join(work, "core_%s.combined.fa" % side)
    F.sh("cat %s %s > %s" % (os.path.join(work, "anchor.fa"), ef, combined))
    core_aln = os.path.join(work, "core_%s.aln.fa" % side)
    F.mafft(combined, core_aln, opts=F.MAFFT_VAL)
    cn, cs = F.read_fa(core_aln)
    ka = next((i for i, n in enumerate(cn) if n.split()[0] == anchor_name), F.consensus_index(cn))
    anchor_row = cs[ka].upper()
    core = {re.sub(r"\([+-]\)$", "", n.split()[0]): s.upper() for i, (n, s) in enumerate(zip(cn, cs)) if i != ka}
    W = len(anchor_row)

    total_extend = 0
    island = None
    stop_reason = None
    verdict = "untested"
    for round_i in range(1, max_rounds + 1):
        need = total_extend + F.MAX_FLANK
        fl = _flank_seqs(bed0, genome, side, need, work, "%s_r%d" % (side, round_i))
        rows, flank_rows = [], []
        for name, crow in core.items():
            f = fl.get(name, "")
            if side == "5prime":
                f = f[-need:] if f else ""
                packed = "-" * (need - len(f)) + f                 # junction-adjacent base in the last column
                rows.append(packed + crow)
                flank_rows.append(packed[::-1])                    # flank_rows: column d = d bp from the junction
            else:
                f = f[:need]
                packed = f + "-" * (need - len(f))
                rows.append(crow + packed)
                flank_rows.append(packed)
        cons = ("-" * need + anchor_row) if side == "5prime" else (anchor_row + "-" * need)
        lo2, hi2, bdiag = F.element_window(cons, rows)
        ext_total = bdiag.get("extended_left" if side == "5prime" else "extended_right", 0)
        ext_bp = max(0, ext_total - total_extend)
        # the flank scanned = MAX_FLANK bp beyond the current element edge (total_extend), shared part included -
        # the same region the old per-round MAFFT version scanned, so the verdicts do not change with the speed-up
        a = total_extend
        block_cols = [list(r[a:a + F.MAX_FLANK]) for r in flank_rows]
        block = F.np.array(block_cols) if block_cols and block_cols[0] else None
        island = F.island_scan_block(block)

        out_log["steps"].append({
            "side": side, "round": round_i, "extend_cols": ext_bp, "extend_bp": ext_bp,
            "extend_total": ext_total, "copy_supported": bdiag, "island": island,
            "extract": "core aligned once (%d cols) + %d bp %s flank packed at the junction, no flank alignment"
                       % (W, need, "left" if side == "5prime" else "right"),
        })

        island_frac = (island or {}).get("island_fraction")
        if island_frac is None:
            # Nothing to scan: no flank columns left beyond the copy-supported edge (the element reaches the end of the
            # window, i.e. the flank is shared all the way out) or under 500 flank bases. A test that could not be done
            # is not a pass (docs/AIM.md). Until 2026-10-17 this case was reported "independent (island_fraction=None)".
            stop_reason = "untested: no flank left to scan -- STILL BAD"
            verdict = "untested"
            break
        if island_frac < island_frac_stop:
            stop_reason = "independent (island_fraction=%s)" % island_frac
            verdict = "independent"
            break
        if ext_bp <= 0:
            stop_reason = ("extension exhausted, still shared (island_fraction=%s) "
                            "-- STILL BAD" % island_frac)
            verdict = "still_shared"
            break
        bed = (F.update_bed_5prime(bed, ext_bp, work) if side == "5prime"
               else update_bed_3prime(bed, ext_bp, work))
        total_extend += ext_bp
        if total_extend >= F.MAX_FLANK:
            stop_reason = ("hit MAX_FLANK cap (%d bp), still shared "
                            "(island_fraction=%s) -- STILL BAD"
                            % (F.MAX_FLANK, island_frac))
            verdict = "still_shared"
            break
    else:
        island_frac = (island or {}).get("island_fraction")
        stop_reason = ("hit max_rounds=%d cap, still shared (island_fraction=%s) "
                        "-- STILL BAD" % (max_rounds, island_frac))
        verdict = "still_shared"

    out_log["stop_reason_%s" % side] = stop_reason
    # read by tools/flank_verdict.py: only "independent" on both sides makes a subfamily interspersed
    out_log["verdict_%s" % side] = verdict
    return bed, total_extend


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
    if final_seq and (ext5 or ext3):      # nothing extended: the searched consensus stays (a rebuild alone shifted TOYB 160 -> 159 bp)
        F.write_fa(reb_path, [(sf, final_seq)])
    else:
        final_seq = None
    adj_bed = os.path.join(out_dir, "loci.adjusted.bed")
    F.sh("cp %s %s" % (bed, adj_bed))
    log["rebuilt_consensus"] = reb_path if final_seq else None
    log["adjusted_bed"] = adj_bed
    log["extend_5prime_bp"] = ext5
    log["extend_3prime_bp"] = ext3
    with open(os.path.join(out_dir, "border_loop.json"), "w") as fh:
        json.dump(log, fh, indent=2)
    still_bad = [s for s in ("5prime", "3prime")
                 if "STILL BAD" in (log.get("stop_reason_%s" % s) or "")]
    print("OK %s extend5=%s extend3=%s consensus=%s%s"
          % (sf, ext5, ext3, reb_path if final_seq else "unchanged",
             ("  STILL_BAD_SIDES=%s" % ",".join(still_bad)) if still_bad else ""))
    print("  5prime: %s" % log.get("stop_reason_5prime"))
    print("  3prime: %s" % log.get("stop_reason_3prime"))


if __name__ == "__main__":
    main()
