#!/usr/bin/env python3
"""List subfamilies that need the element border loop before publish.

Such cases (per HANDOFF §68 + peel step4 border screen):
  - step7 flank test hit 1000 bp without confirming background (undetermined)
  - peel / flank screen flagged both-sided or one-sided (optional TSV)
  - copy-supported element extends past consensus on a 1000 bp validation
    alignment (--scan; top100 only)
"""
import argparse
import glob
import os
import subprocess
import sys

GAPS = set("-.")


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


def from_boundary_tsv(path):
    need = set()
    if not path or not os.path.isfile(path):
        return need
    with open(path) as fh:
        hdr = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            sf, _pop, _side, _bp, status = parts[:5]
            if status in ("undetermined", "insufficient"):
                need.add(sf)
    return need


def from_peel_flags(path):
    need = set()
    if not path or not os.path.isfile(path):
        return need
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sf = line.split("\t")[0].strip()
            if sf:
                need.add(sf)
    return need


def element_extend_scan(run_root, subfams, genome, assigned, flank=1000):
    """Quick MAFFT-free column walk on 1000 bp extracts for top100 copies."""
    disc = os.environ.get("DISC", "/staging/tmp/sinedisc")
    sys.path.insert(0, disc)
    import boundary as B  # noqa

    need = set()
    work = os.path.join(run_root, ".border_scan_tmp")
    os.makedirs(work, exist_ok=True)
    sizes = os.path.join(work, "genome.sizes")
    if not os.path.isfile(sizes):
        subprocess.check_call(
            "cut -f1,2 %s.fai > %s" % (genome, sizes), shell=True)

    for sf in sorted(subfams):
        headers = []
        with open(assigned) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                h = line[1:].strip()
                if "|%s|" % sf in h or h.split("|")[1:2] == [sf]:
                    headers.append(h)
        headers.sort(key=lambda h: float(h.split("|")[-1].replace("(+)", "")
                                         .replace("(-)", "") or 0), reverse=True)
        headers = headers[:100]
        if len(headers) < 8:
            continue
        bed = os.path.join(work, "%s.bed" % sf)
        with open(bed, "w") as fh:
            for h in headers:
                parts = h.split("|")[0]
                strand = "+"
                if parts.endswith("(-)"):
                    strand = "-"
                loc = parts.replace("(+)", "").replace("(-)", "")
                ctg, coords = loc.rsplit(":", 1)
                start_s, end_s = coords.rsplit("-", 1)
                fh.write("%s\t%d\t%d\t%s\t0\t%s\n"
                         % (ctg, int(start_s) - 1, int(end_s), h, strand))
        slop = os.path.join(work, "%s.slop.bed" % sf)
        fa = os.path.join(work, "%s.fa" % sf)
        subprocess.call(
            "bedtools slop -s -l %d -r %d -g %s -i %s > %s 2>/dev/null"
            % (flank, flank, sizes, bed, slop), shell=True)
        subprocess.call(
            "bedtools getfasta -s -nameOnly -fi %s -bed %s > %s 2>/dev/null"
            % (genome, slop, fa), shell=True)
        if not os.path.getsize(fa):
            continue
        cons_seq = None
        cons_fa = os.path.join(run_root, "consensuses.clean.fa")
        if os.path.isfile(cons_fa):
            cnames, cseqs = read_fa(cons_fa)
            for nm, sq in zip(cnames, cseqs):
                if nm.split()[0] == sf:
                    cons_seq = sq
                    break
        names, seqs = read_fa(fa)
        if not seqs:
            continue
        if cons_seq:
            # pad to same length via external mafft is expensive; skip if no cons
            rows = [s.upper() for s in seqs]
            cons = cons_seq.upper().ljust(max(len(s) for s in rows), "-")
            _, _, d = B.element_window(cons, rows)
            if d.get("extended_left", 0) > 0 or d.get("extended_right", 0) > 0:
                need.add(sf)
    return need


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_root")
    ap.add_argument("--scan", action="store_true",
                    help="scan top100 for copy-supported element extension")
    ap.add_argument("--peel-flags",
                    help="TSV: subfamily per line (both/one-sided flags)")
    args = ap.parse_args()
    run_root = args.run_root.rstrip("/")
    tsv = os.path.join(run_root, "step2/step2_output/boundary_refinement.tsv")
    need = from_boundary_tsv(tsv)
    need |= from_peel_flags(args.peel_flags)

    if args.scan:
        assigned = os.path.join(run_root, "step2/step2_output/assigned.fasta")
        genome = os.path.join(run_root, "genome.clean.fa")
        if os.path.isfile(assigned) and os.path.isfile(genome):
            all_sf = sorted({ln.split("|")[1] for ln in open(assigned)
                             if ln.startswith(">") and "|" in ln})
            try:
                need |= element_extend_scan(run_root, all_sf, genome, assigned)
            except Exception as exc:
                print("WARN border scan failed: %s" % exc, file=sys.stderr)

    for sf in sorted(need):
        print(sf)


if __name__ == "__main__":
    main()
