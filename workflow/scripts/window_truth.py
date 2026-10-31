"""Truth sequences of the windows: the matching stretch of each HG002 Q100 haplotype.

Each window boundary is anchored by --anchor-length bp of reference sequence just outside it (left of the start,
right of the end), aligned to the haplotype (minimap2 asm5 via mappy). An anchor must align uniquely (MAPQ >= 30,
>= 90% of its length); if not, it is moved outward by the next --steps offset and the boundary moves with it.
Both anchors must land on the same contig and strand, in order, with a haplotype span of 0.5-2x the reference
span. The stretch between them is written in reference orientation.

Per window two truths are written: the window itself (scoring) and the window +/- --flank (assigning reads to
haplotypes in the oracle mode). Output: <outdir>/<window>.<core|extended>.hap<1|2>.fa and a status table.
"""

import argparse
import os

import mappy
import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--hap1", required=True)
    parser.add_argument("--hap2", required=True)
    parser.add_argument("--flank", type=int, required=True)
    parser.add_argument("--anchor-length", type=int, required=True)
    parser.add_argument("--steps", type=int, nargs="+", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--status", required=True)
    return parser.parse_args()


def unique_hit(aligner, seq):
    hits = [h for h in aligner.map(seq) if h.is_primary]
    if len(hits) != 1:
        return None
    hit = hits[0]
    if hit.mapq < 30 or hit.q_en - hit.q_st < 0.9 * len(seq):
        return None
    return hit


def anchor(aligner, reference, chrom, boundary, side, length, steps):
    """(reference boundary, hit, haplotype position of the boundary) for the first unique anchor, or None.
    side 'start': the anchor is left of the boundary; 'end': right of it."""
    chrom_len = reference.get_reference_length(chrom)
    for step in steps:
        if side == "start":
            pos = boundary - step
            lo, hi = pos - length, pos
        else:
            pos = boundary + step
            lo, hi = pos, pos + length
        if lo < 0 or hi > chrom_len:
            continue
        seq = reference.fetch(chrom, lo, hi).upper()
        if seq.count("N") > 0.1 * length:
            continue
        hit = unique_hit(aligner, seq)
        if hit is None:
            continue
        # haplotype coordinate of the boundary (the anchor's inner end)
        if side == "start":
            hap_pos = hit.r_en + (length - hit.q_en) if hit.strand == 1 else hit.r_st - (length - hit.q_en)
        else:
            hap_pos = hit.r_st - hit.q_st if hit.strand == 1 else hit.r_en + hit.q_st
        return pos, hit, hap_pos
    return None


def extract(aligner, haplotype, reference, chrom, start, end, args):
    left = anchor(aligner, reference, chrom, start, "start", args.anchor_length, args.steps)
    right = anchor(aligner, reference, chrom, end, "end", args.anchor_length, args.steps)
    if left is None or right is None:
        return None, "anchor not unique" + (" (start)" if left is None else " (end)")
    (ref_start, lhit, lpos), (ref_end, rhit, rpos) = left, right
    if lhit.ctg != rhit.ctg or lhit.strand != rhit.strand:
        return None, "anchors on different contigs or strands"
    lo, hi = (lpos, rpos) if lhit.strand == 1 else (rpos, lpos)
    if hi <= lo:
        return None, "anchors out of order"
    ratio = (hi - lo) / (ref_end - ref_start)
    if not 0.5 <= ratio <= 2:
        return None, f"haplotype span {hi - lo} vs reference {ref_end - ref_start}"
    seq = haplotype.fetch(lhit.ctg, max(0, lo), hi)
    if lhit.strand == -1:
        seq = mappy.revcomp(seq)
    return (seq, ref_start, ref_end, lhit.ctg, lo, hi, lhit.strand), "ok"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    reference = pysam.FastaFile(args.reference)
    windows = []
    with open(args.windows) as fh:
        next(fh)
        for line in fh:
            name, chrom, start, end, _ = line.rstrip("\n").split("\t")
            windows.append((name, chrom, int(start), int(end)))
    rows = []
    for hap, path in (("1", args.hap1), ("2", args.hap2)):
        aligner = mappy.Aligner(path, preset="asm5", n_threads=args.threads)
        haplotype = pysam.FastaFile(path)
        for name, chrom, start, end in windows:
            for kind, s, e in (("core", start, end), ("extended", start - args.flank, end + args.flank)):
                result, status = extract(aligner, haplotype, reference, chrom, s, e, args)
                out = os.path.join(args.outdir, f"{name}.{kind}.hap{hap}.fa")
                with open(out, "w") as fh:
                    if result:
                        seq, ref_start, ref_end, ctg, lo, hi, strand = result
                        fh.write(f">{name}_hap{hap}\n")
                        for i in range(0, len(seq), 80):
                            fh.write(seq[i:i + 80] + "\n")
                        rows.append((name, kind, hap, status, ref_start, ref_end, ctg, lo, hi, strand, len(seq)))
                    else:
                        rows.append((name, kind, hap, status, "", "", "", "", "", "", 0))
                print(name, kind, f"hap{hap}", status)
    with open(args.status, "w") as out:
        out.write("window\tkind\thap\tstatus\tref_start\tref_end\tcontig\thap_start\thap_end\tstrand\tlength\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
