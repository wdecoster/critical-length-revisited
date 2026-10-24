"""Pilot regions.

- selection: the pilot chromosomes, plus (--sd-partners) the paralog partners elsewhere of the segmental
  duplications on them (UCSC genomicSuperDups: otherChrom, otherStart, otherEnd), merged. Parent reads with an
  alignment here are used, so reads that could mismap into the pilot's SDs are present.
- evaluation: the GIAB benchmark regions on the pilot chromosomes.
- sizes: selection and evaluation size (bp), for coverage targets.
"""

import argparse
import gzip


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--superdups", required=True)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    parser.add_argument("--sd-partners", action="store_true")
    parser.add_argument("--selection", required=True)
    parser.add_argument("--evaluation", required=True)
    parser.add_argument("--sizes", required=True)
    return parser.parse_args()


def merge(intervals):
    merged = []
    for chrom, start, end in sorted(intervals):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    return merged


def write(path, intervals):
    with open(path, "w") as out:
        for chrom, start, end in intervals:
            out.write(f"{chrom}\t{start}\t{end}\n")
    return sum(end - start for _, start, end in intervals)


def main():
    args = parse_args()
    pilot = set(args.chromosomes)
    lengths = dict(line.split("\t")[:2] for line in open(args.fai))
    selection = [(c, 0, int(lengths[c])) for c in args.chromosomes]
    partners = 0
    if args.sd_partners:
        with gzip.open(args.superdups, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                chrom, other, start, end = f[1], f[7], int(f[8]), int(f[9])
                if chrom in pilot and other not in pilot and other in lengths:
                    selection.append((other, start, end))
                    partners += 1
    selection_bp = write(args.selection, merge(selection))
    evaluation = []
    for line in open(args.benchmark):
        f = line.split("\t")
        if f[0] in pilot:
            evaluation.append((f[0], int(f[1]), int(f[2])))
    evaluation_bp = write(args.evaluation, merge(evaluation))
    with open(args.sizes, "w") as out:
        out.write(f"selection\t{selection_bp}\nevaluation\t{evaluation_bp}\nsd_partner_intervals\t{partners}\n")
    print(f"selection {selection_bp:,} bp ({partners} SD partner intervals), evaluation {evaluation_bp:,} bp")


if __name__ == "__main__":
    main()
