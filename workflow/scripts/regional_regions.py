"""Regions of the regional arm.

- selection: each window +/- --margin, plus the paralog partners elsewhere of the segmental duplications
  (UCSC genomicSuperDups) overlapping a window +/- --flank, merged. Full-length reads with a mapping here are
  cut, so pieces that could mismap into a window's duplications are present.
- evaluation: the windows within the GIAB benchmark regions (where coverage is measured).
- sizes: selection and evaluation size (bp), for coverage targets.
"""

import argparse
import bisect
import gzip
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", required=True, help="TSV: name chrom start end class")
    parser.add_argument("--fai", required=True)
    parser.add_argument("--superdups", required=True)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--margin", type=int, required=True)
    parser.add_argument("--flank", type=int, required=True)
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
    lengths = {c: int(l) for c, l in (line.split("\t")[:2] for line in open(args.fai))}
    windows = []
    with open(args.windows) as fh:
        next(fh)
        for line in fh:
            _, chrom, start, end, _ = line.rstrip("\n").split("\t")
            windows.append((chrom, int(start), int(end)))
    selection = [(c, max(0, s - args.margin), min(lengths[c], e + args.margin)) for c, s, e in windows]
    partners = 0
    with gzip.open(args.superdups, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, start, end = f[1], int(f[2]), int(f[3])
            if f[7] not in lengths:
                continue
            for c, s, e in windows:
                if chrom == c and start < e + args.flank and end > s - args.flank:
                    selection.append((f[7], int(f[8]), int(f[9])))
                    partners += 1
                    break
    selection_bp = write(args.selection, merge(selection))
    benchmark = defaultdict(list)
    for line in open(args.benchmark):
        f = line.split("\t")
        benchmark[f[0]].append((int(f[1]), int(f[2])))
    evaluation = []
    for c, s, e in windows:
        iv = sorted(benchmark.get(c, []))
        i = max(0, bisect.bisect_right([x[0] for x in iv], s) - 1)
        while i < len(iv) and iv[i][0] < e:
            if iv[i][1] > s:
                evaluation.append((c, max(s, iv[i][0]), min(e, iv[i][1])))
            i += 1
    evaluation_bp = write(args.evaluation, merge(evaluation))
    with open(args.sizes, "w") as out:
        out.write(f"selection\t{selection_bp}\nevaluation\t{evaluation_bp}\nsd_partner_intervals\t{partners}\n")
    print(f"selection {selection_bp:,} bp ({partners} SD partner intervals), evaluation {evaluation_bp:,} bp")


if __name__ == "__main__":
    main()
