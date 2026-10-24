"""Names of the reads in a PAF (stdin) with a mapping overlapping the regions (BED), one per line (stdout).
Primary mappings only (tp:A:P), when the tag is present."""

import argparse
import bisect
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", required=True, help="BED, merged (non-overlapping) intervals")
    return parser.parse_args()


def main():
    args = parse_args()
    starts, ends = defaultdict(list), defaultdict(list)
    for line in open(args.regions):
        chrom, start, end = line.split("\t")[:3]
        starts[chrom].append(int(start))
        ends[chrom].append(int(end))
    seen = 0
    for line in sys.stdin:
        f = line.split("\t", 13)
        if "tp:A:S" in line or f[5] not in starts:
            continue
        chrom, start, end = f[5], int(f[7]), int(f[8])
        # the last region starting before the mapping's end overlaps it if it ends after the mapping's start
        i = bisect.bisect_left(starts[chrom], end) - 1
        if i >= 0 and ends[chrom][i] > start:
            sys.stdout.write(f[0] + "\n")
            seen += 1
    print(f"{seen} overlapping mappings", file=sys.stderr)


if __name__ == "__main__":
    main()
