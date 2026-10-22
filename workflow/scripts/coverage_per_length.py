"""Read set summary from read-length histograms: reads, bases, N50, and per cut length the coverage left after
cutting every read into pieces of that length (floor(read length / L) pieces; the remainder is lost).

Input histograms: two columns, read length and number of reads. The first rows of the output are summary
lines (length = "all"); coverage is bases / genome size.
"""

import argparse

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", nargs="+", required=True)
    parser.add_argument("--genome-size", type=float, required=True)
    parser.add_argument("--lengths", nargs="+", type=int, required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    hist = pd.concat(pd.read_table(h, header=None, names=["length", "reads"]) for h in args.histograms)
    hist = hist.groupby("length", as_index=False)["reads"].sum().sort_values("length")
    lengths, reads = hist["length"].to_numpy(float), hist["reads"].to_numpy(float)
    bases = lengths * reads
    total = bases.sum()
    n50 = lengths[np.searchsorted(np.cumsum(bases), total / 2)]
    rows = [{"length": "all", "reads": int(reads.sum()), "gigabases": round(total / 1e9, 2),
             "coverage": round(total / args.genome_size, 2), "n50": int(n50)}]
    for cut in args.lengths:
        pieces = np.floor(lengths / cut) * reads
        kept = (pieces * cut).sum()
        rows.append({"length": cut, "reads": int(pieces.sum()), "gigabases": round(kept / 1e9, 2),
                     "coverage": round(kept / args.genome_size, 2), "n50": cut})
    table = pd.DataFrame(rows)
    table.to_csv(args.output, sep="\t", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
