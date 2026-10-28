#!/usr/bin/env python3
"""All strata of the length distribution runs in one table, with per run the library shape, the target N50 and the
N50 and mean length the pieces achieved (cut_reads.py distribution mode stats). Performance is best plotted
against the achieved N50: at the longest targets the parent reads limit how long the pieces can be."""

import argparse
import os
import re

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strata", nargs="+", required=True)
    parser.add_argument("--pieces", nargs="+", required=True, help="pieces.stats.tsv of each run")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def run_key(path):
    """(read_set, length label, seed) from a run directory .../<read_set>/L<length>/s<seed>/<file>"""
    parts = os.path.normpath(path).split(os.sep)
    return parts[-4], parts[-3][1:], int(parts[-2][1:])


def main():
    args = parse_args()
    table = pd.concat(pd.read_table(f, dtype={"length": str}) for f in args.strata)
    pieces = []
    for path in args.pieces:
        stats = pd.read_table(path).iloc[0]
        read_set, length, seed = run_key(path)
        pieces.append({"read_set": read_set, "length": length, "seed": seed,
                       "achieved_n50": stats["achieved_n50"], "mean_piece": stats["mean_piece"],
                       "pieces_kept": stats["pieces_kept"], "keep_probability": stats["keep_probability"]})
    pieces = pd.DataFrame(pieces)
    match = table["length"].str.extract(r"N(\d+)([a-z]+)$")
    table.insert(1, "shape", match[1])
    table.insert(2, "target_n50", match[0].astype(int))
    table["seed"] = table["seed"].astype(int)
    table = table.merge(pieces, on=["read_set", "length", "seed"], how="left")
    table.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
