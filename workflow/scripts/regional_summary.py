"""Regional arm summary: all per-haplotype evaluations in one table, with the window annotations.

A failed assembly of one window is a result (an empty FASTA, n_contigs 0). An entire condition (read length, seed,
coverage, mode, assembler) without a single contig is not: that is a systematic failure (e.g. the assembler crashing
on every input of that condition), so the summary stops with the list of such conditions. Empty assemblies per
condition are printed either way.
"""
import argparse
import sys

import pandas as pd

CONDITION = ["length", "seed", "coverage", "mode", "assembler"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True)
    p.add_argument("--windows", required=True)
    p.add_argument("--evaluations", nargs="+", required=True)
    a = p.parse_args()
    windows = pd.read_table(a.windows).rename(columns={"name": "window"})
    table = pd.concat(pd.read_table(f) for f in a.evaluations).merge(windows, on="window", how="left")
    empty = table.assign(empty=table["n_contigs"] == 0).groupby(CONDITION)["empty"].agg(["sum", "size"])
    empty = empty.rename(columns={"sum": "empty", "size": "assemblies"})
    print(empty[empty["empty"] > 0].to_string(), file=sys.stderr)
    whole = empty[empty["empty"] == empty["assemblies"]]
    if len(whole):
        sys.exit(f"regional_summary: every assembly is empty in {len(whole)} condition(s), a systematic failure:\n"
                 f"{whole.to_string()}")
    table.to_csv(a.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
