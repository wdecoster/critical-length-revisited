"""Runs of N of at least --min-length bp in a FASTA (plain or gzip/bgzip), as BED (stdout)."""

import argparse
import gzip
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta")
    parser.add_argument("--min-length", type=int, default=100)
    return parser.parse_args()


def open_fasta(path):
    with open(path, "rb") as fh:
        compressed = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if compressed else open(path)


def records(path):
    name, chunks = None, []
    with open_fasta(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.rstrip())
    if name:
        yield name, "".join(chunks)


def main():
    args = parse_args()
    pattern = re.compile(f"[Nn]{{{args.min_length},}}")
    for name, seq in records(args.fasta):
        for m in pattern.finditer(seq):
            sys.stdout.write(f"{name}\t{m.start()}\t{m.end()}\n")


if __name__ == "__main__":
    main()
