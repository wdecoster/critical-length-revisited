#!/usr/bin/env python3
"""Cut long reads into pieces of a fixed length, keeping their real error profile.

Every read of at least --length bases is cut into floor(len / length) pieces of exactly --length bases. The
pieces are placed at random within the read: the leftover bases (len mod length) are spread randomly over the
gaps before, between and after the pieces, so no part of a read is systematically favoured and cut points fall
at random relative to any feature. Reads shorter than --length are dropped. Pieces are named <read>_<k>
(k = 1..n, in read order), so no tool merges evidence across pieces of one read.

With --target-bases, pieces are kept at random with the probability that gives about that many bases in
total (a first pass over the input counts the piece bases). Deterministic for a given --seed.

Input: FASTQ (plain or gzipped), or '-' for stdin (no --target-bases then). Output: FASTQ, gzipped if the name
ends in .gz, or '-' for stdout.
"""

import argparse
import gzip
import sys

import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="FASTQ(.gz), or - for stdin")
    parser.add_argument("--length", type=int, required=True, help="piece length (bp)")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--target-bases", type=float, help="keep pieces at random to about this many bases")
    parser.add_argument("--output", default="-", help="FASTQ(.gz), or - for stdout (default)")
    parser.add_argument("--stats", help="write a summary (TSV) here")
    return parser.parse_args()


def open_text(path, mode):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", compresslevel=1) if "w" in mode else gzip.open(path, mode + "t")
    return open(path, mode)


def fastq(handle):
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline().rstrip("\n")
        handle.readline()
        quality = handle.readline().rstrip("\n")
        yield header[1:].split()[0], sequence, quality


def piece_starts(read_length, length, rng):
    """Start positions of floor(read_length / length) pieces with the leftover spread randomly over the gaps"""
    n = read_length // length
    leftover = read_length - n * length
    # n sorted uniform offsets in [0, leftover]: piece k starts after k pieces and offset k of leftover bases
    offsets = np.sort(rng.integers(0, leftover + 1, size=n))
    return np.arange(n) * length + offsets


def main():
    args = parse_args()
    keep = 1.0
    if args.target_bases:
        if args.input == "-":
            sys.exit("--target-bases needs a file input (two passes)")
        with open_text(args.input, "r") as handle:
            total = sum((len(seq) // args.length) * args.length for _, seq, _ in fastq(handle))
        keep = min(1.0, args.target_bases / total) if total else 1.0
    rng = np.random.default_rng(args.seed)
    counts = {"reads": 0, "reads_cut": 0, "bases_in": 0, "pieces": 0, "pieces_kept": 0, "bases_out": 0}
    with open_text(args.input, "r") as handle, open_text(args.output, "w") as out:
        for name, sequence, quality in fastq(handle):
            counts["reads"] += 1
            counts["bases_in"] += len(sequence)
            if len(sequence) < args.length:
                continue
            counts["reads_cut"] += 1
            starts = piece_starts(len(sequence), args.length, rng)
            # one draw per piece, also when all are kept, so the cuts do not depend on --target-bases
            kept = rng.random(len(starts)) < keep
            counts["pieces"] += len(starts)
            for k, (start, keep_piece) in enumerate(zip(starts, kept), start=1):
                if not keep_piece:
                    continue
                end = start + args.length
                out.write(f"@{name}_{k}\n{sequence[start:end]}\n+\n{quality[start:end]}\n")
                counts["pieces_kept"] += 1
                counts["bases_out"] += args.length
    counts["keep_probability"] = round(keep, 6)
    summary = "\t".join(counts) + "\n" + "\t".join(str(v) for v in counts.values()) + "\n"
    if args.stats:
        open(args.stats, "w").write(summary)
    print(summary, file=sys.stderr, end="")


if __name__ == "__main__":
    main()
