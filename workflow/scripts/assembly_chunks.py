"""Split the sequences of one or more FASTA files (with .fai) into chunks of whole sequences with similar
total length (greedy: longest sequence first, to the chunk with the least sequence so far).

SimLoRD caps a read at the longest sequence it was given, so a chunk without a sequence of at least the
longest simulated read length would silently produce shorter reads: that is an error (--min-longest).
"""

import argparse
import os
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", nargs="+", required=True)
    parser.add_argument("--chunks", type=int, required=True)
    parser.add_argument("--min-longest", type=int, required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    sequences = []
    for fasta in args.fasta:
        for line in open(fasta + ".fai"):
            name, length = line.split("\t")[:2]
            sequences.append((int(length), name, fasta))
    names = [s[1] for s in sequences]
    if len(set(names)) != len(names):
        raise ValueError("sequence names are not unique across the haplotypes")
    chunks = [[] for _ in range(args.chunks)]
    totals = [0] * args.chunks
    for length, name, fasta in sorted(sequences, reverse=True):
        i = totals.index(min(totals))
        chunks[i].append((name, fasta))
        totals[i] += length
    lengths = dict((s[1], s[0]) for s in sequences)
    os.makedirs(args.outdir, exist_ok=True)
    for i, chunk in enumerate(chunks):
        longest = max(lengths[name] for name, _ in chunk)
        if longest < args.min_longest:
            raise ValueError(f"chunk {i} has no sequence of at least {args.min_longest} bp: use fewer chunks")
        with open(os.path.join(args.outdir, f"{i:03d}.fa"), "w") as out:
            for fasta in dict.fromkeys(f for _, f in chunk):
                subprocess.run(["samtools", "faidx", fasta] + [n for n, f in chunk if f == fasta],
                               stdout=out, check=True)
        print(f"chunk {i:03d}: {len(chunk)} sequences, {totals[i]} bp, longest {longest} bp")


if __name__ == "__main__":
    main()
