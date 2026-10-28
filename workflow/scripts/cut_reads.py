#!/usr/bin/env python3
"""Cut long reads into pieces of a fixed length, or of lengths drawn from a library's length distribution, keeping
their real error profile.

Every read of at least --length bases is cut into floor(len / length) pieces of exactly --length bases. The
pieces are placed at random within the read: the leftover bases (len mod length) are spread randomly over the
gaps before, between and after the pieces, so no part of a read is systematically favoured and cut points fall
at random relative to any feature. Reads shorter than --length are dropped. Pieces are named <read>_<k>
(k = 1..n, in read order), so no tool merges evidence across pieces of one read.

Distribution mode (--distribution and --n50 instead of --length): piece lengths follow a real library's read
length distribution (read_lengths histograms: length<TAB>count), rescaled so its N50 is --n50. Each read is filled
with pieces of lengths drawn from it until the next one does not fit; the leftover is spread randomly over the
gaps, as in fixed mode. Long pieces only fit in long reads, so drawing from the target distribution itself would
give too few of them: the drawing weights are calibrated first (--calibration-rounds), by filling the reads'
own lengths (a first pass) and reweighting each length bin by target / achieved share, then rescaling the bins
until the N50 of the pieces is within 1% of --n50 (dropping pieces below --min-piece and the reads' own lengths
both shift it). The stats report the target and achieved N50.

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
    parser.add_argument("--length", type=int, help="piece length (bp), fixed mode")
    parser.add_argument("--distribution", nargs="+", help="distribution mode: read length histogram(s) of a library")
    parser.add_argument("--n50", type=int, help="distribution mode: the library's lengths rescaled to this N50")
    parser.add_argument("--min-piece", type=int, default=200, help="distribution mode: shortest piece")
    parser.add_argument("--calibration-rounds", type=int, default=8)
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


def n50(lengths, counts=None):
    lengths = np.asarray(lengths, dtype=float)
    counts = np.ones(len(lengths)) if counts is None else np.asarray(counts, dtype=float)
    order = np.argsort(lengths)[::-1]
    bases = np.cumsum(lengths[order] * counts[order])
    return float(lengths[order][np.searchsorted(bases, bases[-1] / 2)]) if len(bases) else 0.0


class LengthDistribution:
    """A library's read length distribution rescaled to an N50, in log-spaced bins; lengths are drawn by bin
    (calibrated weights), then uniformly within it"""

    def __init__(self, histograms, target_n50, min_piece, bins=240):
        lengths, counts = [], []
        for path in histograms:
            table = np.loadtxt(path, dtype=float, ndmin=2)
            lengths.append(table[:, 0])
            counts.append(table[:, 1])
        lengths, counts = np.concatenate(lengths), np.concatenate(counts)
        keep = lengths > 0
        lengths, counts = lengths[keep], counts[keep]
        self.source_n50 = n50(lengths, counts)
        scaled = lengths * target_n50 / self.source_n50
        keep = scaled >= min_piece
        scaled, counts = scaled[keep], counts[keep]
        self.edges = np.geomspace(min_piece, scaled.max() * 1.0001, bins + 1)
        target, _ = np.histogram(scaled, bins=self.edges, weights=counts)
        self.target = target / target.sum()  # share of pieces per bin
        self.weights = self.target.copy()
        self.target_n50 = target_n50

    def draw(self, rng, n):
        b = rng.choice(len(self.weights), size=n, p=self.weights / self.weights.sum())
        return np.floor(rng.uniform(self.edges[b], self.edges[b + 1])).astype(np.int64)

    def calibrate(self, read_lengths, rng, rounds):
        """Reweight bins until the pieces cut from these reads follow the target shares"""
        for _ in range(rounds):
            pieces = np.concatenate([fill(r, self, rng)[1] for r in read_lengths] or [np.zeros(0, np.int64)])
            achieved, _ = np.histogram(pieces, bins=self.edges)
            achieved = achieved / max(1, achieved.sum())
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(achieved > 0, self.target / achieved, 1.0)
            self.weights = np.where(self.target > 0, self.weights * np.clip(ratio, 0.2, 5), 0)
        # the pieces the final weights give
        return np.concatenate([fill(r, self, rng)[1] for r in read_lengths] or [np.zeros(0, np.int64)])


def fill(read_length, dist, rng, batch=16):
    """Piece lengths filling one read (drawn until the next does not fit), and their starts with the leftover
    spread randomly over the gaps"""
    lengths, used = [], 0
    while True:
        for length in dist.draw(rng, batch):
            if used + length > read_length:
                lengths = np.array(lengths, dtype=np.int64)
                leftover = read_length - used
                offsets = np.sort(rng.integers(0, leftover + 1, size=len(lengths)))
                starts = np.concatenate([[0], np.cumsum(lengths)[:-1]]).astype(np.int64) + offsets
                return starts, lengths
            lengths.append(int(length))
            used += int(length)


def main_distribution(args):
    """Distribution mode: calibrate the drawing weights on the reads' lengths, then cut"""
    with open_text(args.input, "r") as handle:
        read_lengths = np.array([len(seq) for _, seq, _ in fastq(handle)], dtype=np.int64)
    dist = LengthDistribution(args.distribution, args.n50, args.min_piece)
    calibration_rng = np.random.default_rng(args.seed + 1000)
    sample = calibration_rng.choice(read_lengths, size=min(len(read_lengths), 20000), replace=False)
    # shares first, then the N50 itself: pieces below --min-piece are dropped and the reads cap the longest, both
    # of which shift the N50; rescale the bins until it is within 1% (or 4 tries)
    for _ in range(4):
        pieces = dist.calibrate(sample, calibration_rng, args.calibration_rounds)
        achieved = n50(pieces)
        if not achieved or abs(achieved / args.n50 - 1) < 0.01:
            break
        dist.edges = np.maximum(dist.edges * args.n50 / achieved, args.min_piece)
    # cut: one stream for the cuts, another for keeping (the cuts do not depend on --target-bases)
    cut_rng, keep_rng = np.random.default_rng(args.seed), np.random.default_rng(args.seed + 1)
    plan = [fill(r, dist, cut_rng) for r in read_lengths]
    total = sum(int(lengths.sum()) for _, lengths in plan)
    keep = min(1.0, args.target_bases / total) if args.target_bases and total else 1.0
    counts = {"reads": 0, "reads_cut": 0, "bases_in": 0, "pieces": 0, "pieces_kept": 0, "bases_out": 0}
    kept_lengths = []
    with open_text(args.input, "r") as handle, open_text(args.output, "w") as out:
        for (name, sequence, quality), (starts, lengths) in zip(fastq(handle), plan):
            counts["reads"] += 1
            counts["bases_in"] += len(sequence)
            if not len(lengths):
                continue
            counts["reads_cut"] += 1
            counts["pieces"] += len(lengths)
            kept = keep_rng.random(len(lengths)) < keep
            for k, (start, length, keep_piece) in enumerate(zip(starts, lengths, kept), start=1):
                if not keep_piece:
                    continue
                end = start + length
                out.write(f"@{name}_{k}\n{sequence[start:end]}\n+\n{quality[start:end]}\n")
                counts["pieces_kept"] += 1
                counts["bases_out"] += int(length)
                kept_lengths.append(int(length))
    counts["keep_probability"] = round(keep, 6)
    counts["source_n50"] = round(dist.source_n50)
    counts["target_n50"] = args.n50
    counts["achieved_n50"] = round(n50(kept_lengths)) if kept_lengths else 0
    counts["mean_piece"] = round(float(np.mean(kept_lengths))) if kept_lengths else 0
    return counts


def main():
    args = parse_args()
    if args.distribution:
        if args.n50 is None or args.input == "-":
            sys.exit("--distribution needs --n50 and a file input (two passes)")
        counts = main_distribution(args)
        summary = "\t".join(counts) + "\n" + "\t".join(str(v) for v in counts.values()) + "\n"
        if args.stats:
            open(args.stats, "w").write(summary)
        print(summary, file=sys.stderr, end="")
        return
    if args.length is None:
        sys.exit("give --length (fixed mode) or --distribution and --n50")
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
