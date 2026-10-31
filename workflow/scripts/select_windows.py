"""Windows for the regional assembly arm: a seeded random draw per class from a grid of windows, plus named ones.

Grid: windows of --size bp every --size/2 bp on the autosomes, at least 1 Mb from a chromosome end, without
reference gaps (N runs >= 100 bp) within the window +/- 100 kb. Classes (first match wins, in this order):
- segdup: >= 20 kb of segmental duplication with >= 99% identity (UCSC genomicSuperDups), >= 30% in the
  GIAB benchmark regions (keeps centromeres and satellites out);
- large_insertion: a benchmark insertion >= 5 kb, >= 80% benchmark, <= 10% segmental duplication;
- tandem_repeat: the top 5% of windows by bases in tandem repeats of >= 1 kb, >= 50% benchmark,
  <= 10% segmental duplication;
- plain: >= 90% benchmark, no segmental duplication, tandem repeat bases (>= 1 kb repeats) at most the median.
Chosen windows do not overlap each other. Output: name, chrom, start, end, class (TSV with header) and a BED.
"""

import argparse
import bisect
import gzip
import random
import re
from collections import defaultdict

import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--gaps", required=True, help="BED of reference N runs")
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--segdups", required=True, help="GIAB segdups BED (gz)")
    parser.add_argument("--superdups", required=True, help="UCSC genomicSuperDups.txt.gz")
    parser.add_argument("--tandem-repeats", required=True, help="GIAB tandem repeats BED (gz)")
    parser.add_argument("--truth", required=True, help="GIAB SV benchmark VCF")
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--classes", nargs="+", required=True, help="class=count")
    parser.add_argument("--named", nargs="*", default=[], help="name=chrom:start-end")
    parser.add_argument("--output", required=True)
    parser.add_argument("--bed", required=True)
    return parser.parse_args()


AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def read_bed(path, min_length=0):
    opener = gzip.open if path.endswith(".gz") else open
    intervals = defaultdict(list)
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            start, end = int(start), int(end)
            if end - start >= min_length:
                intervals[chrom].append((start, end))
    return {chrom: merge(iv) for chrom, iv in intervals.items()}


def merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


class Coverage:
    """Bases of merged intervals within a query interval"""

    def __init__(self, intervals):
        self.starts = {c: [s for s, _ in iv] for c, iv in intervals.items()}
        self.intervals = intervals

    def bases(self, chrom, start, end):
        if chrom not in self.intervals:
            return 0
        iv = self.intervals[chrom]
        i = max(0, bisect.bisect_right(self.starts[chrom], start) - 1)
        total = 0
        while i < len(iv) and iv[i][0] < end:
            total += max(0, min(end, iv[i][1]) - max(start, iv[i][0]))
            i += 1
        return total


def main():
    args = parse_args()
    lengths = {c: int(l) for c, l in (line.split("\t")[:2] for line in open(args.fai))}
    gaps = Coverage(read_bed(args.gaps))
    benchmark = Coverage(read_bed(args.benchmark))
    segdups = Coverage(read_bed(args.segdups))
    tandem = Coverage(read_bed(args.tandem_repeats, min_length=1000))
    identical = defaultdict(list)
    with gzip.open(args.superdups, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if float(f[26]) >= 0.99:
                identical[f[1]].append((int(f[2]), int(f[3])))
    identical = Coverage({c: merge(iv) for c, iv in identical.items()})
    insertions = defaultdict(list)
    with pysam.VariantFile(args.truth) as vcf:
        for record in vcf:
            alt = record.alts[0] if record.alts else ""
            if len(alt) - len(record.ref) >= 5000:
                insertions[record.chrom].append(record.pos)
    for positions in insertions.values():
        positions.sort()

    size = args.size
    grid = []
    for chrom in AUTOSOMES:
        for start in range(1_000_000, lengths[chrom] - 1_000_000 - size, size // 2):
            end = start + size
            if gaps.bases(chrom, start - 100_000, end + 100_000):
                continue
            pos = insertions.get(chrom, [])
            i = bisect.bisect_left(pos, start)
            grid.append({
                "chrom": chrom, "start": start, "end": end,
                "benchmark": benchmark.bases(chrom, start, end) / size,
                "segdup": segdups.bases(chrom, start, end) / size,
                "identical_sd": identical.bases(chrom, start, end),
                "tandem": tandem.bases(chrom, start, end),
                "insertion": i < len(pos) and pos[i] < end,
            })
    ranked = sorted(w["tandem"] for w in grid)
    median, top = ranked[len(ranked) // 2], ranked[int(len(ranked) * 0.95)]
    for w in grid:
        if w["identical_sd"] >= 20_000 and w["benchmark"] >= 0.3:
            w["class"] = "segdup"
        elif w["insertion"] and w["benchmark"] >= 0.8 and w["segdup"] <= 0.1:
            w["class"] = "large_insertion"
        elif w["tandem"] >= top and w["benchmark"] >= 0.5 and w["segdup"] <= 0.1:
            w["class"] = "tandem_repeat"
        elif w["benchmark"] >= 0.9 and w["segdup"] == 0 and w["tandem"] <= median:
            w["class"] = "plain"
        else:
            w["class"] = None

    chosen = []

    def free(chrom, start, end):
        return all(c != chrom or end <= s or start >= e for _, c, s, e, _ in chosen)

    for named in args.named:
        name, region = named.split("=")
        chrom, start, end = re.match(r"(\w+):(\d+)-(\d+)$", region).groups()
        chosen.append((name, chrom, int(start), int(end), "named"))
    rng = random.Random(args.seed)
    for spec in args.classes:
        cls, count = spec.split("=")
        candidates = [w for w in grid if w["class"] == cls]
        rng.shuffle(candidates)
        picked = 0
        for w in candidates:
            if picked == int(count):
                break
            if free(w["chrom"], w["start"], w["end"]):
                picked += 1
                chosen.append((f"{cls}_{picked}", w["chrom"], w["start"], w["end"], cls))
        print(f"{cls}: {len(candidates)} candidates, {picked} chosen")
    with open(args.output, "w") as out, open(args.bed, "w") as bed:
        out.write("name\tchrom\tstart\tend\tclass\n")
        for name, chrom, start, end, cls in chosen:
            out.write(f"{name}\t{chrom}\t{start}\t{end}\t{cls}\n")
            bed.write(f"{chrom}\t{start}\t{end}\t{name}\n")


if __name__ == "__main__":
    main()
