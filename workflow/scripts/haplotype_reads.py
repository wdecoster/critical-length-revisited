"""Reads of one window split by haplotype from WhatsHap's haplotag list (read, haplotype H1/H2/none, phase set):
H1 to haplotype 1, H2 to 2, untagged reads (homozygous stretches, no informative variant) to both. Reads shorter
than --min-length are dropped. As in a regional pipeline without block joining: reads of different phase blocks
keep their own block's H1/H2, so a haplotype can switch between blocks."""

import argparse
import gzip

import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", required=True, help="FASTQ of the window's reads (primary records)")
    parser.add_argument("--haplotag-list", required=True)
    parser.add_argument("--min-length", type=int, default=1000)
    parser.add_argument("--output1", required=True)
    parser.add_argument("--output2", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    tags = {}
    opener = gzip.open if args.haplotag_list.endswith(".gz") else open
    with opener(args.haplotag_list, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            tags[f[0]] = f[1]
    counts = {"H1": 0, "H2": 0, "none": 0}
    with gzip.open(args.output1, "wt") as out1, gzip.open(args.output2, "wt") as out2:
        with pysam.FastxFile(args.reads) as fq:
            for read in fq:
                if len(read.sequence) < args.min_length:
                    continue
                tag = tags.get(read.name, "none")
                tag = tag if tag in ("H1", "H2") else "none"
                counts[tag] += 1
                record = f"@{read.name}\n{read.sequence}\n+\n{read.quality or 'I' * len(read.sequence)}\n"
                if tag in ("H1", "none"):
                    out1.write(record)
                if tag in ("H2", "none"):
                    out2.write(record)
    print(counts)


if __name__ == "__main__":
    main()
