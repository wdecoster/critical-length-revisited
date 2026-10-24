"""Per truth SV inside the regions: the number of reads (primary or supplementary alignments, MAPQ >=
--min-mapq) spanning it with --flank bases on both sides. Mechanistic view: recall against spanning reads.

Output (TSV.gz): chrom, pos, end, svtype, size, spanning.
"""

import argparse
import bisect
import gzip
from collections import defaultdict

import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--regions", required=True, help="BED (merged)")
    parser.add_argument("--flank", type=int, default=1000)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-size", type=int, default=50)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    starts, ends = defaultdict(list), defaultdict(list)
    for line in open(args.regions):
        chrom, start, end = line.split("\t")[:3]
        starts[chrom].append(int(start))
        ends[chrom].append(int(end))
    n = 0
    with pysam.VariantFile(args.truth) as vcf, pysam.AlignmentFile(args.bam) as bam, \
            gzip.open(args.output, "wt") as out:
        out.write("chrom\tpos\tend\tsvtype\tsize\tspanning\n")
        for chrom in starts:
            if chrom not in vcf.header.contigs:
                continue
            for record in vcf.fetch(chrom):
                alt = record.alts[0] if record.alts else ""
                if alt.startswith("<") or alt == "*":
                    continue
                diff = len(alt) - len(record.ref)
                if abs(diff) < args.min_size:
                    continue
                start, end = record.pos, record.pos + len(record.ref)
                i = bisect.bisect_right(starts[chrom], start) - 1
                if i < 0 or ends[chrom][i] < end:
                    continue  # not (entirely) inside the regions
                left, right = max(0, start - args.flank), end + args.flank
                spanning = 0
                for read in bam.fetch(chrom, left, left + 1):
                    if read.is_secondary or read.is_unmapped or read.mapping_quality < args.min_mapq:
                        continue
                    if read.reference_start <= left and read.reference_end >= right:
                        spanning += 1
                out.write(f"{chrom}\t{start}\t{end}\t{'INS' if diff > 0 else 'DEL'}\t{abs(diff)}\t{spanning}\n")
                n += 1
    print(f"{n} truth SVs")


if __name__ == "__main__":
    main()
