"""SVs from paftools call VCF output (stdin -> stdout): records whose REF and ALT lengths differ by at least
--min-length become INS or DEL with SVTYPE, SVLEN and END, as the 2020 vcf-to-sv-vcf.py did. The genotype is
1/1: each input is one assembly haplotype, and haplotypes are combined later.
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-length", type=int, default=50)
    parser.add_argument("--sample", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    out = sys.stdout
    out.write("##fileformat=VCFv4.2\n")
    out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="INS or DEL, from the allele lengths">\n')
    out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="ALT minus REF length">\n')
    out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position on the reference">\n')
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n")
    kept = total = 0
    for line in sys.stdin:
        if line.startswith("#"):
            continue
        total += 1
        chrom, pos, _, ref, alt = line.split("\t")[:5]
        alt = alt.split(",")[0]
        svlen = len(alt) - len(ref)
        if abs(svlen) < args.min_length:
            continue
        kept += 1
        svtype = "INS" if svlen > 0 else "DEL"
        end = int(pos) + len(ref) - 1
        out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tSVTYPE={svtype};SVLEN={svlen};END={end}\tGT\t1/1\n")
    print(f"{total} variants, {kept} SVs of at least {args.min_length} bp", file=sys.stderr)


if __name__ == "__main__":
    main()
