"""One diploid sample from a SURVIVOR merge of two haplotypes' SVs: 1/1 where both haplotypes have the SV,
0/1 where one has it (the 2020 truth for a phased assembly, as surpyvor haplomerge). VCF to stdout.
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="SURVIVOR merge output with two samples")
    parser.add_argument("--sample", required=True)
    return parser.parse_args()


def has_variant(field):
    gt = field.split(":")[0].replace("|", "/")
    return any(allele not in ("0", ".") for allele in gt.split("/"))


def main():
    args = parse_args()
    out = sys.stdout
    counts = {"0/1": 0, "1/1": 0}
    for line in open(args.vcf):
        if line.startswith("##"):
            if not line.startswith("##FORMAT") or "ID=GT," in line:
                out.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            out.write("\t".join(fields[:9] + [args.sample]) + "\n")
            continue
        n = sum(has_variant(f) for f in fields[9:11])
        if n == 0:
            continue
        gt = "1/1" if n == 2 else "0/1"
        counts[gt] += 1
        out.write("\t".join(fields[:8] + ["GT", gt]) + "\n")
    print(f"{counts['0/1']} heterozygous, {counts['1/1']} homozygous SVs", file=sys.stderr)


if __name__ == "__main__":
    main()
