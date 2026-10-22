"""Precision, recall and F1 of SV calls against a truth set, computed as surpyvor prf did in 2020.

Both VCFs get every "DUP" replaced by "INS" (plain text, as surpyvor normalize_vcf), then SURVIVOR merge
(distance, 1 caller, same type, strand ignored, no distance estimation, min. length). A merged site counts for
a VCF when that sample's genotype there is heterozygous or homozygous alternative (missing and 0/0 do not).
Sites are identified by chromosome, start and SVTYPE, as sets: TP = sites of both, precision = TP / test
sites, recall = TP / truth sites.
"""

import argparse
import os
import subprocess
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--calls", required=True)
    parser.add_argument("--distance", type=int, default=500)
    parser.add_argument("--min-length", type=int, default=50)
    parser.add_argument("--labels", nargs="*", default=[], help="key=value columns added to the output")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def is_variant(field):
    alleles = field.split(":")[0].replace("|", "/").split("/")
    if "." in alleles:
        return False
    return any(a != "0" for a in alleles)


def site_sets(merged):
    sets = [set(), set()]
    for line in open(merged):
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) for kv in fields[7].split(";") if "=" in kv)
        site = f"{fields[0]}:{int(fields[1]) - 1}-{info.get('SVTYPE')}"
        for index, field in enumerate(fields[9:11]):
            if is_variant(field):
                sets[index].add(site)
    return sets


def main():
    args = parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        paths = []
        for name, vcf in (("truth", args.truth), ("calls", args.calls)):
            path = os.path.join(tmp, f"{name}.vcf")
            with open(vcf) as source, open(path, "w") as out:
                for line in source:
                    out.write(line.replace("DUP", "INS"))
            paths.append(path)
        fofn = os.path.join(tmp, "files.txt")
        with open(fofn, "w") as out:
            out.write("\n".join(paths) + "\n")
        merged = os.path.join(tmp, "merged.vcf")
        subprocess.run(["SURVIVOR", "merge", fofn, str(args.distance), "1", "1", "-1", "-1", str(args.min_length),
                        merged], check=True, stdout=subprocess.DEVNULL)
        truth, test = site_sets(merged)
    tp = len(truth & test)
    precision = tp / len(test) if test else float("nan")
    recall = tp / len(truth) if truth else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else float("nan")
    labels = dict(kv.split("=", 1) for kv in args.labels)
    columns = {**labels, "truth_sites": len(truth), "test_sites": len(test), "tp": tp,
               "precision": round(precision, 4), "recall": round(recall, 4), "f1": round(f1, 4)}
    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n" + "\t".join(map(str, columns.values())) + "\n")
    print(columns)


if __name__ == "__main__":
    main()
