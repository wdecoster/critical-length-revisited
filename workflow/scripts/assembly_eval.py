"""Evaluation of one assembled haplotype of a window against the two truth haplotypes (one row, TSV).

- Contigs: count, total length, largest, N50.
- Per truth haplotype (from minimap2 asm5 PAFs of the contigs against it, primary chains >= 1 kb): the fraction
  of the truth covered by all contigs, by the best single contig, and the identity of the aligned bases (matches
  / alignment block length), as QV. The truth haplotype the assembly matches is the one with the larger
  identity (then single-contig coverage, then coverage); the other gives the haplotype specificity.
- Complete: one contig covers >= --complete-fraction of the matched truth haplotype.
- Phasing: switch and Hamming error rates from yak trioeval with k-mers of the two truth haplotypes (W and H
  lines), when given.
"""

import argparse
import math
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--truth1", required=True)
    parser.add_argument("--truth2", required=True)
    parser.add_argument("--paf1", required=True)
    parser.add_argument("--paf2", required=True)
    parser.add_argument("--yak", help="yak trioeval output")
    parser.add_argument("--complete-fraction", type=float, default=0.95)
    parser.add_argument("--min-block", type=int, default=1000)
    parser.add_argument("--labels", nargs="*", default=[], help="key=value columns")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def fasta_lengths(path):
    lengths, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name:
                lengths[name] += len(line.strip())
    return lengths


def n50(lengths):
    total, running = sum(lengths), 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= total / 2:
            return length
    return 0


def union(intervals):
    total, end = 0, -1
    for s, e in sorted(intervals):
        if e <= end:
            continue
        total += e - max(s, end)
        end = e
    return total


def against(paf, truth_len, min_block):
    per_contig, matches, block = defaultdict(list), 0, 0
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or "tp:A:P" not in f[12:] or int(f[10]) < min_block:
                continue
            per_contig[f[0]].append((int(f[7]), int(f[8])))
            matches += int(f[9])
            block += int(f[10])
    if not truth_len:
        return None
    covered = union([iv for ivs in per_contig.values() for iv in ivs]) / truth_len
    single = max((union(ivs) for ivs in per_contig.values()), default=0) / truth_len
    identity = matches / block if block else float("nan")
    return {"covered": covered, "single": single, "identity": identity}


def qv(identity):
    if identity != identity:
        return float("nan")
    return 60.0 if identity >= 1 else min(60.0, -10 * math.log10(1 - identity))


def yak_rates(path):
    rates = {"switch_error": float("nan"), "hamming_error": float("nan")}
    if not path:
        return rates
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "W" and len(f) >= 4:
            rates["switch_error"] = float(f[3])
        elif f[0] == "H" and len(f) >= 4:
            rates["hamming_error"] = float(f[3])
    return rates


def main():
    args = parse_args()
    contigs = list(fasta_lengths(args.contigs).values())
    truth_lengths = [sum(fasta_lengths(p).values()) for p in (args.truth1, args.truth2)]
    row = dict(label.split("=", 1) for label in args.labels)
    row.update({
        "n_contigs": len(contigs), "assembly_length": sum(contigs), "largest_contig": max(contigs, default=0),
        "contig_n50": n50(contigs), "truth1_length": truth_lengths[0], "truth2_length": truth_lengths[1],
    })
    results = [against(p, t, args.min_block) for p, t in ((args.paf1, truth_lengths[0]), (args.paf2, truth_lengths[1]))]
    if all(r is not None for r in results):
        # the two haplotypes are nearly identical in coverage terms; identity separates them (then coverage)
        key = lambda i: (results[i]["identity"] if results[i]["identity"] == results[i]["identity"] else -1,
                         results[i]["single"], results[i]["covered"])
        best = max((0, 1), key=key)
        other = 1 - best
        b, o = results[best], results[other]
        row.update({
            "truth_ok": True, "matched_haplotype": best + 1,
            "covered": round(b["covered"], 4), "single_contig": round(b["single"], 4),
            "complete": b["single"] >= args.complete_fraction,
            "identity": round(b["identity"], 6), "qv": round(qv(b["identity"]), 2),
            "identity_other": round(o["identity"], 6),
        })
    else:
        row.update({"truth_ok": False, "matched_haplotype": "", "covered": "", "single_contig": "",
                    "complete": "", "identity": "", "qv": "", "identity_other": ""})
    row.update(yak_rates(args.yak))
    with open(args.output, "w") as out:
        out.write("\t".join(row) + "\n")
        out.write("\t".join(str(v) for v in row.values()) + "\n")


if __name__ == "__main__":
    main()
