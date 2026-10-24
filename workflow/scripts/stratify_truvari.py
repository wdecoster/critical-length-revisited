"""Recall, precision and F1 per SV stratum from truvari ga4gh output (base and comp VCFs with FORMAT/BD).

Strata: all; type (DEL, INS, INV, DUP, other); size bin (lower edges --size-bins, bp); context (first of
--contexts whose BED overlaps the variant, else 'other'); and type x size, type x context. A variant's span is
POS..POS+len(REF) (an insertion: its position); its size is |SVLEN|, else the REF/ALT length difference.
Recall counts base records (TP / (TP + FN)), precision comp records (TP / (TP + FP)).
"""

import argparse
import bisect
import gzip
from collections import defaultdict

import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True, help="ga4gh base VCF")
    parser.add_argument("--query", required=True, help="ga4gh comp VCF")
    parser.add_argument("--contexts", nargs="*", default=[], help="name=BED(.gz), in priority order")
    parser.add_argument("--size-bins", nargs="+", type=int, required=True)
    parser.add_argument("--chromosomes", nargs="*", help="only load context intervals on these")
    parser.add_argument("--labels", nargs="*", default=[], help="key=value columns added to every row")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


class Intervals:
    def __init__(self, path, chromosomes=None):
        raw = defaultdict(list)
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")):
                    continue
                chrom, start, end = line.split("\t")[:3]
                if chromosomes and chrom not in chromosomes:
                    continue
                raw[chrom].append((int(start), int(end)))
        self.starts, self.ends = {}, {}
        for chrom, intervals in raw.items():
            merged = []
            for start, end in sorted(intervals):
                if merged and start <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            self.starts[chrom] = [s for s, _ in merged]
            self.ends[chrom] = [e for _, e in merged]

    def overlaps(self, chrom, start, end):
        starts = self.starts.get(chrom)
        if not starts:
            return False
        i = bisect.bisect_left(starts, end) - 1
        return i >= 0 and self.ends[chrom][i] > start


def size_bin(size, bins):
    i = bisect.bisect_right(bins, size) - 1
    if i < 0:
        return f"<{bins[0]}"
    return f"{bins[i]}-{bins[i + 1]}" if i + 1 < len(bins) else f">={bins[i]}"


def describe(record, bins, contexts):
    svtype = record.info.get("SVTYPE") if "SVTYPE" in record.header.info else None
    alt = record.alts[0] if record.alts else ""
    svlen = record.info.get("SVLEN") if "SVLEN" in record.header.info else None
    if isinstance(svlen, tuple):
        svlen = svlen[0]
    diff = len(alt) - len(record.ref) if not alt.startswith("<") else None
    size = abs(svlen) if svlen is not None else abs(diff) if diff is not None else 0
    if svtype not in ("DEL", "INS", "INV", "DUP"):
        svtype = "INS" if diff and diff > 0 else "DEL" if diff and diff < 0 else "other"
    start = record.pos
    end = record.pos + len(record.ref) if svtype != "INS" else record.pos + 1
    if svtype in ("INV", "DUP") and record.stop > start:
        end = record.stop
    context = next((name for name, iv in contexts if iv.overlaps(record.chrom, start - 1, end)), "other")
    return svtype, size_bin(size, bins), context


def tally(path, bins, contexts, positive, negative):
    counts = defaultdict(lambda: [0, 0])
    with pysam.VariantFile(path) as vcf:
        sample = list(vcf.header.samples)[0]
        for record in vcf:
            decision = record.samples[sample].get("BD")
            if decision not in (positive, negative):
                continue
            svtype, size, context = describe(record, bins, contexts)
            for kind, stratum in (("all", "all"), ("type", svtype), ("size", size), ("context", context),
                                  ("type_size", f"{svtype}|{size}"), ("type_context", f"{svtype}|{context}")):
                counts[(kind, stratum)][0 if decision == positive else 1] += 1
    return counts


def main():
    args = parse_args()
    chromosomes = set(args.chromosomes) if args.chromosomes else None
    contexts = [(spec.split("=", 1)[0], Intervals(spec.split("=", 1)[1], chromosomes)) for spec in args.contexts]
    base = tally(args.truth, args.size_bins, contexts, "TP", "FN")
    comp = tally(args.query, args.size_bins, contexts, "TP", "FP")
    labels = dict(kv.split("=", 1) for kv in args.labels)
    with open(args.output, "w") as out:
        out.write("\t".join(list(labels) + ["kind", "stratum", "tp_base", "fn", "tp_comp", "fp", "recall",
                                             "precision", "f1"]) + "\n")
        for key in sorted(set(base) | set(comp)):
            tp_base, fn = base.get(key, [0, 0])
            tp_comp, fp = comp.get(key, [0, 0])
            recall = tp_base / (tp_base + fn) if tp_base + fn else float("nan")
            precision = tp_comp / (tp_comp + fp) if tp_comp + fp else float("nan")
            f1 = (2 * recall * precision / (recall + precision) if recall + precision > 0 else 0.0) \
                if tp_base + fn and tp_comp + fp else float("nan")
            out.write("\t".join(list(labels.values()) + [key[0], key[1], str(tp_base), str(fn), str(tp_comp),
                                                         str(fp), f"{recall:.4f}", f"{precision:.4f}",
                                                         f"{f1:.4f}"]) + "\n")


if __name__ == "__main__":
    main()
