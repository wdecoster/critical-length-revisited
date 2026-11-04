"""Oracle haplotype assignment: each read goes to the truth haplotype (window +/- flank) it aligns to best
(matching bases minus edit distance, over its primary alignments; minimap2 lr:hq via mappy). Ties, and reads that
align to neither, go to both haplotypes (as untagged reads in the realistic mode). Reads shorter than
--min-length are dropped.
Windows without a truth sequence for a haplotype (truth anchors not found: window_truth.py) have no oracle: both
read sets are written empty and the assignments file says why, so the oracle assemblies of that window are empty
and reported as such instead of failing the run."""

import argparse
import gzip

import mappy


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", required=True)
    parser.add_argument("--hap1", required=True)
    parser.add_argument("--hap2", required=True)
    parser.add_argument("--min-length", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--output1", required=True)
    parser.add_argument("--output2", required=True)
    parser.add_argument("--assignments", required=True)
    return parser.parse_args()


def score(aligner, seq):
    hits = [h for h in aligner.map(seq) if h.is_primary]
    return sum(h.mlen - h.NM for h in hits) if hits else None


def main():
    args = parse_args()
    missing = [path for path in (args.hap1, args.hap2) if not open(path).read(1)]
    if missing:  # no truth for a haplotype: no oracle possible for this window
        for output in (args.output1, args.output2):
            gzip.open(output, "wt").close()
        with open(args.assignments, "w") as table:
            table.write(f"# no oracle: empty truth {', '.join(missing)}\n")
        print(f"no oracle: empty truth {', '.join(missing)}")
        return
    aligners = [mappy.Aligner(path, preset="lr:hq", n_threads=args.threads) for path in (args.hap1, args.hap2)]
    counts = {"1": 0, "2": 0, "both": 0}
    with gzip.open(args.output1, "wt") as out1, gzip.open(args.output2, "wt") as out2, \
            open(args.assignments, "w") as table:
        for name, seq, qual in mappy.fastx_read(args.reads):
            if len(seq) < args.min_length:
                continue
            s1, s2 = (score(a, seq) for a in aligners)
            if s1 is not None and (s2 is None or s1 > s2):
                hap = "1"
            elif s2 is not None and (s1 is None or s2 > s1):
                hap = "2"
            else:
                hap = "both"
            counts[hap] += 1
            table.write(f"{name}\t{hap}\t{s1}\t{s2}\n")
            record = f"@{name}\n{seq}\n+\n{qual or 'I' * len(seq)}\n"
            if hap in ("1", "both"):
                out1.write(record)
            if hap in ("2", "both"):
                out2.write(record)
    print(counts)


if __name__ == "__main__":
    main()
