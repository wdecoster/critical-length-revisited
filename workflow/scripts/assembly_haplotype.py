"""Write one haplotype of an assembly as plain FASTA, streamed from a URL (gzipped or not).

With --report and --unit, only the sequences of that assembly unit (e.g. phase0 of a FALCON-Phase GenBank
assembly, where both phases are in one FASTA) are kept, selected by GenBank accession from the NCBI assembly
report. Names are the first word of the header.
"""

import argparse
import gzip
import io
import urllib.request


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--report", help="NCBI assembly report (URL)")
    parser.add_argument("--unit", help="assembly unit to keep, e.g. phase0")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def unit_accessions(report_url, unit):
    """GenBank accessions of one assembly unit (column 8 of the report)"""
    accessions = set()
    with urllib.request.urlopen(report_url) as handle:
        for line in io.TextIOWrapper(handle, encoding="utf-8"):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[7] == unit:
                accessions.add(fields[4])
    if not accessions:
        raise ValueError(f"no sequences of unit {unit} in {report_url}")
    return accessions


def main():
    args = parse_args()
    keep = unit_accessions(args.report, args.unit) if args.unit else None
    kept = written = 0
    with urllib.request.urlopen(args.url) as handle, open(args.output, "w") as out:
        stream = gzip.GzipFile(fileobj=handle) if args.url.endswith(".gz") else handle
        write = False
        for line in io.TextIOWrapper(stream, encoding="ascii"):
            if line.startswith(">"):
                name = line[1:].split()[0]
                write = keep is None or name in keep
                if write:
                    kept += 1
                    out.write(f">{name}\n")
            elif write:
                out.write(line)
                written += len(line) - 1
    print(f"{kept} sequences, {written} bp")
    if keep is not None and kept != len(keep):
        raise ValueError(f"expected {len(keep)} sequences of unit {args.unit}, found {kept}")


if __name__ == "__main__":
    main()
