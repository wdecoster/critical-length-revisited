# Data inventory (milestone 1): read-length distributions of the HG002 read sets, and the coverage left after
# cutting reads into pieces of a fixed length. Read sets marked download are kept for Arm A.


rule download_reads:
    output:
        os.path.join(RES, "reads", "{read_set}", "{name}"),
    params:
        url=lambda wc: next(u for u in READ_SETS[wc.read_set]["urls"] if os.path.basename(u) == wc.name),
    log:
        os.path.join(LOGS, "inventory", "download.{read_set}.{name}.log"),
    shell:
        # resumable; ONT open data has no published checksums, so the BAM is checked for truncation
        """
        curl -sSfL -C - --retry 5 -o {output} {params.url} &> {log}
        """


rule read_lengths:
    """Histogram of read lengths (primary records and unmapped reads, so every read once) of one source"""
    input:
        lambda wc: [] if read_set_file(wc.read_set, wc.i).startswith("http") else read_set_file(wc.read_set, wc.i),
    output:
        os.path.join(OUT, "inventory", "{read_set}", "{i}.lengths.tsv"),
    params:
        source=lambda wc: read_set_file(wc.read_set, wc.i),
    threads: 4
    conda:
        os.path.join(ENVS, "reads.yaml")
    log:
        os.path.join(LOGS, "inventory", "{read_set}.{i}.lengths.log"),
    shell:
        # streamed sources are read without their index (some published .bai files do not match the BAM)
        """
        samtools quickcheck {params.source} 2> {log} || [ ! -f {params.source} ]
        samtools view -@ {threads} -F 0x900 {params.source} 2>> {log} \
          | awk -F'\\t' '{{n[length($10)]++}} END {{for (l in n) print l"\\t"n[l]}}' \
          | sort -k1,1n > {output}
        """


rule coverage_per_length:
    """Per read set: N50 and total coverage, and the coverage left when every read is cut into pieces of
    each length (the remainder of a read is lost)"""
    input:
        lambda wc: expand(os.path.join(OUT, "inventory", wc.read_set, "{i}.lengths.tsv"),
                          i=range(len(READ_SETS[wc.read_set]["urls"]))),
    output:
        os.path.join(OUT, "inventory", "{read_set}.coverage_per_length.tsv"),
    params:
        genome_size=config["inventory"]["genome_size"],
        lengths=config["inventory"]["lengths"],
    conda:
        os.path.join(ENVS, "python.yaml")
    log:
        os.path.join(LOGS, "inventory", "{read_set}.coverage_per_length.log"),
    shell:
        """
        python {SCRIPTS}/coverage_per_length.py --histograms {input} --genome-size {params.genome_size} \
          --lengths {params.lengths} --output {output} &> {log}
        """


def inventory_targets():
    return [os.path.join(OUT, "inventory", f"{read_set}.coverage_per_length.tsv") for read_set in READ_SETS]
