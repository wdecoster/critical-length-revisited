# Then and now, steps 1-2 (docs/plan.md; called the ladder in file and rule names): the 2020 pipeline, unchanged, on HG00733 (the 2020 assembly) and on
# HG002-Q100 v1.1. Reads: SimLoRD 1.0.2 fixed-length reads with its default (multi-pass PacBio) error model,
# simulated from both haplotypes (--without-ns), aligned with minimap2 2.14 map-pb, SVs called with
# Sniffles 1.0.10. Truth: the assembly haplotypes aligned to GRCh38 (minimap2 2.14 asm5, as dipcall) and
# called with paftools, SVs >= 50 bp, haplotypes merged. Matching: SURVIVOR merge as surpyvor prf.

LADDER_DIR = os.path.join(OUT, "ladder")


def maybe_temp(path):
    return path if LADDER["keep_bams"] else temp(path)


def ladder_bam(wildcards):
    """The reads at a coverage: the simulated alignments, or a subsample of them"""
    name = "simulated" if int(wildcards.coverage) == LADDER["simulated_coverage"][wildcards.genome] else f"{wildcards.coverage}x"
    return os.path.join(LADDER_DIR, wildcards.genome, f"L{wildcards.length}", f"{name}.bam")


rule assembly_chunks:
    """Both haplotypes split into chunks of whole sequences of similar total length, simulated in parallel.
    Every chunk gets one of the longest sequences first, so each can provide reads of the longest length."""
    input:
        fai=expand(os.path.join(RES, "assemblies", "{{genome}}", "{hap}.fa.fai"), hap=["hap1", "hap2"]),
        fasta=expand(os.path.join(RES, "assemblies", "{{genome}}", "{hap}.fa"), hap=["hap1", "hap2"]),
    output:
        expand(os.path.join(RES, "assemblies", "{{genome}}", "chunks", "{chunk}.fa"), chunk=CHUNKS),
    params:
        outdir=lambda wc: os.path.join(RES, "assemblies", wc.genome, "chunks"),
        chunks=LADDER["chunks"],
        max_length=max(LADDER["lengths"]),
    conda:
        os.path.join(ENVS, "reads.yaml")
    log:
        os.path.join(LOGS, "ladder", "{genome}.chunks.log"),
    shell:
        """
        python {SCRIPTS}/assembly_chunks.py --fasta {input.fasta} --chunks {params.chunks} \
          --min-longest {params.max_length} --outdir {params.outdir} &> {log}
        """


rule simlord:
    """SimLoRD as in 2020 (-fl, --without-ns, default error model); -c is per chunk, so coverage is per
    haplotype: half the haploid-equivalent coverage"""
    input:
        os.path.join(RES, "assemblies", "{genome}", "chunks", "{chunk}.fa"),
    output:
        temp(os.path.join(LADDER_DIR, "{genome}", "L{length}", "reads", "{chunk}.fastq.gz")),
    params:
        coverage=lambda wc: LADDER["simulated_coverage"][wc.genome] // 2,
        prefix=lambda wc, output: output[0].removesuffix(".fastq.gz"),
    conda:
        os.path.join(ENVS, "simlord.yaml")
    log:
        os.path.join(LOGS, "ladder", "{genome}.L{length}.{chunk}.simlord.log"),
    shell:
        """
        simlord --read-reference {input} -c {params.coverage} --without-ns -fl {wildcards.length} --no-sam \
          {params.prefix} &> {log}
        gzip -1 {params.prefix}.fastq  # -> {output}
        """


rule ladder_align:
    """minimap2 2.14 -ax map-pb --MD, as in 2020"""
    input:
        expand(os.path.join(LADDER_DIR, "{{genome}}", "L{{length}}", "reads", "{chunk}.fastq.gz"), chunk=CHUNKS),
    output:
        bam=maybe_temp(os.path.join(LADDER_DIR, "{genome}", "L{length}", "simulated.bam")),
        bai=maybe_temp(os.path.join(LADDER_DIR, "{genome}", "L{length}", "simulated.bam.bai")),
    threads: 32
    priority: 10  # finish alignments (and delete their reads) before simulating more
    conda:
        os.path.join(ENVS, "minimap2_2014.yaml")
    log:
        os.path.join(LOGS, "ladder", "{genome}.L{length}.align.log"),
    shell:
        """
        minimap2 --MD -ax map-pb -t {threads} {REF} <(zcat {input}) 2> {log} \
          | samtools sort -@ 4 -m 2G -o {output.bam} - 2>> {log}
        samtools index {output.bam}
        """


rule ladder_subsample:
    """A coverage below the simulated one, subsampled from the simulated reads"""
    input:
        bam=os.path.join(LADDER_DIR, "{genome}", "L{length}", "simulated.bam"),
        bai=os.path.join(LADDER_DIR, "{genome}", "L{length}", "simulated.bam.bai"),
    output:
        bam=maybe_temp(os.path.join(LADDER_DIR, "{genome}", "L{length}", "{coverage}x.bam")),
        bai=maybe_temp(os.path.join(LADDER_DIR, "{genome}", "L{length}", "{coverage}x.bam.bai")),
    params:
        fraction=lambda wc: int(wc.coverage) / LADDER["simulated_coverage"][wc.genome],
    threads: 4
    conda:
        os.path.join(ENVS, "minimap2_2014.yaml")
    shell:
        """
        samtools view -@ {threads} -b --subsample {params.fraction} --subsample-seed 42 -o {output.bam} {input.bam}
        samtools index {output.bam}
        """


rule sniffles1:
    """Sniffles 1.0.10 as in 2020 (--genotype, defaults), or with another setting"""
    input:
        bam=ladder_bam,
        bai=lambda wc: ladder_bam(wc) + ".bai",
    output:
        os.path.join(LADDER_DIR, "{genome}", "L{length}", "{coverage}x.sniffles1_{setting}.vcf"),
    params:
        args=lambda wc: LADDER["sniffles1"][wc.setting]["args"],
    threads: 4
    conda:
        os.path.join(ENVS, "sniffles1.yaml")
    log:
        os.path.join(LOGS, "ladder", "{genome}.L{length}.{coverage}x.sniffles1_{setting}.log"),
    shell:
        "sniffles -m {input.bam} -v {output} -t {threads} --genotype {params.args} &> {log}"


rule truth_2020_haplotype:
    """One haplotype against GRCh38 as dipcall does it (minimap2 asm5, paftools call), SVs >= 50 bp as INS/DEL"""
    input:
        os.path.join(RES, "assemblies", "{genome}", "{hap}.fa"),
    output:
        paf=os.path.join(LADDER_DIR, "{genome}", "truth_2020", "{hap}.paf.gz"),
        vcf=os.path.join(LADDER_DIR, "{genome}", "truth_2020", "{hap}.svs.vcf"),
    params:
        min_length=LADDER["min_sv_length"],
    threads: 16
    conda:
        os.path.join(ENVS, "minimap2_2014.yaml")
    log:
        os.path.join(LOGS, "ladder", "{genome}.truth_2020.{hap}.log"),
    shell:
        """
        minimap2 -cx asm5 --cs -r2k -t {threads} {REF} {input} 2> {log} | gzip > {output.paf}
        zcat {output.paf} | sort -k6,6 -k8,8n | paftools.js call -f {REF} - 2>> {log} \
          | python {SCRIPTS}/sv_from_paftools.py --min-length {params.min_length} --sample {wildcards.hap} \
          > {output.vcf} 2>> {log}
        """


rule truth_2020:
    """Both haplotypes' SVs merged into one diploid truth (SURVIVOR merge, as surpyvor haplomerge)"""
    input:
        expand(os.path.join(LADDER_DIR, "{{genome}}", "truth_2020", "{hap}.svs.vcf"), hap=["hap1", "hap2"]),
    output:
        os.path.join(LADDER_DIR, "{genome}", "truth_2020", "truth.vcf"),
    params:
        distance=LADDER["match_distance"],
        min_length=LADDER["min_sv_length"],
        tmp=lambda wc: os.path.join(LADDER_DIR, wc.genome, "truth_2020", "tmp"),
    conda:
        os.path.join(ENVS, "survivor.yaml")
    log:
        os.path.join(LOGS, "ladder", "{genome}.truth_2020.log"),
    shell:
        """
        mkdir -p {params.tmp}
        ls {input} > {params.tmp}/files.txt
        SURVIVOR merge {params.tmp}/files.txt {params.distance} 1 1 -1 -1 {params.min_length} \
          {params.tmp}/merged.vcf &> {log}
        python {SCRIPTS}/combine_haplotypes.py --vcf {params.tmp}/merged.vcf --sample {wildcards.genome} \
          > {output} 2>> {log}
        rm -rf {params.tmp}
        """


rule evaluate_2020:
    """Precision, recall and F1 as surpyvor prf computed them: SURVIVOR merge of truth and calls (distance
    500, same type, min. length 50, after DUP -> INS), counted over merged sites"""
    input:
        truth=os.path.join(LADDER_DIR, "{genome}", "truth_2020", "truth.vcf"),
        calls=os.path.join(LADDER_DIR, "{genome}", "L{length}", "{coverage}x.sniffles1_{setting}.vcf"),
    output:
        os.path.join(LADDER_DIR, "{genome}", "L{length}", "{coverage}x.sniffles1_{setting}.eval_2020.tsv"),
    params:
        distance=LADDER["match_distance"],
        min_length=LADDER["min_sv_length"],
    conda:
        os.path.join(ENVS, "survivor.yaml")
    log:
        os.path.join(LOGS, "ladder", "{genome}.L{length}.{coverage}x.sniffles1_{setting}.eval_2020.log"),
    shell:
        """
        python {SCRIPTS}/evaluate_2020.py --truth {input.truth} --calls {input.calls} \
          --distance {params.distance} --min-length {params.min_length} \
          --labels genome={wildcards.genome} length={wildcards.length} coverage={wildcards.coverage} \
          setting={wildcards.setting} --output {output} &> {log}
        """


rule ladder_summary:
    input:
        [os.path.join(LADDER_DIR, g, f"L{l}", f"{c}x.sniffles1_{s}.eval_2020.tsv") for g, l, c, s in ladder_runs()],
    output:
        os.path.join(LADDER_DIR, "summary_steps1-2.tsv"),
    conda:
        os.path.join(ENVS, "python.yaml")
    shell:
        """
        python -c "import sys, pandas as pd; pd.concat(pd.read_table(f) for f in sys.argv[2:]).to_csv(sys.argv[1], sep='\\t', index=False)" \
          {output} {input}
        """
