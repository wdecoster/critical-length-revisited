# Stage 2, the pilot (docs/plan.md, milestone 3): Arm A on the pilot chromosomes. Real HG002 R10 reads are
# selected by where they map (pilot chromosomes plus SD paralog partners), cut into pieces of a fixed length
# (cut_reads.py), aligned to the whole genome (minimap2 lr:hq), subsampled to exact coverages in the
# evaluation regions, called (Sniffles2; cuteSV at a few lengths) and scored against GIAB v5.0q with truvari
# (strict: GIAB's setup + refine; lenient), stratified by SV type, size and context.

PILOT = config["pilot"]
PILOT_DIR = os.path.join(OUT, "pilot")
PILOT_RES = os.path.join(RES, "pilot")
TRUTH_VCF = os.path.join(RES, "giab", "HG002_GRCh38_v5.0q_stvar.vcf.gz")
BENCHMARK_BED = os.path.join(RES, "giab", "HG002_GRCh38_v5.0q_stvar.benchmark.bed")


def pilot_coverages(length):
    """Coverages a length is evaluated at, highest first"""
    return sorted((c for c, max_length in PILOT["coverages"].items() if int(length) <= max_length), reverse=True)


def pilot_seeds(read_set, length):
    seeds = list(PILOT["seeds"])
    rep = PILOT.get("replicates") or {}
    if read_set in rep.get("read_sets", []) and int(length) in rep.get("lengths", []):
        seeds += rep["seeds"]
    return seeds


def pilot_callers(length):
    return [c for c, s in PILOT["callers"].items() if s["lengths"] == "all" or int(length) in s["lengths"]]


def pilot_runs():
    """(read_set, length, seed, coverage, caller, scoring) of the pilot"""
    return [
        (read_set, length, seed, coverage, caller, scoring)
        for read_set, spec in PILOT["read_sets"].items()
        for length in spec["lengths"]
        for seed in pilot_seeds(read_set, length)
        for coverage in pilot_coverages(length)
        for caller in pilot_callers(length)
        for scoring in PILOT["scoring"]
    ]


def run_dir(wc):
    return os.path.join(PILOT_DIR, wc.read_set, f"L{wc.length}", f"s{wc.seed}")


RUN = os.path.join(PILOT_DIR, "{read_set}", "L{length}", "s{seed}")


wildcard_constraints:
    seed=r"\d+",
    caller="|".join(map(re.escape, PILOT["callers"])),
    scoring="|".join(map(re.escape, PILOT["scoring"])),


rule superdups:
    output:
        os.path.join(PILOT_RES, "genomicSuperDups.txt.gz"),
    params:
        url=PILOT["superdups_url"],
    log:
        os.path.join(LOGS, "pilot", "superdups.log"),
    shell:
        "curl -sSfL -o {output} {params.url} &> {log}"


rule pilot_regions:
    """Selection regions (pilot chromosomes + SD paralog partners elsewhere) and evaluation regions (GIAB
    benchmark regions on the pilot chromosomes), with their sizes"""
    input:
        fai=REF + ".fai",
        superdups=rules.superdups.output,
        benchmark=BENCHMARK_BED,
    output:
        selection=os.path.join(PILOT_RES, "selection.bed"),
        evaluation=os.path.join(PILOT_RES, "evaluation.bed"),
        sizes=os.path.join(PILOT_RES, "sizes.tsv"),
    params:
        chromosomes=" ".join(PILOT["chromosomes"]),
        partners="--sd-partners" if PILOT["sd_partners"] else "",
    conda:
        os.path.join(ENVS, "python.yaml")
    log:
        os.path.join(LOGS, "pilot", "regions.log"),
    shell:
        """
        python {SCRIPTS}/pilot_regions.py --fai {input.fai} --superdups {input.superdups} \
          --benchmark {input.benchmark} --chromosomes {params.chromosomes} {params.partners} \
          --selection {output.selection} --evaluation {output.evaluation} --sizes {output.sizes} &> {log}
        """


def source_is_aligned(read_set):
    return not READ_SETS[read_set]["download"]


rule parent_names:
    """Names of the reads of one source with an alignment (primary or supplementary) overlapping the selection
    regions. Aligned sources (streamed) use their own alignments; unaligned ones are mapped first (minimap2
    lr:hq without base-level alignment: positions only)"""
    input:
        regions=rules.pilot_regions.output.selection,
        source=lambda wc: [] if source_is_aligned(wc.read_set) else read_set_file(wc.read_set, wc.i),
    output:
        os.path.join(PILOT_DIR, "parents", "{read_set}.{i}.names.txt.gz"),
    params:
        source=lambda wc: read_set_file(wc.read_set, wc.i),
        aligned=lambda wc: "1" if source_is_aligned(wc.read_set) else "",
    threads: 32
    conda:
        os.path.join(ENVS, "minimap2.yaml")
    log:
        os.path.join(LOGS, "pilot", "parents.{read_set}.{i}.log"),
    shell:
        """
        if [ -n "{params.aligned}" ]; then
            # streamed through stdin, so -L filters the stream instead of using the (mismatched) remote index
            samtools view -@ 4 -u -F 0x104 {params.source} 2> {log} \
              | samtools view -L {input.regions} - 2>> {log} | cut -f 1
        else
            samtools fastq -@ 4 {params.source} 2> {log} \
              | minimap2 -x lr:hq -t {threads} --secondary=no {REF} - 2>> {log} \
              | python {SCRIPTS}/paf_overlap.py --regions {input.regions} 2>> {log}
        fi | sort -u -S 4G | gzip > {output}
        """


rule parent_reads:
    """The selected parent reads of a read set as FASTQ, in sequencing orientation (primary records)"""
    input:
        names=lambda wc: expand(os.path.join(PILOT_DIR, "parents", "{{read_set}}.{i}.names.txt.gz"),
                                i=range(len(READ_SETS[wc.read_set]["urls"]))),
        sources=lambda wc: [] if source_is_aligned(wc.read_set) else
            [read_set_file(wc.read_set, i) for i in range(len(READ_SETS[wc.read_set]["urls"]))],
    output:
        fastq=os.path.join(PILOT_DIR, "parents", "{read_set}.fastq.gz"),
        names=temp(os.path.join(PILOT_DIR, "parents", "{read_set}.names.txt")),
    params:
        sources=lambda wc: " ".join(read_set_file(wc.read_set, i) for i in range(len(READ_SETS[wc.read_set]["urls"]))),
    threads: 8
    conda:
        os.path.join(ENVS, "minimap2.yaml")
    log:
        os.path.join(LOGS, "pilot", "parents.{read_set}.fastq.log"),
    shell:
        """
        zcat {input.names} | sort -u > {output.names}
        for source in {params.sources}; do
            samtools view -@ 4 -u -F 0x900 -N {output.names} $source 2>> {log} | samtools fastq -@ 4 - 2>> {log}
        done | gzip -1 > {output.fastq}
        """


rule cut_reads:
    """Parent reads cut into pieces of one length (one placement seed), kept at random to about oversample x
    the highest coverage x the selection regions"""
    input:
        fastq=os.path.join(PILOT_DIR, "parents", "{read_set}.fastq.gz"),
        sizes=rules.pilot_regions.output.sizes,
    output:
        fastq=temp(os.path.join(RUN, "pieces.fastq.gz")),
        stats=os.path.join(RUN, "pieces.stats.tsv"),
    params:
        target=lambda wc, input: PILOT["oversample"] * max(pilot_coverages(wc.length)) * pilot_size(input.sizes, "selection"),
    conda:
        os.path.join(ENVS, "analysis.yaml")
    log:
        os.path.join(LOGS, "pilot", "{read_set}.L{length}.s{seed}.cut.log"),
    shell:
        """
        python {SCRIPTS}/cut_reads.py {input.fastq} --length {wildcards.length} --seed {wildcards.seed} \
          --target-bases {params.target} --output {output.fastq} --stats {output.stats} &> {log}
        """


def pilot_size(sizes, name):
    return float(dict(line.split() for line in open(sizes))[name])


rule pilot_align:
    """Pieces aligned to the whole genome (minimap2 lr:hq), and their mean coverage in the evaluation regions"""
    input:
        fastq=os.path.join(RUN, "pieces.fastq.gz"),
        evaluation=rules.pilot_regions.output.evaluation,
    output:
        bam=temp(os.path.join(RUN, "aligned.bam")),
        bai=temp(os.path.join(RUN, "aligned.bam.bai")),
        coverage=os.path.join(RUN, "aligned.coverage.txt"),
    params:
        prefix=lambda wc, output: output.coverage.removesuffix(".coverage.txt") + ".mosdepth",
    threads: 32
    priority: 10  # align (and free the pieces) before cutting more
    conda:
        os.path.join(ENVS, "minimap2.yaml")
    log:
        os.path.join(LOGS, "pilot", "{read_set}.L{length}.s{seed}.align.log"),
    shell:
        """
        minimap2 -ax lr:hq -t {threads} {REF} {input.fastq} 2> {log} \
          | samtools sort -@ 4 -m 2G -o {output.bam} - 2>> {log}
        samtools index {output.bam}
        mosdepth -n --fast-mode -t 4 --by {input.evaluation} {params.prefix} {output.bam} 2>> {log}
        awk '$1 == "total_region" {{print $4}}' {params.prefix}.mosdepth.summary.txt > {output.coverage}
        rm -f {params.prefix}.*
        """


rule pilot_subsample:
    """Exact coverage in the evaluation regions, subsampled from the alignment (all of it if it falls short:
    the achieved coverage is recorded)"""
    input:
        bam=os.path.join(RUN, "aligned.bam"),
        bai=os.path.join(RUN, "aligned.bam.bai"),
        coverage=os.path.join(RUN, "aligned.coverage.txt"),
    output:
        bam=temp(os.path.join(RUN, "{coverage}x.bam")),
        bai=temp(os.path.join(RUN, "{coverage}x.bam.bai")),
        coverage=os.path.join(RUN, "{coverage}x.coverage.txt"),
    threads: 4
    conda:
        os.path.join(ENVS, "minimap2.yaml")
    log:
        os.path.join(LOGS, "pilot", "{read_set}.L{length}.s{seed}.{coverage}x.subsample.log"),
    shell:
        """
        full=$(cat {input.coverage})
        fraction=$(awk -v c={wildcards.coverage} -v f=$full 'BEGIN {{x = c / f; print (x > 1 ? 1 : x)}}')
        if [ "$fraction" = "1" ]; then
            echo "only $full x available for {wildcards.coverage}x: all reads used" > {log}
            ln -f {input.bam} {output.bam}
        else
            samtools view -@ {threads} -b --subsample $fraction --subsample-seed {wildcards.seed} \
              -o {output.bam} {input.bam} 2> {log}
        fi
        samtools index {output.bam}
        awk -v c={wildcards.coverage} -v f=$full 'BEGIN {{print (f < c ? f : c)}}' > {output.coverage}
        """


rule sniffles2:
    input:
        bam=os.path.join(RUN, "{coverage}x.bam"),
        bai=os.path.join(RUN, "{coverage}x.bam.bai"),
    output:
        vcf=os.path.join(RUN, "{coverage}x.sniffles2.vcf.gz"),
        tbi=os.path.join(RUN, "{coverage}x.sniffles2.vcf.gz.tbi"),
    params:
        args=PILOT["callers"]["sniffles2"]["args"],
    threads: 8
    conda:
        os.path.join(ENVS, "sniffles2.yaml")
    log:
        os.path.join(LOGS, "pilot", "{read_set}.L{length}.s{seed}.{coverage}x.sniffles2.log"),
    shell:
        """
        sniffles --input {input.bam} --vcf {output.vcf} --reference {REF} --threads {threads} \
          --sample-id HG002 --allow-overwrite {params.args} &> {log}
        """


rule cutesv:
    input:
        bam=os.path.join(RUN, "{coverage}x.bam"),
        bai=os.path.join(RUN, "{coverage}x.bam.bai"),
    output:
        vcf=os.path.join(RUN, "{coverage}x.cutesv.vcf.gz"),
        tbi=os.path.join(RUN, "{coverage}x.cutesv.vcf.gz.tbi"),
    params:
        args=PILOT["callers"]["cutesv"]["args"],
        support=lambda wc: max(PILOT["callers"]["cutesv"]["min_support"],
                               round(int(wc.coverage) * PILOT["callers"]["cutesv"]["support_fraction"])),
        workdir=lambda wc: os.path.join(run_dir(wc), f"{wc.coverage}x.cutesv_tmp"),
        vcf=lambda wc, output: output.vcf.removesuffix(".gz"),
    threads: 8
    conda:
        os.path.join(ENVS, "cutesv.yaml")
    log:
        os.path.join(LOGS, "pilot", "{read_set}.L{length}.s{seed}.{coverage}x.cutesv.log"),
    shell:
        """
        mkdir -p {params.workdir}
        cuteSV {input.bam} {REF} {params.vcf} {params.workdir} --threads {threads} --sample HG002 --genotype \
          -s {params.support} -l 50 -L -1 {params.args} &> {log}
        bcftools sort -Oz -o {output.vcf} {params.vcf} 2>> {log}
        tabix -f -p vcf {output.vcf}
        rm -rf {params.workdir} {params.vcf}
        """


rule pilot_truth:
    """GIAB v5.0q without ALT='*' records (GIAB's recommendation for truvari)"""
    input:
        TRUTH_VCF,
    output:
        vcf=os.path.join(PILOT_RES, "truth.vcf.gz"),
        tbi=os.path.join(PILOT_RES, "truth.vcf.gz.tbi"),
    conda:
        os.path.join(ENVS, "truvari.yaml")
    log:
        os.path.join(LOGS, "pilot", "truth.log"),
    shell:
        """
        bcftools view -e 'ALT="*"' -Oz -o {output.vcf} {input} 2> {log}
        tabix -f -p vcf {output.vcf}
        """


rule truvari:
    """truvari bench in the evaluation regions (strict: + refine), consolidated with truvari ga4gh into truth
    and query VCFs with per-record TP/FN/FP"""
    input:
        truth=rules.pilot_truth.output.vcf,
        calls=os.path.join(RUN, "{coverage}x.{caller}.vcf.gz"),
        evaluation=rules.pilot_regions.output.evaluation,
    output:
        truth=os.path.join(RUN, "{coverage}x.{caller}.{scoring}.truth.vcf.gz"),
        query=os.path.join(RUN, "{coverage}x.{caller}.{scoring}.query.vcf.gz"),
        summary=os.path.join(RUN, "{coverage}x.{caller}.{scoring}.summary.json"),
    params:
        args=lambda wc: PILOT["scoring"][wc.scoring]["args"],
        refine=lambda wc: "1" if PILOT["scoring"][wc.scoring]["refine"] else "",
        bench=lambda wc: os.path.join(run_dir(wc), f"{wc.coverage}x.{wc.caller}.{wc.scoring}.truvari"),
    threads: 8
    conda:
        os.path.join(ENVS, "truvari.yaml")
    log:
        os.path.join(LOGS, "pilot", "{read_set}.L{length}.s{seed}.{coverage}x.{caller}.{scoring}.truvari.log"),
    shell:
        """
        rm -rf {params.bench}
        truvari bench -b {input.truth} -c {input.calls} -f {REF} --includebed {input.evaluation} \
          -o {params.bench} {params.args} &> {log}
        if [ -n "{params.refine}" ]; then
            truvari refine -f {REF} --use-original-vcfs --align mafft -t {threads} {params.bench} &>> {log}
            truvari ga4gh -i {params.bench} -o {params.bench}/ga4gh &>> {log}
        else
            truvari ga4gh --no-refine -i {params.bench} -o {params.bench}/ga4gh &>> {log}
        fi
        # ga4gh: base and comp records with FORMAT/BD = TP/FN/FP (refined where refine ran)
        cp {params.bench}/ga4gh.base.vcf.gz {output.truth}
        cp {params.bench}/ga4gh.comp.vcf.gz {output.query}
        cp {params.bench}/ga4gh.summary.json {output.summary}
        """


rule stratify:
    """Recall, precision and F1 per stratum (type, size, context and their combinations)"""
    input:
        truth=os.path.join(RUN, "{coverage}x.{caller}.{scoring}.truth.vcf.gz"),
        query=os.path.join(RUN, "{coverage}x.{caller}.{scoring}.query.vcf.gz"),
        coverage=os.path.join(RUN, "{coverage}x.coverage.txt"),
        strata=lambda wc: [os.path.join(RES, "stratifications", f"{s}.bed.gz") for s in PILOT["contexts"]],
    output:
        os.path.join(RUN, "{coverage}x.{caller}.{scoring}.strata.tsv"),
    params:
        contexts=" ".join(f"{s}={os.path.join(RES, 'stratifications', s + '.bed.gz')}" for s in PILOT["contexts"]),
        bins=" ".join(map(str, PILOT["size_bins"])),
        chromosomes=" ".join(PILOT["chromosomes"]),
    conda:
        os.path.join(ENVS, "analysis.yaml")
    log:
        os.path.join(LOGS, "pilot", "{read_set}.L{length}.s{seed}.{coverage}x.{caller}.{scoring}.strata.log"),
    shell:
        """
        python {SCRIPTS}/stratify_truvari.py --truth {input.truth} --query {input.query} \
          --contexts {params.contexts} --size-bins {params.bins} --chromosomes {params.chromosomes} \
          --labels read_set={wildcards.read_set} length={wildcards.length} seed={wildcards.seed} \
          coverage={wildcards.coverage} achieved_coverage=$(cat {input.coverage}) caller={wildcards.caller} \
          scoring={wildcards.scoring} --output {output} &> {log}
        """


rule spanning_reads:
    """Per truth SV in the evaluation regions: reads spanning it with the flank on both sides"""
    input:
        bam=os.path.join(RUN, "{coverage}x.bam"),
        bai=os.path.join(RUN, "{coverage}x.bam.bai"),
        truth=rules.pilot_truth.output.vcf,
        evaluation=rules.pilot_regions.output.evaluation,
    output:
        os.path.join(RUN, "{coverage}x.spanning.tsv.gz"),
    params:
        flank=PILOT["spanning_flank"],
        mapq=PILOT["spanning_mapq"],
    conda:
        os.path.join(ENVS, "analysis.yaml")
    log:
        os.path.join(LOGS, "pilot", "{read_set}.L{length}.s{seed}.{coverage}x.spanning.log"),
    shell:
        """
        python {SCRIPTS}/spanning_reads.py --bam {input.bam} --truth {input.truth} --regions {input.evaluation} \
          --flank {params.flank} --min-mapq {params.mapq} --output {output} &> {log}
        """


rule keep_bam:
    """The alignments of keep_bams_lengths outlive their temporary originals (a hard link)"""
    input:
        bam=os.path.join(RUN, "{coverage}x.bam"),
        bai=os.path.join(RUN, "{coverage}x.bam.bai"),
    output:
        bam=os.path.join(RUN, "{coverage}x.kept.bam"),
        bai=os.path.join(RUN, "{coverage}x.kept.bam.bai"),
    shell:
        "ln -f {input.bam} {output.bam} && ln -f {input.bai} {output.bai}"


rule pilot_summary:
    """All strata of all pilot runs in one table"""
    input:
        strata=[os.path.join(PILOT_DIR, r, f"L{l}", f"s{s}", f"{c}x.{caller}.{scoring}.strata.tsv")
                for r, l, s, c, caller, scoring in pilot_runs()],
        spanning=sorted({os.path.join(PILOT_DIR, r, f"L{l}", f"s{s}", f"{c}x.spanning.tsv.gz")
                         for r, l, s, c, _, _ in pilot_runs()}),
        kept=sorted({os.path.join(PILOT_DIR, r, f"L{l}", f"s{s}", f"{c}x.kept.bam")
                     for r, l, s, c, _, _ in pilot_runs() if l in PILOT["keep_bams_lengths"]}),
    output:
        os.path.join(PILOT_DIR, "summary.tsv"),
    conda:
        os.path.join(ENVS, "analysis.yaml")
    shell:
        """
        python -c "import sys, pandas as pd; pd.concat(pd.read_table(f) for f in sys.argv[2:]).to_csv(sys.argv[1], sep='\\t', index=False)" \
          {output} {input.strata}
        """
