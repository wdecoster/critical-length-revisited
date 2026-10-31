# Assembly follow-up, regional arm (docs/assembly_plan.md, milestone 1). Windows are assembled per haplotype from
# the reads that land in them, as a regional pangenome pipeline does, at a grid of library N50s and coverages:
#   windows (checkpoint) -> selection regions (windows +/- margin, SD partners) -> full-length parent reads
#   (from a one-time position-only alignment of the read set) -> cut (N50 or fixed length) -> aligned to the whole
#   genome -> subsampled -> per window: the reads in window +/- flank
#   realistic: DeepVariant SNVs -> WhatsHap phase + haplotag -> H1/H2, untagged to both
#   oracle:    each read to the truth haplotype it aligns to best (ceiling)
#   diploid:   hifiasm on all reads of the window
#   -> Flye / hifiasm per haplotype -> evaluated against the window's HG002 Q100 haplotypes.

REG = config["regional"]
REG_DIR = os.path.join(OUT, "regional")
REG_RES = os.path.join(RES, "regional")
REG_RUN = os.path.join(REG_DIR, "L{length}", "s{seed}")
WIN = os.path.join(REG_RUN, "{coverage}x", "{window}")
Q100 = os.path.join(RES, "assemblies", "HG002_Q100")


wildcard_constraints:
    window=r"[A-Za-z0-9_]+",
    mode="realistic|oracle|diploid",
    assembler="flye|hifiasm",
    h="1|2",


def regional_lengths():
    return [f"N{n}{REG['shape']}" for n in REG["n50s"]] + [str(l) for l in REG["lengths"]]


def regional_windows():
    """Window names from the checkpoint's table"""
    with open(checkpoints.regional_windows.get().output.tsv) as fh:
        next(fh)
        return [line.split("\t")[0] for line in fh]


def regional_assemblies():
    """(mode, assembler) pairs"""
    pairs = [(m, a) for m in REG["modes"] for a in REG["assemblers"]]
    return pairs + ([("diploid", "hifiasm")] if REG["diploid"] else [])


def regional_evaluations(wildcards):
    return [
        os.path.join(REG_DIR, f"L{length}", f"s{seed}", f"{coverage}x", window, mode, f"{assembler}.hap{h}.eval.tsv")
        for length in regional_lengths()
        for seed in REG["seeds"]
        for coverage in REG["coverages"]
        for window in regional_windows()
        for mode, assembler in regional_assemblies()
        for h in (1, 2)
    ]


WINDOW_REGION = """$(awk -F '\\t' -v w={wildcards.window} -v f={params.flank} \
    '$1 == w {{s = $3 - f; if (s < 0) s = 0; print $2 ":" s "-" $4 + f}}' {input.windows})"""


rule reference_gaps:
    input:
        REF,
    output:
        os.path.join(REG_RES, "reference_gaps.bed"),
    conda:
        os.path.join(ENVS, "python.yaml")
    log:
        os.path.join(LOGS, "regional", "reference_gaps.log"),
    shell:
        "python {SCRIPTS}/reference_gaps.py {input} > {output} 2> {log}"


checkpoint regional_windows:
    input:
        fai=REF + ".fai",
        gaps=rules.reference_gaps.output,
        benchmark=BENCHMARK_BED,
        segdups=os.path.join(RES, "stratifications", "segdups.bed.gz"),
        superdups=rules.superdups.output,
        tandem=os.path.join(RES, "stratifications", "tandem_repeats.bed.gz"),
        truth=TRUTH_VCF,
    output:
        tsv=os.path.join(REG_RES, "windows.tsv"),
        bed=os.path.join(REG_RES, "windows.bed"),
    params:
        classes=" ".join(f"{c}={n}" for c, n in REG["windows"]["classes"].items()),
        named=" ".join(f"{n}={r}" for n, r in (REG["windows"].get("named") or {}).items()),
        size=REG["windows"]["size"],
        seed=REG["windows"]["seed"],
    conda:
        os.path.join(ENVS, "analysis.yaml")
    log:
        os.path.join(LOGS, "regional", "windows.log"),
    shell:
        """
        python {SCRIPTS}/select_windows.py --fai {input.fai} --gaps {input.gaps} --benchmark {input.benchmark} \
          --segdups {input.segdups} --superdups {input.superdups} --tandem-repeats {input.tandem} \
          --truth {input.truth} --size {params.size} --seed {params.seed} --classes {params.classes} \
          --named {params.named} --output {output.tsv} --bed {output.bed} &> {log}
        """


rule regional_regions:
    input:
        windows=rules.regional_windows.output.tsv,
        fai=REF + ".fai",
        superdups=rules.superdups.output,
        benchmark=BENCHMARK_BED,
    output:
        selection=os.path.join(REG_RES, "selection.bed"),
        evaluation=os.path.join(REG_RES, "evaluation.bed"),
        sizes=os.path.join(REG_RES, "sizes.tsv"),
    params:
        margin=REG["selection_margin"],
        flank=REG["flank"],
    conda:
        os.path.join(ENVS, "python.yaml")
    log:
        os.path.join(LOGS, "regional", "regions.log"),
    shell:
        """
        python {SCRIPTS}/regional_regions.py --windows {input.windows} --fai {input.fai} \
          --superdups {input.superdups} --benchmark {input.benchmark} --margin {params.margin} \
          --flank {params.flank} --selection {output.selection} --evaluation {output.evaluation} \
          --sizes {output.sizes} &> {log}
        """


rule read_positions:
    """Where every read of a read set maps (minimap2 lr:hq, positions only, no secondary): computed once, any
    window selection is then a lookup"""
    input:
        lambda wc: read_set_file(wc.read_set, wc.i),
    output:
        os.path.join(RES, "positions", "{read_set}.{i}.paf.gz"),
    threads: 16
    resources:
        mem_mb=24000,
    conda:
        os.path.join(ENVS, "minimap2.yaml")
    log:
        os.path.join(LOGS, "regional", "positions.{read_set}.{i}.log"),
    benchmark:
        os.path.join(LOGS, "regional", "positions.{read_set}.{i}.benchmark.tsv")
    shell:
        """
        samtools fastq -@ 4 {input} 2> {log} \
          | minimap2 -x lr:hq -t {threads} --secondary=no {REF} - 2>> {log} | gzip -1 > {output}
        """


rule regional_parents:
    """Full-length reads with a mapping in the selection regions, as FASTQ"""
    input:
        positions=lambda wc: expand(os.path.join(RES, "positions", "{read_set}.{i}.paf.gz"), read_set=REG["read_set"],
                                    i=range(len(READ_SETS[REG["read_set"]]["urls"]))),
        sources=lambda wc: [read_set_file(REG["read_set"], i) for i in range(len(READ_SETS[REG["read_set"]]["urls"]))],
        selection=rules.regional_regions.output.selection,
    output:
        fastq=os.path.join(REG_DIR, "parents.fastq.gz"),
        names=temp(os.path.join(REG_DIR, "parents.names.txt")),
    threads: 8
    conda:
        os.path.join(ENVS, "minimap2.yaml")
    log:
        os.path.join(LOGS, "regional", "parents.log"),
    shell:
        """
        zcat {input.positions} | python {SCRIPTS}/paf_overlap.py --regions {input.selection} 2> {log} \
          | sort -u -S 4G > {output.names}
        for source in {input.sources}; do
            samtools view -@ 4 -u -F 0x900 -N {output.names} $source 2>> {log} | samtools fastq -@ 4 - 2>> {log}
        done | gzip -1 > {output.fastq}
        """


def regional_cut_args(wc, input):
    if is_distribution(wc.length):
        return f"--distribution {' '.join(input.shape)} --n50 {nominal(wc.length)}"
    return f"--length {wc.length}"


rule regional_cut:
    """Parent reads cut into pieces (a library N50 of the configured shape, or a fixed length), kept at random to
    about oversample x the highest coverage x the selection regions"""
    input:
        fastq=rules.regional_parents.output.fastq,
        sizes=rules.regional_regions.output.sizes,
        shape=lambda wc: expand(os.path.join(OUT, "inventory", REG["shape"], "{i}.lengths.tsv"),
                                i=range(len(READ_SETS[REG["shape"]]["urls"]))) if is_distribution(wc.length) else [],
    output:
        fastq=temp(os.path.join(REG_RUN, "pieces.fastq.gz")),
        stats=os.path.join(REG_RUN, "pieces.stats.tsv"),
    params:
        args=regional_cut_args,
        target=lambda wc, input: REG["oversample"] * max(REG["coverages"]) * pilot_size(input.sizes, "selection"),
    conda:
        os.path.join(ENVS, "analysis.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.cut.log"),
    shell:
        """
        python {SCRIPTS}/cut_reads.py {input.fastq} {params.args} --seed {wildcards.seed} \
          --target-bases {params.target} --output {output.fastq} --stats {output.stats} &> {log}
        """


rule regional_align:
    """Pieces aligned to the whole genome (minimap2 lr:hq), and their mean coverage in the evaluation regions"""
    input:
        fastq=rules.regional_cut.output.fastq,
        evaluation=rules.regional_regions.output.evaluation,
    output:
        bam=temp(os.path.join(REG_RUN, "aligned.bam")),
        bai=temp(os.path.join(REG_RUN, "aligned.bam.bai")),
        coverage=os.path.join(REG_RUN, "aligned.coverage.txt"),
    params:
        prefix=lambda wc, output: output.coverage.removesuffix(".coverage.txt") + ".mosdepth",
    threads: 16
    conda:
        os.path.join(ENVS, "minimap2.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.align.log"),
    benchmark:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.align.benchmark.tsv")
    shell:
        """
        minimap2 -ax lr:hq -t {threads} {REF} {input.fastq} 2> {log} \
          | samtools sort -@ 4 -m 2G -o {output.bam} - 2>> {log}
        samtools index {output.bam}
        mosdepth -n --fast-mode -t 4 --by {input.evaluation} {params.prefix} {output.bam} 2>> {log}
        awk '$1 == "total_region" {{print $4}}' {params.prefix}.mosdepth.summary.txt > {output.coverage}
        rm -f {params.prefix}.*
        """


rule regional_subsample:
    """Exact coverage in the evaluation regions (all reads if it falls short: the achieved coverage is recorded)"""
    input:
        bam=rules.regional_align.output.bam,
        bai=rules.regional_align.output.bai,
        coverage=rules.regional_align.output.coverage,
    output:
        bam=os.path.join(REG_RUN, "{coverage}x.bam"),
        bai=os.path.join(REG_RUN, "{coverage}x.bam.bai"),
        coverage=os.path.join(REG_RUN, "{coverage}x.coverage.txt"),
    threads: 4
    conda:
        os.path.join(ENVS, "minimap2.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.subsample.log"),
    shell:
        """
        full=$(cat {input.coverage})
        fraction=$(awk -v c={wildcards.coverage} -v f=$full 'BEGIN {{x = c / f; print (x > 1 ? 1 : x)}}')
        if [ "$fraction" = "1" ]; then
            echo "only $full x available for {wildcards.coverage}x: all reads used" > {log}
            samtools view -b -o {output.bam} {input.bam} 2>> {log}
        else
            samtools view -@ {threads} -b --subsample $fraction --subsample-seed {wildcards.seed} \
              -o {output.bam} {input.bam} 2> {log}
        fi
        samtools index {output.bam}
        awk -v c={wildcards.coverage} -v f=$full 'BEGIN {{print (f < c ? f : c)}}' > {output.coverage}
        """


rule window_reads:
    """The alignments in the window +/- flank, and the reads touching it (primary records, full sequence)"""
    input:
        bam=rules.regional_subsample.output.bam,
        bai=rules.regional_subsample.output.bai,
        windows=rules.regional_windows.output.tsv,
    output:
        bam=os.path.join(WIN, "reads.bam"),
        bai=os.path.join(WIN, "reads.bam.bai"),
        fastq=os.path.join(WIN, "reads.fastq.gz"),
        names=temp(os.path.join(WIN, "reads.names.txt")),
    params:
        flank=REG["flank"],
    threads: 2
    conda:
        os.path.join(ENVS, "minimap2.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.reads.log"),
    shell:
        """
        region=""" + WINDOW_REGION + """
        samtools view -b -o {output.bam} {input.bam} $region 2> {log}
        samtools index {output.bam}
        samtools view {output.bam} | cut -f 1 | sort -u > {output.names}
        samtools view -u -F 0x900 -N {output.names} {input.bam} 2>> {log} | samtools fastq - 2>> {log} \
          | gzip -1 > {output.fastq}
        """


rule window_deepvariant:
    input:
        bam=rules.window_reads.output.bam,
        bai=rules.window_reads.output.bai,
        windows=rules.regional_windows.output.tsv,
    output:
        vcf=os.path.join(WIN, "realistic", "deepvariant.vcf.gz"),
    params:
        flank=REG["flank"],
        model=REG["deepvariant_model"],
        tmp=lambda wc, output: os.path.join(os.path.dirname(output.vcf), "deepvariant_tmp"),
    threads: 4
    resources:
        mem_mb=16000,
    container:
        REG["deepvariant_container"]
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.deepvariant.log"),
    benchmark:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.deepvariant.benchmark.tsv")
    shell:
        """
        region=""" + WINDOW_REGION + """
        run_deepvariant --model_type={params.model} --ref={REF} --reads={input.bam} --regions=$region \
          --sample_name=HG002 --output_vcf={output.vcf} --intermediate_results_dir={params.tmp} \
          --num_shards={threads} --call_variants_extra_args='allow_empty_examples=true' &> {log}
        rm -rf {params.tmp}
        """


rule window_phase:
    """WhatsHap: PASS SNVs phased with the window's reads, then the reads haplotagged (H1/H2/none per read)"""
    input:
        vcf=rules.window_deepvariant.output.vcf,
        bam=rules.window_reads.output.bam,
        bai=rules.window_reads.output.bai,
    output:
        vcf=os.path.join(WIN, "realistic", "phased.vcf.gz"),
        tags=os.path.join(WIN, "realistic", "haplotag.tsv.gz"),
        tagged=temp(os.path.join(WIN, "realistic", "tagged.bam")),
    threads: 1
    conda:
        os.path.join(ENVS, "phasing.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.whatshap.log"),
    shell:
        """
        bcftools view -f PASS -v snps -Oz -o {output.vcf}.input.vcf.gz {input.vcf} 2> {log}
        tabix -f -p vcf {output.vcf}.input.vcf.gz
        whatshap phase --ignore-read-groups --sample HG002 --only-snvs -r {REF} -o {output.vcf} \
          {output.vcf}.input.vcf.gz {input.bam} &>> {log}
        tabix -f -p vcf {output.vcf}
        whatshap haplotag --ignore-read-groups --sample HG002 -r {REF} --output-haplotag-list {output.tags} \
          -o {output.tagged} {output.vcf} {input.bam} &>> {log}
        rm -f {output.vcf}.input.vcf.gz {output.vcf}.input.vcf.gz.tbi
        """


rule realistic_reads:
    input:
        reads=rules.window_reads.output.fastq,
        tags=rules.window_phase.output.tags,
    output:
        hap1=os.path.join(WIN, "realistic", "hap1.fastq.gz"),
        hap2=os.path.join(WIN, "realistic", "hap2.fastq.gz"),
    params:
        minlen=REG["min_read_length"],
    conda:
        os.path.join(ENVS, "phasing.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.realistic_reads.log"),
    shell:
        """
        python {SCRIPTS}/haplotype_reads.py --reads {input.reads} --haplotag-list {input.tags} \
          --min-length {params.minlen} --output1 {output.hap1} --output2 {output.hap2} &> {log}
        """


rule window_truth:
    """HG002 Q100 haplotype sequence of every window (core) and window +/- flank (extended)"""
    input:
        windows=rules.regional_windows.output.tsv,
        hap1=os.path.join(Q100, "hap1.fa"),
        hap2=os.path.join(Q100, "hap2.fa"),
        fai=REF + ".fai",
    output:
        truth=directory(os.path.join(REG_RES, "truth")),
        status=os.path.join(REG_RES, "truth.status.tsv"),
    params:
        flank=REG["flank"],
        anchor=REG["anchor_length"],
        steps=" ".join(map(str, REG["anchor_steps"])),
    threads: 8
    resources:
        mem_mb=24000,
    conda:
        os.path.join(ENVS, "evaluation.yaml")
    log:
        os.path.join(LOGS, "regional", "truth.log"),
    shell:
        """
        python {SCRIPTS}/window_truth.py --windows {input.windows} --reference {REF} --hap1 {input.hap1} \
          --hap2 {input.hap2} --flank {params.flank} --anchor-length {params.anchor} --steps {params.steps} \
          --threads {threads} --outdir {output.truth} --status {output.status} &> {log}
        for f in {output.truth}/*.core.hap*.fa; do
            [ -s $f ] && yak count -k31 -t {threads} -o ${{f%.fa}}.yak $f 2>> {log}
        done
        """


rule oracle_reads:
    input:
        reads=rules.window_reads.output.fastq,
        truth=rules.window_truth.output.truth,
    output:
        hap1=os.path.join(WIN, "oracle", "hap1.fastq.gz"),
        hap2=os.path.join(WIN, "oracle", "hap2.fastq.gz"),
        assignments=os.path.join(WIN, "oracle", "assignments.tsv"),
    params:
        truth=lambda wc, input: os.path.join(input.truth, f"{wc.window}.extended.hap"),
        minlen=REG["min_read_length"],
    threads: 2
    conda:
        os.path.join(ENVS, "evaluation.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.oracle_reads.log"),
    shell:
        """
        python {SCRIPTS}/oracle_haplotypes.py --reads {input.reads} --hap1 {params.truth}1.fa \
          --hap2 {params.truth}2.fa --min-length {params.minlen} --threads {threads} --output1 {output.hap1} \
          --output2 {output.hap2} --assignments {output.assignments} &> {log}
        """


rule assemble_haploid:
    """One haplotype's reads assembled on their own (Flye --nano-hq, or hifiasm --ont haploid). A failed
    assembly (too few reads) gives an empty FASTA, reported as such."""
    input:
        reads=os.path.join(WIN, "{mode}", "hap{h}.fastq.gz"),
        windows=rules.regional_windows.output.tsv,
    output:
        os.path.join(WIN, "{mode}", "{assembler}.hap{h}.fa"),
    wildcard_constraints:
        mode="realistic|oracle",
    params:
        flank=REG["flank"],
        flye=REG["flye_extra"],
        hifiasm=REG["hifiasm_extra"],
        workdir=lambda wc, output: output[0].removesuffix(".fa") + "_work",
    threads: 4
    resources:
        mem_mb=16000,
    conda:
        os.path.join(ENVS, "assembly.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.{mode}.{assembler}.hap{h}.log"),
    benchmark:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.{mode}.{assembler}.hap{h}.benchmark.tsv")
    shell:
        """
        size=$(awk -F '\\t' -v w={wildcards.window} -v f={params.flank} '$1 == w {{print $4 - $3 + 2 * f}}' {input.windows})
        rm -rf {params.workdir} && mkdir -p {params.workdir}
        if [ "{wildcards.assembler}" = "flye" ]; then
            flye --nano-hq {input.reads} --out-dir {params.workdir} --genome-size $size --threads {threads} \
              {params.flye} &> {log} || true
            cp {params.workdir}/assembly.fasta {output} 2>> {log} || touch {output}
        else
            hifiasm --ont -o {params.workdir}/asm -t {threads} {params.hifiasm} {input.reads} &> {log} || true
            awk '/^S/ {{print ">"$2"\\n"$3}}' {params.workdir}/asm.bp.p_ctg.gfa > {output} 2>> {log} || touch {output}
        fi
        rm -rf {params.workdir}
        """


rule assemble_diploid:
    """All reads of the window assembled by hifiasm --ont in its diploid mode (its own haplotype separation)"""
    input:
        reads=rules.window_reads.output.fastq,
    output:
        hap1=os.path.join(WIN, "diploid", "hifiasm.hap1.fa"),
        hap2=os.path.join(WIN, "diploid", "hifiasm.hap2.fa"),
    params:
        workdir=lambda wc, output: os.path.join(os.path.dirname(output.hap1), "hifiasm_work"),
    threads: 4
    resources:
        mem_mb=16000,
    conda:
        os.path.join(ENVS, "assembly.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.diploid.hifiasm.log"),
    benchmark:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.diploid.hifiasm.benchmark.tsv")
    shell:
        """
        rm -rf {params.workdir} && mkdir -p {params.workdir}
        hifiasm --ont -o {params.workdir}/asm -t {threads} {input.reads} &> {log} || true
        for h in 1 2; do
            out={output.hap1}; [ $h = 2 ] && out={output.hap2}
            awk '/^S/ {{print ">"$2"\\n"$3}}' {params.workdir}/asm.bp.hap$h.p_ctg.gfa > $out 2>> {log} || touch $out
        done
        rm -rf {params.workdir}
        """


rule evaluate_haplotype:
    """One assembled haplotype against both truth haplotypes of the window (assembly_eval.py). yak trioeval
    with -c 1 -d 1: the truth k-mers come from assembled sequence (each seen once), not from reads"""
    input:
        contigs=os.path.join(WIN, "{mode}", "{assembler}.hap{h}.fa"),
        truth=rules.window_truth.output.truth,
        coverage=rules.regional_subsample.output.coverage,
    output:
        paf1=os.path.join(WIN, "{mode}", "{assembler}.hap{h}.truth1.paf"),
        paf2=os.path.join(WIN, "{mode}", "{assembler}.hap{h}.truth2.paf"),
        yak=os.path.join(WIN, "{mode}", "{assembler}.hap{h}.yak.txt"),
        tsv=os.path.join(WIN, "{mode}", "{assembler}.hap{h}.eval.tsv"),
    params:
        truth=lambda wc, input: os.path.join(input.truth, f"{wc.window}.core.hap"),
        complete=REG["complete_fraction"],
    threads: 2
    conda:
        os.path.join(ENVS, "evaluation.yaml")
    log:
        os.path.join(LOGS, "regional", "L{length}.s{seed}.{coverage}x.{window}.{mode}.{assembler}.hap{h}.eval.log"),
    shell:
        """
        for t in 1 2; do
            paf={output.paf1}; [ $t = 2 ] && paf={output.paf2}
            if [ -s {input.contigs} ] && [ -s {params.truth}$t.fa ]; then
                minimap2 -cx asm5 -t {threads} {params.truth}$t.fa {input.contigs} > $paf 2>> {log}
            else
                touch $paf
            fi
        done
        if [ -s {input.contigs} ] && [ -s {params.truth}1.yak ] && [ -s {params.truth}2.yak ]; then
            yak trioeval -c 1 -d 1 -t {threads} {params.truth}1.yak {params.truth}2.yak {input.contigs} > {output.yak} 2>> {log}
        else
            touch {output.yak}
        fi
        python {SCRIPTS}/assembly_eval.py --contigs {input.contigs} --truth1 {params.truth}1.fa \
          --truth2 {params.truth}2.fa --paf1 {output.paf1} --paf2 {output.paf2} --yak {output.yak} \
          --complete-fraction {params.complete} --output {output.tsv} \
          --labels length={wildcards.length} seed={wildcards.seed} coverage={wildcards.coverage} \
          achieved_coverage=$(cat {input.coverage}) window={wildcards.window} mode={wildcards.mode} \
          assembler={wildcards.assembler} haplotype={wildcards.h} &>> {log}
        """


rule regional_summary:
    input:
        evaluations=regional_evaluations,
        windows=rules.regional_windows.output.tsv,
        truth=rules.window_truth.output.status,
    output:
        os.path.join(REG_DIR, "summary.tsv"),
    conda:
        os.path.join(ENVS, "analysis.yaml")
    log:
        os.path.join(LOGS, "regional", "summary.log"),
    shell:
        """
        python -c "import sys, pandas as pd; w = pd.read_table(sys.argv[2]).rename(columns={{'name': 'window'}}); \
pd.concat(pd.read_table(f) for f in sys.argv[3:]).merge(w, on='window', how='left').to_csv(sys.argv[1], sep='\\t', index=False)" \
          {output} {input.windows} {input.evaluations} &> {log}
        """
