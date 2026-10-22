rule giab_file:
    """A GIAB v5.0q benchmark file, checked against its published md5"""
    output:
        os.path.join(RES, "giab", "{path}"),
    params:
        url=lambda wc: f"{config['truth']['giab_url']}/{wc.path}",
        md5=lambda wc: TRUTH_FILES[wc.path],
    wildcard_constraints:
        path="|".join(map(re.escape, TRUTH_FILES)),
    log:
        os.path.join(LOGS, "resources", "giab.{path}.log"),
    shell:
        """
        curl -sSfL -o {output} {params.url} &> {log}
        echo "{params.md5}  {output}" | md5sum -c - >> {log}
        """


rule giab_index:
    input:
        os.path.join(RES, "giab", "{path}.vcf.gz"),
    output:
        os.path.join(RES, "giab", "{path}.vcf.gz.tbi"),
    conda:
        os.path.join(ENVS, "bcftools.yaml")
    shell:
        "tabix -f -p vcf {input}"


rule stratification:
    output:
        os.path.join(RES, "stratifications", "{strat}.bed.gz"),
    params:
        url=lambda wc: f"{config['truth']['stratifications_url']}/{STRATS[wc.strat]}",
    log:
        os.path.join(LOGS, "resources", "stratification.{strat}.log"),
    shell:
        "curl -sSfL -o {output} {params.url} &> {log}"


def assembly_args(wildcards):
    """Where a haplotype comes from: its own file, or one unit (phase) of a diploid GenBank assembly"""
    assembly = ASSEMBLIES[wildcards.genome]
    if wildcards.hap in assembly:
        return f"--url {assembly[wildcards.hap]}"
    unit = assembly["units"][int(wildcards.hap[-1]) - 1]
    return f"--url {assembly['diploid']} --report {assembly['report']} --unit {unit}"


rule assembly_haplotype:
    """One haplotype of an assembly as plain FASTA (+ .fai)"""
    output:
        fasta=os.path.join(RES, "assemblies", "{genome}", "{hap}.fa"),
        fai=os.path.join(RES, "assemblies", "{genome}", "{hap}.fa.fai"),
    params:
        source=assembly_args,
    conda:
        os.path.join(ENVS, "reads.yaml")
    log:
        os.path.join(LOGS, "resources", "assembly.{genome}.{hap}.log"),
    shell:
        """
        python {SCRIPTS}/assembly_haplotype.py {params.source} --output {output.fasta} &> {log}
        samtools faidx {output.fasta}
        """
