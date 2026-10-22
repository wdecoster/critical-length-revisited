import os
import re


def output_path(path):
    """Relative output and resource paths are under config output_root (e.g. a study folder), absolute ones
    as given."""
    return path if os.path.isabs(path) else os.path.join(config.get("output_root", "."), path)


OUT = output_path(config["outdir"])
RES = output_path("resources")  # downloaded truth sets, stratifications, assemblies, reads
LOGS = os.path.join(OUT, "logs")
ENVS = os.path.join(workflow.basedir, "envs")
SCRIPTS = os.path.join(workflow.basedir, "scripts")
REF = output_path(config["reference"])

TRUTH_FILES = config["truth"]["files"]
STRATS = config["truth"]["stratifications"]
READ_SETS = config["read_sets"]
ASSEMBLIES = config["assemblies"]

LADDER = config["ladder"]
LADDER_GENOMES = LADDER["genomes"]
CHUNKS = [f"{i:03d}" for i in range(LADDER["chunks"])]


def sniffles1_settings(length):
    """Sniffles 1 settings run at a read length: the default everywhere, the others from their min_length"""
    return [name for name, s in LADDER["sniffles1"].items() if int(length) >= s["min_length"]]


def ladder_runs():
    """(genome, length, coverage, sniffles setting) of ladder steps 1-2"""
    return [
        (genome, length, coverage, setting)
        for genome in LADDER_GENOMES
        for length in LADDER["lengths"]
        for coverage in LADDER["coverages"]
        if coverage <= LADDER["simulated_coverage"][genome]
        for setting in sniffles1_settings(length)
    ]


def read_set_file(read_set, i):
    """A read set's i-th source: the downloaded copy, or the URL (streamed)"""
    url = READ_SETS[read_set]["urls"][int(i)]
    if READ_SETS[read_set]["download"]:
        return os.path.join(RES, "reads", read_set, os.path.basename(url))
    return url


wildcard_constraints:
    genome="|".join(map(re.escape, ASSEMBLIES)),
    hap="hap1|hap2",
    chunk=r"\d{3}",
    length=r"\d+",
    coverage=r"\d+",
    setting="|".join(map(re.escape, LADDER["sniffles1"])),
    read_set="|".join(map(re.escape, READ_SETS)),
    i=r"\d+",
    strat="|".join(map(re.escape, STRATS)),
