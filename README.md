# critical-length-revisited

This repo revisits [De Coster et al. 2020, *Critical length in long-read resequencing*](https://doi.org/10.1093/nargab/lqz027) (NAR Genom Bioinform 2(1):lqz027) with current data and tools, for a follow-up post on [Gigabase or gigabyte](https://wdecoster.github.io/gigabaseorgigabyte/). The first focus is structural variant discovery; phasing comes later.

- [docs/assessment.md](docs/assessment.md): a critical assessment of the original SV methods.
- [docs/plan.md](docs/plan.md): the replication plan. It uses the HG002 Q100 truth, ONT R10 reads (real reads cut to length, plus simulation), Sniffles2 and cuteSV, and stratified truvari evaluation.

The original code is at https://github.com/wdecoster/read_length_SV_discovery.

## Workflow (stage 1)

The Snakemake workflow currently covers milestones 1–2 of the plan:

- **Data inventory:**
  - downloads the HG002 ultra-long R10 read set;
  - makes read-length histograms of the ultra-long and standard read sets;
  - computes the coverage left after cutting reads to each length (`results/inventory/`).
- **Truth:** downloads and checksums the GIAB HG002 v5.0q SV benchmark, the GIAB dipcall output and the stratifications.
- **Ladder steps 1–2:** the 2020 pipeline, unchanged, on the 2020 HG00733 assembly and on HG002-Q100 v1.1:
  - SimLoRD 1.0.2;
  - minimap2 2.14 `map-pb`;
  - Sniffles 1.0.10;
  - a paftools truth set;
  - surpyvor-style SURVIVOR matching.

  Output: `results/ladder/summary_steps1-2.tsv`.

Setup on a server (conda/mamba required; GRCh38 no-alt analysis set as reference):

```bash
scripts/setup_server.sh -r /path/to/GRCh38_no_alt.fa -o /path/to/output_root
mamba activate critical-length-revisited
snakemake --configfile config/server.yaml --profile profiles/server -n
```

Useful targets:
- `inventory`: the read-set inventory only.
- `ladder_summary`: the ladder.

The ladder simulates ~1.3 Tb of reads in total. SimLoRD runs at about 2 Mb/s per job, and the CLR-like alignments with minimap2 2.14 are the slowest step. Alignments are deleted once called, unless `ladder.keep_bams` is set.

