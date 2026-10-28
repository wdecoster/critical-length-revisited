# critical-length-revisited

This repo revisits [De Coster et al. 2020, *Critical length in long-read resequencing*](https://doi.org/10.1093/nargab/lqz027) (NAR Genom Bioinform 2(1):lqz027) with current data and tools, for a follow-up post on [Gigabase or gigabyte](https://wdecoster.github.io/gigabaseorgigabyte/). The first focus is structural variant discovery; phasing comes later.

- [docs/assessment.md](docs/assessment.md): notes on the 2020 SV methods, and what to change for the replication.
- [docs/plan.md](docs/plan.md): the replication plan. It uses the HG002 Q100 truth, ONT R10 reads (real reads cut to length, plus simulation), Sniffles2 and cuteSV, and stratified truvari evaluation.

The original code is at https://github.com/wdecoster/read_length_SV_discovery.

## Workflow (stage 1)

The Snakemake workflow currently covers milestones 1–2 of the plan:

- **Data inventory:**
  - downloads the HG002 ultra-long R10 read set;
  - makes read-length histograms of the ultra-long and standard read sets;
  - computes the coverage left after cutting reads to each length (`results/inventory/`).
- **Truth:** downloads and checksums the GIAB HG002 v5.0q SV benchmark, the GIAB dipcall output and the stratifications.
- **Then and now, steps 1–2:** the 2020 pipeline, unchanged, on the 2020 HG00733 assembly and on HG002-Q100 v1.1:
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


## Workflow (stage 2: the pilot)

`snakemake ... pilot` runs Arm A of the plan on the pilot chromosomes (chr1, chr7, chr15, chr16, chr20):

- **Parent reads:** HG002 R10 reads with an alignment on the pilot chromosomes, or on the paralog partners elsewhere of the pilot's segmental duplications (UCSC genomicSuperDups). The ultra-long set is mapped first (minimap2 `lr:hq`, positions only); the standard set brings its own alignments.
- **Cutting:** `workflow/scripts/cut_reads.py` cuts every read of at least L into floor(len/L) pieces of exactly L, placed at random within the read (the leftover spread over the gaps), named `<read>_<k>`. It is a standalone tool:

  ```bash
  python workflow/scripts/cut_reads.py reads.fastq.gz --length 20000 --seed 1 --output pieces.fastq.gz
  ```
- **Alignment and coverage:** pieces are aligned to the whole genome (minimap2 2.31 `lr:hq`) and subsampled to exact coverages (30x up to 50 kb, 15x up to 150 kb) measured in the GIAB benchmark regions of the pilot chromosomes. A length with less coverage available uses all reads and records what it got.
- **Calling:** Sniffles2 2.8.1 at every length, cuteSV 2.1.4 at 5, 20 and 100 kb (ONT clustering settings, minimum support scaled to coverage).
- **Scoring:** truvari 5.5.0 against GIAB v5.0q (ALT=`*` removed): strict (GIAB's setup, then `truvari refine`) and lenient (`--dup-to-ins --pctseq 0`). Recall, precision and F1 per stratum (type, size, context: SD, TR, low mappability, other), and per truth SV the reads spanning it with 1 kb flanks.
- **Replicates:** three cut placements at 5, 20 and 100 kb; the standard (`giab_2025.01`) set at 10, 20 and 30 kb as calibration of the ultra-long library.

Output: `results/pilot/summary.tsv` (one row per run and stratum), plus per run the truvari results and spanning counts. Alignments are deleted once called (`pilot.keep_bams_lengths` keeps some as `*.kept.bam`).


## Workflow (milestone 4: library N50)

`snakemake ... pilot_distributions` cuts the same ultra-long parent reads into pieces whose lengths follow a real library's read length distribution (the inventory histograms of the standard and the ultra-long library), rescaled to N50 values of 5, 8, 10, 15, 20, 30, 50 and 100 kb. The runs (`L N<n50><shape>`, e.g. `N20000standard`) then go through the pilot's alignment, subsampling, calling and scoring.

- **Cutting:** `cut_reads.py --distribution <histograms> --n50 <bp>` fills each read with pieces of drawn lengths until the next does not fit, and spreads the leftover over the gaps. Long pieces only fit in long reads, so the drawing weights are first calibrated on the reads' own lengths (per length bin: target / achieved share), and the bins rescaled until the pieces' N50 is within 1% of the target. On the ultra-long reads the standard shape reaches every N50; the ultra-long shape at 100 kb reaches about 90 kb (the parents run out of longer reads), which the stats record.

  ```bash
  python workflow/scripts/cut_reads.py reads.fastq.gz --distribution lengths.tsv --n50 20000 --seed 1 --output pieces.fastq.gz
  ```
- **Coverage:** 30x up to an N50 of 50 kb, 15x at every N50.

Output: `results/pilot/distributions.tsv`, the strata of every run with its shape, target N50 and the N50 and mean length achieved. Plot against the achieved N50.
