# Replication plan: SV discovery vs read length, 2026 edition

Scope:
- SV discovery only; phasing comes later.
- One genome, HG002.
- ONT R10.4.1 only: no HiFi comparator, and no R9.
- The repo goes public together with the blog post on [Gigabase or gigabyte](https://wdecoster.github.io/gigabaseorgigabyte/) (source in `~/repos/gigabaseorgigabyte`, Hugo).

## Questions

1. **Critical length per stratum (main question).** For each SV stratum, find the read length beyond which SV discovery no longer improves, both in real ONT R10 data and in simulation. Strata are type, size, and context (TR, SD, low mappability, other).
2. **Which old conclusions were artifacts?** How much of the 2020 curve came from the truth set, the CLR error model, the tools, and the lenient matching? This is the "one change at a time" ladder, built for the blog.
3. **Sequence completeness of insertions.** For insertions of a given size, what read length is needed to get the *full inserted sequence* right, not just the breakpoint? This links to the ~18.8 kb RA insertion in TMEM106B.
4. **The >150 kb decline.** Does it reproduce, and is it caused by the simulator's coverage near gaps or by Sniffles 1's `--max_num_splits`?
5. **Translation to practice.** What do the curves mean for the read-length distributions we actually get? That means our cohort, standard and ultra-long libraries, and degraded DNA. See "Cohort read lengths" below.

## Truth: HG002

- **HG002-Q100 v1.1:** a T2T diploid assembly with no gaps. This removes issue 1 (read source and truth sharing the same holes).
- **GIAB v5.0q SV benchmark (`v5.0q_stvar`) on GRCh38, with its confident-region BED.** This is the primary truth. It is still labelled a *draft*, so check the release notes and pin the exact version. It excludes very large and complex SVs and CNVs.
- **GIAB genomic stratifications (v3.6 or later) for GRCh38.** These cover tandem repeats and homopolymers, segmental duplications, low mappability, the MHC, difficult regions, and others.
- **Secondary truth outside the confident regions:** our own dipcall/PAV-style calls of Q100 against GRCh38. These are needed to say *anything* about segmental duplications, where long reads should matter most. They are reported separately and labelled as less certain.

## Reads: two arms

### Arm A, real reads cut to length (main evidence)

1. Take the longest HG002 R10.4.1 sup dataset(s) available. Candidates:
   - ONT open data `giab_2025.01`, which is already used in the Tier 2 benchmark;
   - an ONT ultra-long Kit 14 HG002 run;
   - any HG002 UL R10 data used for Q100 / HPRC.
   
   I still need to check what exists, its coverage above 100 kb, and the basecaller version.
2. For a target length L, cut every read of at least L into consecutive pieces of L, starting at a random offset, and discard the remainder. Then subsample to the target coverage.
   - This keeps the real error profile, real coverage bias, and real chimeras, and changes only the length.
   - A length filter would instead select reads of at least L and confound length with coverage.
3. **Limit:** coverage at L is bounded by the total read length in reads of at least L. It must be measured to see how far the real arm reaches, probably 100–150 kb at 30x if an UL dataset is used. Beyond that we rely on Arm B.
4. **Real length distributions:** the same cutting can reproduce empirical distributions as well as fixed lengths. Proposed profiles:
   - our cohort's R10 runs (see below);
   - a standard LSK114 run;
   - an UL run;
   - a "degraded" profile with N50 around 8 kb.

### Cohort read lengths (`cramino_metrics.tsv`, 489 samples)

The cramino summary gives yield, N50, median and mean length, and yield above 25 kb. "R10-like" below means a modal identity of at least 98.5%.

| | R10-like (n=132) | R9-like (n=357) |
|---|---|---|
| N50, P10 / median / P90 (kb) | 8.9 / 16.1 / 21.1 | 14.4 / 19.9 / 24.7 |
| Median read length, median (kb) | 6.8 | 11.6 |
| Yield above 25 kb, median (Gb) | 17.6 | 28.8 |
| Mean coverage, median | 38.7x | 34.1x |

- **Our cohort sits right at or below the 2020 threshold.** The R10 runs, mostly frontal cortex, have shorter reads than the R9 runs. A tenth of them have an N50 below 9 kb. So the replication answers a concrete question for our own data.
- **In the cohort, yield does not correlate with N50 (r = -0.01).** Within standard preparations, longer reads did not cost yield. Our data has no ultra-long runs, so the trade-off for ultra-long libraries cannot be estimated from it. If needed, we would take ultra-long yields from ONT or literature numbers and treat them as assumptions.
- **The summary statistics are not enough to reproduce the shape of a distribution.**
  - A lognormal fitted from median and mean underestimates N50 by about 24% (the observed tail is heavier).
  - The plan is to define three profiles from the summary (cohort R10 at the P10, median and P90 of N50).
  - For one representative sample per profile, we get the full read-length histogram with a read-only `cramino --hist` on the server, and replay it in Arm A.
  - Fallback: a lognormal with parameters from the median and N50.

### Arm B, simulation from the Q100 assembly (extends the length range, gives control)

- **Simulator:** PBSIM3 with an ONT R10 error model, or Badread. Choose after a quick check that the simulated read identity and error mix (identity, homopolymer indels) match the Arm A reads after alignment.
- **Lengths:** the same grid as Arm A plus 300 kb, 500 kb and 1 Mb. Both fixed lengths and lognormal distributions, simulated from both haplotypes.
- **Gaps:** because the Q100 assembly has none, coverage near contig ends is only an issue at chromosome ends, which tests the explanation for the >150 kb decline in issue 2.

### Length grid and coverage

- **Lengths:** 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 50, 75, 100, 150 and 200 kb. Arm B adds 300 kb, 500 kb and 1 Mb.
- **Coverage:** 30x (main) and 15x (clinical and cohort-like). Downsampling from the 30x BAM avoids realignment.
- **Replicates:** 3 random offsets or seeds at 1, 5, 20, 100 and 500 kb, and 1 elsewhere. Uncertainty comes from bootstrapping over 1 Mb blocks.

## Analysis

- **Alignment:** a current minimap2 with `-x lr:hq`.
- **Calling:** Sniffles2 is the main caller (current release, coverage-adaptive defaults). cuteSV is the second caller, to show the conclusion does not depend on one tool. Severus is optional.
- **Assembly arm (undecided; my recommendation is to leave it out of the first post):** hifiasm in ONT mode or Verkko on Arm A data at a few lengths (10, 30, 100 kb), with SVs called from the assembly. This answers at what read length assembling beats calling from alignments, and for which classes.

## Evaluation

- **truvari bench** in the GIAB confident regions, with `--dup-to-ins` and default sequence and size similarity. Repeats are harmonised with `truvari refine`. Genotype concordance is reported too.
- **The 2020 matcher, run on the same calls:** type only, 500 bp, no size or sequence check. This shows the leniency effect directly.
- **Strata:**
  - type: DEL, INS, INV, DUP→INS;
  - size: 50–100, 100–300 (Alu), 300 bp–1 kb, 1–6 kb, 6–10 kb (L1), 10–50 kb, ≥ 50 kb;
  - context: TR, SD, low mappability, other.
- **Insertion completeness:** for true-positive insertions, the sequence identity of the called ALT to the truth, and the fraction that are full-length (identity ≥ 0.95 and length within 5%), plotted against insertion size / L.
- **Breakpoint and size accuracy** of true positives.
- **Critical length:** the smallest L with an F1 of at least 99% (and 95%) of that stratum's plateau, with a bootstrap 95% CI. Reported as a table of L per stratum, which is the main figure for the blog.

## The ladder (question 2)

Each step changes one thing, at 5 selected lengths (1, 5, 20, 100, 500 kb):

1. **The original pipeline on HG00733**, reproducing Fig. 1. The original assembly, GCA_003634875.1, is still available. This is a sanity check that we can reproduce 2020.
2. **Same pipeline, HG002-Q100 instead of HG00733,** scored with the 2020 matcher. This step changes the truth and the read source.
3. **Change the evaluation to truvari**, with confident regions and sequence matching.
4. **Change the reads to the R10 error model** (Arm B simulation).
5. **Change the tools** to current minimap2 `lr:hq` and Sniffles2.
6. **Change to real reads** (Arm A).

At each step, report the critical length (all SVs) and the >150 kb behaviour. Also test step 1 with Sniffles 1 `--max_num_splits 1000`, and plot per-SV coverage against distance to a contig gap.

## Compute (rough, to check on a pilot)

- **Main grid:** about 15 lengths × 2 arms × 2 coverages (downsampled) × (1–3 replicates), so roughly 80 alignments of 30x genome-wide data.
- **Cost of a full genome:** at about 100 core-hours per 30x alignment, the full genome is about 8–10k core-hours.
- **Pilot first, on 20% of the genome:** chr1, chr7, chr15, chr16 (SD-rich) and chr20 (well-behaved).
  - In Arm A, select reads whose *parent* read maps to these chromosomes, but align the pieces to the whole genome so that mismapping is still possible.
  - Arm B simulates only these chromosomes but aligns to the whole genome.
- **Decision point:** after the pilot, decide whether the genome-wide run is needed for the blog. Probably only for the SD stratum.
- **Storage:** CRAM, and delete the cut FASTQs after alignment, as in the Tier 2 benchmark.

## Workflow

A Snakemake workflow in this repo, following the regional-pangenome and ont-phasing-benchmark conventions:

- server profile;
- `output_root` in the study folder;
- `setup_server.sh`;
- a pinned conda env per tool.

Proposed layout:

```
config/        config.yaml (grid, data sources, truth versions)
workflow/      rules: data, cut_reads, simulate, align, call, truth, evaluate, ladder, plots
docs/          assessment.md, plan.md, blog draft
notebooks/     figures for the blog
```

## Milestones

1. **Data inventory:**
   - HG002 R10 datasets, their read-length distributions, and the achievable coverage per L.
   - Pin the truth versions.
   - Get three cohort read-length histograms.
2. **Ladder steps 1–2:** reproduce 2020, then swap the truth. This is a quick, cheap way to check that the critique holds.
3. **Pilot grid, Arms A and B, on 20% of the genome.**
4. **Assess, then decide** on running genome-wide and on the assembly arm.
5. **Figures and blog draft.**

## Proposed blog outline

1. What the 2020 paper said and why the question still matters.
2. What we can do now that we couldn't then: a T2T diploid truth, R10 accuracy, and current tools.
3. The ladder: which parts of the old curve were real and which were artifacts.
4. The new answer: critical length per SV class, plus the classes that keep improving.
5. Practical guidance on library choice, ultra-long or not, and degraded samples.

## Decisions (2026-10-22)

- **HG002 only, ONT only.**
- **Cohort read lengths come from `cramino_metrics.tsv`,** plus histograms of three representative samples.
- **Assembly-based calling is undecided.** It is a different question (what read length is needed to *assemble* rather than *call*), and expensive. Decide after the pilot; it could become a follow-up post.
- **The blog appears on Gigabase or gigabyte, and the repo becomes public.** Because the repo will be public:
  - only public data;
  - no cohort sample IDs, only aggregate statistics of read length;
  - no server paths in committed config (keep those in a `server.yaml` that is not committed).

## Still open

- **Which HG002 R10 datasets to use** (milestone 1). Arm A needs enough coverage in reads of 100 kb and longer.
- **Whether a genome-wide run is needed** after the pilot.
