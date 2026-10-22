# Replication plan: SV discovery vs read length, 2026 edition

## Scope

- SV discovery only; phasing comes later.
- One genome, HG002.
- ONT R10.4.1 only: no HiFi comparator, no R9.
- The repo goes public together with a post on [Gigabase or gigabyte](https://wdecoster.github.io/gigabaseorgigabyte/).

## Questions

1. **Critical length per stratum (main question).** For each SV stratum, find the read length beyond which SV discovery no longer improves, in real ONT R10 data, with simulation as support. Strata are type, size, and context (TR, SD, low mappability, other).
2. **Which old conclusions were artifacts?** How much of the 2020 curve came from:
   - the truth set;
   - the length-dependent CLR error model;
   - the tools;
   - the lenient matching?

   This is the "ladder".
3. **Sequence completeness of insertions.** For insertions of a given size, what read length is needed to get the *full inserted sequence* right, not just the breakpoint?
4. **The >150 kb decline.** Does it reproduce, and is it caused by the simulator's coverage near contig ends and gaps, or by Sniffles 1's split filter?
5. **Translation to practice: when is long long enough?** How well does SV discovery do as a function of library N50, using realistic distribution shapes at a fixed coverage? Readers can act on this part most directly.

## Truth

- **HG002-Q100 v1.1**, a T2T diploid assembly with no gaps.
  - The GIAB v5.0q benchmark is built from v1.1, so we simulate from v1.1 as well. The current version is v1.2 (Aug 2026), but reads and truth must agree.
- **GIAB HG002 v5.0q SV benchmark, the primary truth.** Files: `HG002_GRCh38_v5.0q_stvar.vcf.gz` with `HG002_GRCh38_v5.0q_stvar.benchmark.bed`.
  - It is labelled a *Draft* and was made with dipcall (`-z200000,10000`) from Q100 v1.1. Pin the exact files with checksums.
  - It does **not** exclude complex SVs.
  - It excludes SDs, satellites and TRs > 10 kb *only where the assembly-to-reference alignment breaks*, plus inversion regions (~37 Mb), gaps and their flanks.
  - So part of the SD content is inside the benchmark regions, and SDs can be evaluated as a stratum, with the n reported.
- **GIAB genomic stratifications (v3.6) for GRCh38.** These cover tandem repeats and homopolymers, segmental duplications, low mappability, the MHC, and difficult regions.
- **Outside the benchmark regions:** GIAB's own `GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.{vcf.gz,bed}`, reported as "assembly concordance", not as truth. We do not run our own dipcall/PAV.
- **chrX and chrY:** included, and evaluated as haploid, as in the GIAB benchmark. Report them separately if they behave differently.

## Reads

### Arm A: real reads cut to length (main evidence)

**Data:**

| Dataset | Use | Notes |
|---|---|---|
| ONT London Calling 2024 HG002 ULK114, 2 flow cells, sup v5.0.0 (`s3://ont-open-data/londoncalling2024/`) | Main source | Reported by the review as about 165 Gb mapped, read N50 about 107 kb, roughly 28x in reads ≥ 100 kb. *To verify in milestone 1.* |
| ONT `giab_2025.01` HG002, standard LSK114, 2 flow cells | Calibration at L ≤ 30 kb (see below) | |
| `gm24385_2023.12` UL R10 (~40x, N50 ~91 kb) | Backup only | Older basecaller. Never merge with the main source. |

**Cutting:**
- For target length L, cut each read of at least L into ⌊len/L⌋ pieces of length L.
- Place the pieces at random within the parent read, spreading the leftover bases randomly over the gaps. This wastes fewer bases than cutting from a random offset and discarding the remainder.
- Give pieces unique names (`parent_k`), so no tool merges evidence across pieces.
- Subsample to the target coverage, measured as aligned bases within the benchmark BED.

What this does and does not change:
- It keeps the real error profile, coverage bias and chimeras, and changes only the length. Interior cut points are placed at random relative to breakpoints, just like natural fragmentation.
- Pieces lack *real* read-end artifacts: low-quality starts, adapter or transposase remnants, and fragmentation sequence bias. This slightly favours short L, and the post should say so.

**One molecule pool:**
- All L come from the same UL dataset, so library or flow cell is not confounded with L.
- To check that the UL library behaves like a standard one, compare it against `giab_2025.01` at 10, 20 and 30 kb: cut both, and compare the curves.

**Coverage limits:**
- About 30x reaches roughly 50 kb.
- 15–20x reaches roughly 100 kb.
- The grid is run at 30x up to 50 kb and at 15x over the full range. That gives one consistent 15x curve, plus a 30x curve that is reliable up to 50 kb.
- Exact limits are computed in milestone 1 from the read-length distribution.

**Length distributions (question 5):**
- Cut to an empirical length distribution instead of a fixed L. The shapes come from public data: the standard `giab_2025.01` library and the UL library.
- Each shape is rescaled to a grid of N50 values: 5, 8, 10, 15, 20, 30, 50 and 100 kb.
- Plot SV performance against library N50 at fixed coverage. This is the general "when is long long enough" recommendation, so no cohort-specific profiles.
- If a sentence helps place the recommendation, quote typical N50 values of real runs (public or aggregate, e.g. "standard brain-tissue libraries often have an N50 of 10–20 kb"). There is no separate analysis of our cohort.

### Arm B: simulation from Q100 v1.1 (support and control)

Purpose:
- an exact fixed-length grid without coverage limits;
- the ladder;
- the >150 kb question on a gap-free genome.

Simulator:
- PBSIM3 (check that it ships an R10.4.1 error model) or Badread (`nanopore2023` model). Badread is slow at genome scale, and its chimera, junk and glitch rates must be set deliberately.
- Validate against Arm A **per context**: identity, plus indel rates in homopolymers and TRs compared with non-repeat sequence. Neither simulator fully reproduces R10's systematic errors, and this is the main reason Arm A is the primary evidence.

Lengths: the Arm A grid plus 300 kb. 500 kb and 1 Mb are used only in the ladder.

### Length grid, coverage, replicates

- **Lengths:** 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 50, 75, 100, 150 kb (Arm A, as far as coverage allows), and 200 and 300 kb (Arm B).
- **Coverage:** 15x for the full grid, and 30x up to 50 kb. Downsample rather than realign.
- **Replicates:** 3 random cut placements or seeds at 3 lengths (5, 20, 100 kb). Other uncertainty comes from a block bootstrap (1 Mb blocks).

## Analysis

- **Alignment:** minimap2 (≥ 2.28, pinned) with `-x lr:hq`.
- **Calling:**
  - Sniffles2 is the main caller, with pinned version and coverage-adaptive defaults. Its MAPQ filter shapes behaviour in SDs, so record it.
  - cuteSV runs at 3 lengths, as a check that the result does not depend on one tool. It uses ONT-recommended clustering settings, and a minimum support scaled to coverage (its default `-s` is fixed, the 2020 flaw again).
- **Assembly arm:** undecided. My recommendation is to leave it out of the first post; it could become a follow-up.

## Evaluation

- **Primary scoring:** GIAB's recommended setup, `truvari bench --pick ac --passonly -r 2000 -C 5000` followed by `truvari refine`. First remove `ALT="*"` records. Genotype concordance is reported too.
- **Lenient secondary scoring:** for questions where only the event matters, `--dup-to-ins --pctseq 0`. DUP calls have no sequence, so `--dup-to-ins` only makes sense without sequence matching.
- **The 2020 matcher on the same calls:** SURVIVOR site-based, type only, 500 bp. This shows the leniency directly.
- **Strata:**
  - type: DEL, INS, INV;
  - size: 50–100, 100–300 (Alu), 300 bp–1 kb, 1–6 kb, 6–10 kb (L1), 10–50 kb, ≥ 50 kb;
  - context: TR, SD, low mappability, other.

  Report n per stratum. Strata with few events (≥ 50 kb, INV, SD) are merged or shown with CIs only.
- **Recall and precision separately,** not only F1.
- **Mechanistic view:** recall against the number of reads that fully span the SV plus 1 kb of flank. This is cheap and explains *why* each stratum saturates where it does.
- **Insertion completeness:**
  - Strict matching censors truncated insertions: they become an FN plus an FP. So completeness is computed from a relaxed match (`-P 0 -O 0`, or the nearest insertion within 1 kb).
  - Measures: ALT length relative to the truth, and sequence identity.
  - Shown as a function of insertion size relative to L.
- **Critical length:**
  - Fit a saturating curve per stratum (isotonic, or log-logistic in log L).
  - The critical length is the smallest L whose recall (and F1) is not significantly below the plateau. Test this with a *paired* block bootstrap across L: the same blocks for every L.
  - Report it as a table with 95% CIs. A plain "maximum" plateau would bias the critical length upwards.

## The ladder (question 2)

Each step changes one thing.
- **Coverage:** constant at 20x throughout.
- **Lengths:** 1, 5, 20, 50 and 100 kb (step 6 is limited by real-read coverage), plus 300 and 500 kb for steps 1–5.
- **Output:** recall, precision and F1 per step. Five lengths are too few for a critical length.

1. **The original pipeline on HG00733,** reproducing Fig. 1 as a sanity check. Verify that GCA_003634875.1 is the phased haplotype set used in 2020.
2. **Same pipeline, HG002-Q100 v1.1 instead of HG00733,** with the 2020 matcher. This step changes the truth and the read source. The genome is gap-free, which tests the gap explanation for the >150 kb decline.
3. **Evaluation changed to truvari,** with benchmark regions and one-to-one, sequence-aware matching.
4. **Reads changed to an R10 error model,** with constant error across lengths. This removes both the CLR error model and the confounding of error with length.
5. **Tools changed** to minimap2 `lr:hq` and Sniffles2.
6. **Real reads** (Arm A).

Notes:
- **Attribution depends on the order of the steps.** To show how large that effect is, a 2×2 (error model × tools) is run at 5 and 50 kb.
- **Side tests at step 1:**
  - Sniffles 1 with `--max_num_splits 1000`;
  - per-SV coverage against the distance to a contig end or gap.

## Pilot and compute

**Pilot region:** chr1, chr7, chr15, chr16 (SD-rich) and chr20 (well-behaved), about 20% of the genome.

**SD-stratum bias:** selecting reads by where their parent maps removes reads from paralogous copies on other chromosomes. Those are exactly the reads that would mismap into the pilot's SDs, so false positives in the SD stratum would be underestimated.
- In the pilot, also include parent reads that map to the SD paralog partners of the pilot SDs (from the GIAB SD stratification).
- Run genome-wide at 4 lengths (2, 10, 30, 100 kb) for the SD stratum.
- Pieces are always aligned to the whole genome.
- Arm B uses the same approach.

**Estimate:**
- About 100 core-hours per 30x genome-wide alignment, less at 15x and in the pilot.
- `truvari refine` can take hours per call set.
- Badread can be slow, so prefer PBSIM3 if it validates.
- Storage: about 30–50 GB per 30x CRAM, roughly 3–4 TB in total if everything is kept. Keep VCFs and summaries; delete the cut FASTQs and the BAMs/CRAMs after calling, except at a few lengths kept for figures.

## Workflow

A Snakemake workflow in this repo, with:
- a server profile;
- an `output_root` outside the repo;
- `scripts/setup_server.sh`;
- a pinned conda env per tool.

Server paths only go in an uncommitted `config/server.yaml`. The read-cutting step is written as a small standalone tool and released with the post.

```
config/        config.yaml (grid, public data sources, truth versions + checksums)
workflow/      rules: data, cut_reads, simulate, align, call, evaluate, ladder, plots
docs/          assessment.md, plan.md, blog draft
notebooks/     figures for the blog
```

## Milestones

1. **Data inventory.**
   - Verify the LC2024 UL and `giab_2025.01` datasets: read-length distributions, basecaller versions, and achievable coverage per L. Fix the coverage grid from these.
   - Download and checksum the truth files and stratifications.
2. **Ladder steps 1–2.** Reproduce 2020, then swap in HG002. This cheaply checks the main critique and the >150 kb explanation.
3. **Pilot:** Arm A at 15x over the grid, the SD paralog extension, simulator validation, and Arm B at selected lengths.
4. **Assess,** then the genome-wide SD runs, the remaining ladder steps, and the distribution profiles.
5. **Figures and blog draft.**

## Blog: core figures

1. **The ladder:** which parts of the old curve were real and which were artifacts, including the confounding of error with length.
2. **Real-read recall against L per stratum,** with the spanning-reads explanation.
3. **A table of critical lengths with CIs per stratum,** plus the classes that keep improving.
4. **SV performance against library N50 at fixed coverage:** the general recommendation of when reads are long enough, per SV class.

## Decisions (2026-10-22)

- **HG002 only, ONT only.**
- **No cohort-specific analysis.** The recommendation is general, as a function of N50. Aggregate read-length statistics may be quoted for context; publishing them is fine.
- **Assembly-based calling undecided;** decide after the pilot.
- **Public repo.** It contains:
  - only public data;
  - no cohort sample IDs;
  - no references to unpublished results or internal projects;
  - no server paths.

  Before publication, grep the full git history for sample-ID prefixes.

## Still open

- Verification of the UL dataset and the final coverage grid (milestone 1).
- The assembly arm.
