# Replication plan: SV discovery vs read length, 2026 edition

## Scope

- SV discovery only; phasing comes later.
- One genome, HG002.
- ONT R10.4.1 only: no HiFi comparator, no R9.
- The repo goes public together with a post on [Gigabase or gigabyte](https://wdecoster.github.io/gigabaseorgigabyte/).

## Questions

1. **Critical length per stratum (main question).** For each SV stratum, find the read length beyond which SV discovery no longer improves, in real ONT R10 data (simulation only as a fallback beyond 150 kb). Strata are type, size, and context (TR, SD, low mappability, other).
2. **Sequence completeness of insertions.** For insertions of a given size, what read length is needed to get the *full inserted sequence* right, not just the breakpoint?
3. **Translation to practice: when is long long enough?** How well does SV discovery do as a function of library N50, using realistic distribution shapes at a fixed coverage? Readers can act on this part most directly.

**Context, not a focus:** a short "then and now" comparison with the 2020 result (see below). It briefly sets the record straight where needed, such as the decline above 150 kb, without dwelling on the 2020 methods.

## Truth

- **HG002-Q100 v1.1**, a T2T diploid assembly. It is gap-free except for the centres of 9 of the 10 rDNA arrays.
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

**Coverage limits** (from ONT's QC length histograms of the two ULK flow cells: 165 Gb, read N50 107 kb, about 53x; to be confirmed by the inventory rule):

| Cut length | 20 kb | 30 kb | 50 kb | 75 kb | 100 kb | 150 kb | 200 kb | 300 kb |
|---|---|---|---|---|---|---|---|---|
| Coverage after cutting | 40x | 37x | 31x | 26x | 22x | 16x | 12x | 7x |

- **The grid:**
  - 30x up to 50 kb;
  - 15x up to 150 kb.
- **This changes the plan:** real reads cover lengths up to 150 kb at 15x, so simulation is only a fallback (see below).

**Length distributions (question 5):**
- Cut to an empirical length distribution instead of a fixed L. The shapes come from public data: the standard `giab_2025.01` library and the UL library.
- Each shape is rescaled to a grid of N50 values: 5, 8, 10, 15, 20, 30, 50 and 100 kb.
- Plot SV performance against library N50 at fixed coverage. This is the general "when is long long enough" recommendation, so no cohort-specific profiles.
- If a sentence helps place the recommendation, quote typical N50 values of real runs (public or aggregate, e.g. "standard brain-tissue libraries often have an N50 of 10–20 kb"). There is no separate analysis of our cohort.

### Fallback: simulation from Q100 v1.1 (only if needed)

The recommendation rests on real reads only. Simulation is a fallback, used only if an SV class (most likely segmental duplications) has not saturated at 150 kb in real data. In that case it extends that class's curve to 200–300 kb, as an extrapolation that is labelled as such.

- **Simulator, if needed:** squigulator (raw signal) followed by basecalling with Dorado. This gives the most realistic R10 errors, at the cost of GPU time. Alternatives:
  - Badread (`nanopore2023` model; slow at genome scale);
  - PBSIM3 (check for an R10.4.1 model).

  Any of them keeps the error rate constant across lengths.
- **Validation against real reads, per context:** identity, plus indel rates in homopolymers and TRs compared with non-repeat sequence, at lengths both arms share.
- **Simulate only the affected regions** (e.g. SDs with paralog partners), not the whole genome.
- **SimLoRD is used only for the then-and-now reproduction of 2020.**

### Length grid, coverage, replicates

- **Lengths:** 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 50, 75, 100, 150 kb (Arm A, as far as coverage allows), and 200 and 300 kb only through the simulation fallback.
- **Coverage:** 15x for the full grid (up to 150 kb), and 30x up to 50 kb. Downsample rather than realign.
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

## Then and now (context)

A small comparison that places the new recommendation next to the 2020 one. It is not a step-by-step attribution.

1. **The 2020 pipeline on HG00733,** as a check that the 2020 curve is reproduced (40x and 20x). Verify that GCA_003634875.1 is the phased haplotype set used in 2020.
2. **The same pipeline on HG002-Q100 v1.1 (20x),** which is gap-free apart from the rDNA arrays. This shows whether the decline above 150 kb came from the simulation. A side test runs Sniffles 1 with `--max_num_splits 1000` at 100 kb and above.
3. **The 2026 curve** from the main analysis (real R10 reads, current tools, GIAB truth) at the same coverage.

**Output:** one figure with the three curves, and one or two sentences in the post. These already run in stage 1 (steps 1–2); step 3 comes free with the main analysis.

## Pilot and compute

**Pilot region:** chr1, chr7, chr15, chr16 (SD-rich) and chr20 (well-behaved), about 20% of the genome.

**SD-stratum bias:** selecting reads by where their parent maps removes reads from paralogous copies on other chromosomes. Those are exactly the reads that would mismap into the pilot's SDs, so false positives in the SD stratum would be underestimated.
- In the pilot, also include parent reads that map to the SD paralog partners of the pilot SDs (from the GIAB SD stratification).
- Run genome-wide at 4 lengths (2, 10, 30, 100 kb) for the SD stratum.
- Pieces are always aligned to the whole genome.
- The simulation fallback, if used, follows the same approach.

**Estimate:**
- About 100 core-hours per 30x genome-wide alignment, less at 15x and in the pilot.
- `truvari refine` can take hours per call set.
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
docs/          plan.md, assembly_plan.md
notebooks/     figures for the blog
```

## Milestones

1. **Data inventory.**
   - Verify the LC2024 UL and `giab_2025.01` datasets: read-length distributions, basecaller versions, and achievable coverage per L. Fix the coverage grid from these.
   - Download and checksum the truth files and stratifications.
2. **Then and now, steps 1–2.** Reproduce 2020, then run the same pipeline on HG002.
3. **Pilot:** Arm A at 15x over the grid, plus the SD paralog extension.
4. **Assess,** then the genome-wide SD runs and the distribution profiles.
5. **Figures and blog draft.**

## Blog: core figures

The post leads with the recommendation; 2020 is context.

1. **SV performance against library N50 at fixed coverage:** the general recommendation of when reads are long enough, per SV class.
2. **A table of critical lengths with CIs per stratum,** plus the classes that keep improving.
3. **Real-read recall against L per stratum,** with the spanning-reads explanation.
4. **Then and now:** the 2020 curve, the same pipeline on HG002, and the 2026 curve, in one figure with a short caption.

## Decisions (2026-10-22)

- **HG002 only, ONT only.**
- **Focus on the new recommendation.** The 2020 paper is context: set the record straight briefly where needed, without emphasising its limitations.
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
- Whether the simulation fallback is needed: decided after the pilot, only for SV classes that have not saturated at 150 kb.
