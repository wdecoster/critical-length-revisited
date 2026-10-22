# Replication plan: SV discovery vs read length, 2026 edition

Scope: SV discovery only (phasing comes later). ONT R10.4.1 only, since R9 is no longer relevant for new work.

## Questions

1. **Critical length per stratum (main question).** For each SV stratum, find the read length beyond which SV discovery no longer improves, both in real ONT R10 data and in simulation. Strata are type, size, and context (TR, SD, low mappability, other).
2. **Which old conclusions were artifacts?** How much of the 2020 curve came from the truth set, the CLR error model, the tools, and the lenient matching? This is the "one change at a time" ladder, built for the blog.
3. **Sequence completeness of insertions.** For insertions of a given size, what read length is needed to get the *full inserted sequence* right, not just the breakpoint? This links to the ~18.8 kb RA insertion in TMEM106B.
4. **The >150 kb decline.** Does it reproduce, and is it caused by the simulator's coverage near gaps or by Sniffles 1's `--max_num_splits`?
5. **Translation to practice (optional).** What do the curves mean for real read-length distributions (sheared, standard, ultra-long, degraded DNA)? Does the lower yield of ultra-long libraries cancel out their benefit at a fixed flow-cell budget?

## Truth: HG002

- **HG002-Q100 v1.1:** a T2T diploid assembly with no gaps. This removes issue 1 (read source and truth sharing the same holes).
- **GIAB v5.0q SV benchmark (`v5.0q_stvar`) on GRCh38, with its confident-region BED.** This is the primary truth. It is still labelled a *draft*, so check the release notes and pin the exact version. It excludes very large and complex SVs and CNVs.
- **GIAB genomic stratifications (v3.6 or later) for GRCh38.** These cover tandem repeats and homopolymers, segmental duplications, low mappability, the MHC, difficult regions, and others.
- **Secondary truth outside the confident regions:** our own dipcall/PAV-style calls of Q100 against GRCh38. These are needed to say *anything* about segmental duplications, where long reads should matter most. They are reported separately and labelled as less certain.
- **Possible second genome:** HG00733 (HPRC release 2 assembly, plus ONT UL data). It gives continuity with the original paper. There is no GIAB benchmark for it, so its assembly-derived truth would be weaker. *To discuss.*

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
   - our cohort's R10 runs;
   - a standard LSK114 run;
   - an UL run;
   - a "degraded" profile with N50 around 8 kb.

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
- **Assembly arm (stretch goal, expensive):** hifiasm in ONT mode or Verkko on Arm A data at a few lengths (10, 30, 100 kb), with SVs called from the assembly. This answers at what read length assembling beats calling from alignments, and for which classes.

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

1. **Data inventory:** HG002 R10 datasets, their read-length distributions, and the achievable coverage per L. Pin the truth versions.
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

## Open points to discuss

- Is HG002 alone enough, or do we add HG00733 for continuity?
- What real data do we use? Can we include our own cohort's read-length distributions and our flow-cell yields per protocol? These are needed for question 5.
- Is a HiFi comparator in or out? It is not ONT, but it isolates the effect of the error rate from read length.
- Is the assembly arm in scope?
- Where should the blog appear? Should this repo become public at publication time?
