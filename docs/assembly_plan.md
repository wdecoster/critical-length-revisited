# Follow-up plan: de novo assembly vs read length

Status: plan only (2026-10-31). Waits until the server has room; see "Order and gating".

The SV post answered "how long should reads be for read-based SV discovery?". This follow-up asks the same
question for de novo assembly, which is a different problem. Assembly also needs coverage for error correction.
It needs reads that span repeats to stay contiguous and in the right copy. And it needs variant density along
the reads to separate the haplotypes. Read length should therefore matter more, and for longer, than it did for
SV calling.

Same principles as the SV study:
- Real ONT R10.4.1 HG002 reads cut to length.
- Fixed lengths for mechanism, library N50 for the recommendation.
- HG002 Q100 v1.1 as truth.
- Results reported per genomic context.

## Questions

1. **Contiguity.** How do the contig count, contig N50/NG50, auN and the largest contig change with library N50
   and coverage? Where does the gain level off?
2. **Correctness.** How do misassemblies, NGA50 (contiguity after breaking at misassemblies), base accuracy (QV)
   and collapsed or duplicated sequence change with read length? Are short-read assemblies just more fragmented,
   or also more often wrong?
3. **Haplotype resolution.** Without trio or Hi-C data, how far do the reads alone phase the assembly? Measures:
   switch and Hamming error rates, and the phased block N50, all against the two Q100 haplotypes.
4. **Assembly-based SV discovery.** Score SVs called from the assemblies against the same GIAB benchmark and
   strata as the read-based calls. At the same N50 and coverage, when does assembly beat read-based calling?
   Insertions over 1 kb and segmental duplications are the obvious candidates. This also answers the open
   follow-up in the first post (completeness of inserted sequence).
5. **Targeted (regional) assembly.** For a window of about 0.5 Mb, what N50 and coverage give a complete,
   correct assembly of each haplotype from the reads that map there? This is how regional pangenome studies
   assemble loci across many samples (reads extracted per region and per haplotype, then assembled). Per
   haplotype, the coverage is half the sample's.

Q1–Q4 are about the genome (chromosome) scale; Q5 is about the regional scale. Q5 is the cheapest arm and the
one most relevant for regional pangenomics, so it should go first.

## Reads and conditions

- **Source:** the parent reads of the pilot (ONT London Calling 2024 HG002 ultra-long, 53x, N50 106 kb), cut
  with `cut_reads.py`:
  - **N50 mode (main arm):** the standard shape, N50 5, 10, 20, 30 and 50 kb.
  - **Fixed lengths (mechanism):** 10, 20, 50 and 100 kb, where the coverage allows.
- **Coverage:**
  - **Genome-scale arms (Q1–Q4):** 30x as the main coverage. hifiasm `--ont` relies on read overlaps for error
    correction, so 15x is expected to be poor; run it once as a reference point. Run 40x at N50 ≤ 30 kb if the
    parent supply allows; check it with the inventory.
  - **Regional arm (Q5):** 15x and 30x per sample, which is 7.5x and 15x per haplotype. These are the coverages
    regional cohort data typically have.
- **Coverage ceiling:** a 100 kb N50, and fixed lengths of 100 kb or more, cannot reach 30x from these data. As
  in the first post, the ultra-long end is only covered at 15–20x. State this as a limitation; do not stretch
  the data.
- **Replicates:** a second cut seed at 20 kb N50 for every arm, to measure the noise.

### Read selection: the main design choice

The pilot's parent reads were selected by aligning the full-length read to the pilot chromosomes and their
segmental duplication partners. For assembly this introduces a bias. A short piece of a long read that maps
uniquely is kept, even if the piece on its own would have mapped to a paralog elsewhere. Short reads therefore
get pre-cleaned input, which makes them look better in repeats.

Options:

| Design | Cost | Bias |
|---|---|---|
| A. Whole genome: cut all reads, assemble everything | ~1 day per assembly at 32 threads, 100–200 GB RAM (estimate, to measure) | none; numbers comparable to the literature |
| B. Chromosome subset, pieces selected by *their own* alignment (cut the whole data set, align the pieces, keep those overlapping chr7 + chr20 + SD partners) | one whole-data alignment per condition (hours at 32 threads), then ~1 h per assembly | small: mimics what a targeted analysis does |
| C. Chromosome subset from the existing parents | cheapest | short pieces pre-cleaned by long-read mapping (above) |

Plan:
- **B for the grid.**
- **A for three points** (10, 30 and 50 kb N50 at 30x), to anchor the absolute numbers.
- **C only for development and debugging.**

Chromosomes chr7 and chr20:
- chr20 is small and well behaved.
- chr7 has segmental duplication-rich and tandem repeat-rich regions (the hotspots seen in the SV post) and a
  large centromere. The Q100 truth covers them.
- Optionally add chr1 for its large segmental duplications.

### Regional arm (Q5) design

- **Windows:** 50–100 windows of 500 kb, in three groups:
  - random windows on the pilot chromosomes, inside the GIAB benchmark regions;
  - windows chosen for segmental duplication content (long, near-identical copies) and for tandem repeat hotspots;
  - a few windows with known large insertions in HG002.

  The selection is scripted and seeded.
- **Reads per window:** the pieces whose alignment overlaps the window plus 50 kb flanks, from the same aligned
  pieces as the SV study (one alignment per condition).
- **Haplotype assignment, two modes:**
  - *Oracle:* assign each piece to the Q100 haplotype it aligns to better; ambiguous pieces go to both. This
    isolates assembly from phasing.
  - *Realistic:* phase SNVs with LongPhase at that coverage, then haplotag the pieces; untagged pieces go to both
    haplotypes. Phasing error then shows up in the assembly.
- **Assemblers, per haplotype:**
  - Flye `--nano-hq` (what regional pipelines often use);
  - hifiasm `--ont` (unphased per haplotype, and diploid on all reads of the window).
- **Truth per window:** the Q100 haplotype sequence matching the window. Get it by mapping each Q100 haplotype to
  GRCh38 (minimap2 asm5; dipcall's alignment if it is kept) and lifting the window boundaries.
- **Measures per haplotype:**
  - complete (one contig spanning the window, ≥ 99% of the truth length);
  - the number of contigs;
  - the largest contig as a fraction of the window;
  - identity to the truth (edit distance per Mb, QV);
  - SVs against the truth haplotype (alignment, then calls, then truvari with the GIAB setup);
  - misjoins: contig segments mapping outside the window or to the wrong haplotype.

## Tools (pin versions at implementation)

| Step | Tool |
|---|---|
| Genome-scale assembly | hifiasm `--ont` (current release), outputs `*.bp.hap1/hap2.p_ctg.gfa`, primary and unitigs |
| Regional assembly | Flye 2.9.x `--nano-hq`; hifiasm `--ont` |
| GFA to FASTA, contig stats | `gfatools gfa2fa`, `seqkit stats -a` / gfastats (count, N50, NG50, auN, largest, total) |
| Correctness | QUAST (`--large`, against the matching Q100 haplotype): misassemblies, NGA50, duplication ratio; or minimap2 `-x asm5` + `paftools.js asmstat/misjoin` as a lighter cross-check |
| Base accuracy | yak QV against Q100 k-mers (`yak count` on the Q100 haplotypes); Merqury as an alternative |
| Phasing | yak `trioeval` with k-mers from the two Q100 haplotypes (stand-ins for the parents): switch and Hamming errors; phased block N50 from the contig assignments |
| Gene completeness (optional) | compleasm on the chromosome subset or the whole genome |
| Assembly SVs | dipcall (`-z200000,10000`, as GIAB) on hap1/hap2, then truvari with the GIAB setup + refine and the same stratification script |
| Read phasing (Q5, realistic) | LongPhase (SNVs), `longphase haplotag` |

New conda environments are needed for hifiasm, flye, quast, yak, gfatools/seqkit and dipcall.

## Workflow (in this repository)

- Add `workflow/rules/assembly.smk` (genome-scale) and `workflow/rules/regional.smk` (Q5), plus a config block
  `assembly:` with:
  - conditions (N50s, fixed lengths, coverages, seeds);
  - the design (A/B/C) and chromosomes;
  - the window selection;
  - the assemblers and their arguments.
- Reuse the existing pieces:
  - the cutting rules (`cut_reads.py`, both modes) and the inventory for the coverage checks;
  - the alignment rule (minimap2 `lr:hq`) for design B and the regional reads;
  - `stratify_truvari.py` for the assembly SVs, so the strata are identical to the SV post.
- New scripts:
  - `select_windows.py`;
  - `assign_haplotype.py` (oracle);
  - `window_truth.py`;
  - `assembly_stats.py`, which collects all measures into one `results/assembly/summary.tsv` and
    `results/regional/summary.tsv`;
  - `blog_figures_assembly.py`.
- Resources:
  - hifiasm whole genome: `threads: 32`, a declared `mem_mb` (set after the timing test);
  - the other jobs at ≤ 16 threads;
  - all large intermediates marked `temp()`, except the final assemblies, which are small.
- Targets: `assembly_regional`, `assembly_subset`, `assembly_wholegenome`, `assembly_all`.

## Compute and server budget (estimates, to be measured in milestone 1)

| Arm | Jobs | Cost |
|---|---|---|
| Q5 regional | ~75 windows × 2 haplotypes × ~10 conditions × 2 assemblers ≈ 3,000 small jobs (minutes each, 4 threads) | ~1–2 days at 16 cores |
| B alignments | ~12 conditions × whole-data alignment | ~12 × 3–6 h at 32 threads |
| B assemblies + evaluation | ~12–15 hifiasm runs on chr7+chr20 | ~1 h each at 16 threads, plus QUAST/dipcall |
| A whole genome | 3 hifiasm runs | ~1 day each at 32 threads, run one at a time; peak memory to check against the server first |

- **Core budget:** stay within the shared-server budget, together with the other runs. hifiasm jobs run one at
  a time, niced.
- **Memory:** check free memory before every whole-genome run (an earlier crash came from one memory-heavy
  process).
- **Disk:**
  - one full set of pieces at 30x is about 100 Gb of FASTQ (gzip about 50–60 GB), so keep pieces `temp()` and
    process the conditions in series;
  - assemblies are small (a few GB in total).

## Order and gating

0. **Precondition (now):** do not delete the pilot's parent reads (the bulk of the 88 GB of pilot
   intermediates) when the SV post is published. This plan reuses them. The rest of the pilot intermediates can
   go.
1. **Milestone 1, timing and plumbing (design C, development):**
   - chr20 only, two conditions (10 and 30 kb N50, 30x);
   - hifiasm + Flye + all evaluations end to end;
   - measure time and peak memory, then settle the grid.
2. **Milestone 2, Q5 regional arm:** the full window set, oracle and realistic phasing. This gives the first
   blog-worthy result, and the most relevant one for regional pangenomics.
3. **Milestone 3, Q1–Q4 on the chromosome subset (design B).**
4. **Milestone 4, three whole-genome anchors (design A).** Only if milestone 3 shows a pattern worth anchoring,
   and only when the server is quiet.
5. **Blog post(s):** one post, or two (regional, then genome) if the regional result stands on its own.

Start only after the current long runs on the server have finished, and with the user's go for each launch.

## Blog: core figures (draft)

1. Contig NG50 against library N50 at 30x (and 15x), with a dashed line for the Q100 truth.
2. Correctness: NGA50 and misassemblies per Gb against N50. Fragmented, or wrong?
3. Regional: the fraction of haplotypes assembled completely against N50, per window class (plain, SD, TR) and
   per phasing mode (oracle, realistic).
4. Assembly-based against read-based SV F1 at the same N50, for insertions over 1 kb and for SDs.
5. A compact table: N50 × coverage × (NG50, QV, switch error, complete haplotypes).

## Risks and open questions

- **hifiasm `--ont` with short reads:** its error correction may perform very poorly below about 10 kb N50.
  That is a result, but don't spend compute below 5 kb.
- **Assembling a chromosome subset** leaves out the reads of the rest of the genome. Report numbers relative to
  the matching whole-genome anchor, and never compare subset numbers with published genome assemblies.
- **Haplotype naming:** use Q100 mat/pat consistently. hifiasm hap1/hap2 are arbitrary, so assign them by
  alignment before scoring.
- **Oracle phasing is optimistic by design.** Present it next to the realistic mode, never alone.
- **Q5 assemblers:** Flye's per-haplotype output can contain haplotigs from mis-tagged reads. Count contigs
  mapping to the other haplotype as an error class.
- **Repository goes public:** no sample IDs, no server paths, and no mention of any in-house cohort. Use generic
  wording only (regional pangenome studies).
