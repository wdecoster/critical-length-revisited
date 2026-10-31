# Follow-up plan: de novo assembly vs read length

Status: plan agreed in outline (2026-10-31, decisions below). Waits until the server has room; see "Order and
gating".

## Decisions (2026-10-31)

1. **The regional arm (Q5) comes first, and is the focus for a long time.** The genome-scale arms follow later.
2. **Realistic phasing** is the main mode of the regional arm: what a regional pipeline does in practice
   (below). The oracle mode stays only as a ceiling.
3. **Named locus:** MAPT (17q21.31), as the hard case study (inversion polymorphism, large near-identical
   segmental duplications). The statistical backbone is the set of window classes; see "Regional arm" for the
   choice of an easy reference locus.
4. **Read selection: design B** (pieces selected by their own alignment).
5. **Assemblers:** Flye and hifiasm for the regional arm, hifiasm for the genome scale. Verkko and hifiasm's
   ultra-long mode are out of scope.
6. **Coverage:** 20x and 30x per sample, for both arms. 15x per sample is too low for reliable regional assembly
   in practice. 20x and 30x also show when read length beats coverage and when it does not.
7. **Data:** the pilot data stay (no deletions).

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
  - **Both arms:** 20x and 30x per sample (10x and 15x per haplotype). The 20x runs are subsampled from the 30x
    alignment, as in the pilot. hifiasm `--ont` relies on read overlaps for error correction, so lower
    coverages are not tested.
- **Coverage ceiling:** 30x reaches an N50 of 50 kb; 20x also reaches an N50 of 100 kb (about 28x of the data is
  in reads of 100 kb or more). So "20x at 100 kb against 30x at 30 kb" is a direct test of read length against
  coverage. Fixed lengths of 150 kb only reach about 16x; leave them out.
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
- **B for the grid** (decided).
- **A for three points** (10, 30 and 50 kb N50 at 30x), to anchor the absolute numbers.
- **C only for development and debugging.**

Chromosomes chr7 and chr20:
- chr20 is small and well behaved.
- chr7 has segmental duplication-rich and tandem repeat-rich regions (the hotspots seen in the SV post) and a
  large centromere. The Q100 truth covers them.
- Optionally add chr1 for its large segmental duplications.

### Regional arm (Q5) design

- **Windows:** 50–100 windows of about 450–500 kb, anywhere in the genome (design B aligns everything), in
  classes:
  - *plain:* random windows inside the GIAB benchmark regions, with few repeats. These are the easy loci of a
    typical regional study, and the reference for everything else;
  - *segmental duplications:* windows with long, near-identical copies;
  - *tandem repeat hotspots;*
  - *large insertions:* windows with a known large insertion in HG002;
  - *MAPT:* the whole 17q21.31 region as a named case study (about 1–1.5 Mb, coordinates set from the segmental
    duplication annotation and the Q100 alignment; check HG002's H1/H2 status first).

  The selection is scripted and seeded. For an easy named example in the post, pick one plain window after the
  fact; no separate named simple locus.
- **Reads per window:** the pieces whose own alignment overlaps the window plus 50 kb flanks. Cheaper equivalent
  of design B for this arm: select the full-length reads aligned to the window ± 1 Mb *and* to its segmental
  duplication partners elsewhere, cut them, align the pieces to the whole genome, and keep the pieces that land
  in the window. The only difference from full B: pieces from unannotated similar sequence elsewhere that would
  mismap into the window are missed. Decide at implementation whether to use full B (shared with the genome
  arm) or this partner version.
- **Haplotype assignment (realistic, main mode):** what a regional pipeline does:
  - SNV calling (DeepVariant for R10) and SV calling (Sniffles2) in the window plus flanks;
  - phasing with LongPhase `--ont`, with the SVs;
  - haplotagging the pieces;
  - untagged pieces (homozygous stretches) go to both haplotypes.

  A regional pipeline across many samples can join phase blocks with population phasing. A single sample cannot,
  so block switches remain. Optionally test joining blocks with a reference panel (SHAPEIT5 with the 1000 Genomes
  panel; HG002 is not in it).
- **Oracle mode (ceiling only):** assign each piece to the Q100 haplotype it aligns to better; ambiguous pieces
  go to both.
- **Assemblers, per haplotype:**
  - Flye `--nano-hq --no-alt-contigs` (what regional pipelines often use);
  - hifiasm `--ont` in haploid mode on the haplotype's reads, and diploid on all reads of the window.
- **Truth per window:** the Q100 haplotype sequence matching the window. Get it by mapping each Q100 haplotype to
  GRCh38 (minimap2 asm5; dipcall's alignment if it is kept) and lifting the window boundaries.
- **Measures per haplotype:**
  - complete: the contigs cover ≥ 95% of the window, and one contig spans it;
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
| Read phasing (Q5, realistic) | DeepVariant (R10 model) + Sniffles2 in the window; LongPhase `--ont` with SVs; `longphase haplotag` |

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

0. **Data:** keep the pilot data, including the parent reads, and the full ultra-long read set.
1. **Milestone 1, regional plumbing and timing:**
   - 5 windows (2 plain, 1 segmental duplication, 1 tandem repeat, MAPT), 2 conditions (10 and 30 kb N50, 30x);
   - the realistic phasing chain, Flye and hifiasm, and all regional measures, end to end;
   - measure time, peak memory and disk; settle the read selection (full B or the partner version) and the grid.
2. **Milestone 2, the regional arm:** the full window set, N50 grid and fixed lengths at 20x and 30x,
   realistic phasing (oracle as the ceiling). The first blog post.
3. **Milestone 3, Q1–Q4 on the chromosome subset (design B), at 20x and 30x.** Later.
4. **Milestone 4, three whole-genome anchors (design A).** Only if milestone 3 shows a pattern worth anchoring,
   and only when the server is quiet.
5. **Blog posts:** regional first; genome scale as a separate post later.

Start only after the current long runs on the server have finished, and with the user's go for each launch.

## Blog: core figures (draft)

Regional post:
1. The fraction of haplotypes assembled completely against library N50, at 20x and 30x, per window class.
2. Read length against coverage: 20x at a long N50 next to 30x at a shorter one.
3. Realistic phasing against the oracle ceiling: how much is lost to phasing.
4. MAPT as the worked example: what each N50 gets you there.
5. A compact table: N50 × coverage × window class (complete haplotypes, contigs, identity).

Genome-scale post (later):
1. Contig NG50 against library N50 at 20x and 30x, with a dashed line for the Q100 truth.
2. Correctness: NGA50 and misassemblies per Gb against N50. Fragmented, or wrong?
3. Assembly-based against read-based SV F1 at the same N50, for insertions over 1 kb and for SDs.

## Risks and open questions

- **hifiasm `--ont` with short reads:** its error correction may perform very poorly below about 10 kb N50.
  That is a result, but don't spend compute below 5 kb.
- **Assembling a chromosome subset** leaves out the reads of the rest of the genome. Report numbers relative to
  the matching whole-genome anchor, and never compare subset numbers with published genome assemblies.
- **MAPT may not assemble at any N50.** Then that is the result: report how far each N50 gets (largest contig,
  which copies collapse) rather than only complete or incomplete.
- **Haplotype naming:** use Q100 mat/pat consistently. hifiasm hap1/hap2 are arbitrary, so assign them by
  alignment before scoring.
- **Oracle phasing is optimistic by design.** Present it next to the realistic mode, never alone.
- **Q5 assemblers:** Flye's per-haplotype output can contain haplotigs from mis-tagged reads. Count contigs
  mapping to the other haplotype as an error class.
- **Repository goes public:** no sample IDs, no server paths, and no mention of any in-house cohort. Use generic
  wording only (regional pangenome studies).
