# af_fish

Analysis code for:

> Fox, J. A., Hunt, D. A. G. A., Hendry, A. P., Chapman, L. J., & Barrett, R. D. H. (2024).
> Counter-gradient variation in gene expression between fish populations facilitates
> colonization of low-dissolved oxygen environments. *Molecular Ecology*, 33, e17419.
> https://doi.org/10.1111/mec.17419

**Raw data:** NCBI SRA, BioProject [PRJNA1041828](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1041828)

**Full GO results and gene-level tables** (counts, log2FC, FDR-corrected *p*-values) are in the
paper's supplementary materials.

---

## What this study did

We compared gill gene expression between low-dissolved-oxygen (low-DO, swamp) and high-DO
(stream) populations of two closely related cyprinid fishes in the Rwembaita Swamp System,
Kibale National Park, Uganda:

- ***Enteromius apleurogramma*** — a recent range expander into the system (first recorded 2015)
- ***Enteromius neumayeri*** — a long-established native species in the same drainage

Fish from each population were either sampled immediately on capture or held for 14 days in
aerated outdoor ponds before sampling. Crossing *population* with *sampling time* separates
candidate-evolved differences from plastic ones. Because there is no reference genome for either
species, the pipeline builds a de novo transcriptome per species with Trinity, quantifies against
it, and then tests for differential expression and for signatures of selection in SNPs called
from the same reads.

## Naming key

The code predates the current genus name and uses short codes throughout. **This table is the
key to reading every filename, variable and contrast in this repository.**

| Code in scripts | Paper | Meaning |
| --- | --- | --- |
| `BA`, `EA` | *E. apleurogramma* | Recent range-expanding species |
| `BN`, `EN` | *E. neumayeri* | Long-established native species |
| `stream`, `St` | high-DO (H) | Normoxic tributary stream population (~6 mg O₂/L) |
| `swamp`, `Sw` | low-DO (L) | Hypoxic papyrus swamp population (~1 mg O₂/L) |
| `field`, `F` | immediate (I) | Gills taken within ~10 min of capture |
| `pond`, `P` | acclimation (A) | Gills taken after 14 days in aerated holding ponds |
| `bf` | — | "before removal" — reads *before* the unused overrepresented-sequence filter |

So the four sample groups per species are:

| Group in filenames | Paper label |
| --- | --- |
| `stream_field` | H-I |
| `stream_pond` | H-A |
| `swamp_field` | L-I |
| `swamp_pond` | L-A |

### Contrasts

The paper uses three a priori planned contrasts (Methods §2.7, Figure 1):

| Contrast in filenames | Paper | Interprets as |
| --- | --- | --- |
| `swamp_field_vs_swamp_pond` | L-I vs L-A | Plasticity **+** laboratory effects |
| `stream_field_vs_stream_pond` | H-I vs H-A | Laboratory effects alone |
| `stream_pond_vs_swamp_pond` | L-A vs H-A | **Candidate-evolved** differences |

A fourth contrast, `stream_field_vs_swamp_field` (H-I vs L-I), is the raw field divergence. It is
not one of the three planned contrasts but is computed and used — for the heatmap clusters of
interest (§2.8) and for GO enrichment on field-divergent genes.

Gene sets built from these in `22_compareEvP.R` (§2.9):

- **evolved** = DEGs in L-A vs H-A
- **plastic** = DEGs in L-I vs L-A that are *not* DEGs in H-I vs H-A (i.e. plasticity with
  laboratory effects subtracted)
- **shared** = the intersection of the two; log2FC in these is correlated to test for
  counter-gradient variation

## Sample sizes

70 of 80 extractions passed QC and were sequenced; one *E. neumayeri* (low-DO, immediate) sample
was later dropped for poor quality, leaving 69 for analysis.

| | low-DO | high-DO |
| --- | --- | --- |
| Immediate | EN 8, EA 10 | EN 9, EA 9 |
| Acclimation | EN 7, EA 11 | EN 8, EA 7 |

The Trinity assemblies were built from all pre-QC samples per species (37 EA, 33 EN); see
`scripts/07_snpAnalysis/BA_files.txt` and `BN_files.txt` for the sample ID lists.

---

## Repository layout

```
scripts/
├── 01_preProcessing/    read QC and cleaning             (paper §2.5)
├── 02_trinityAssembly/  de novo transcriptome assembly   (§2.6)
├── 03_trinityQC/        assembly quality assessment      (§2.6)
├── 04_transcriptQuant/  quantification and edgeR DEGs    (§2.7)
├── 05_annotation/       ORF prediction and annotation    (§2.6)
├── 06_diffAnalysis/     downstream stats, clustering, GO (§2.8, 2.9, 2.11)
└── 07_snpAnalysis/      SNP calling, Fst, outlier scan   (§2.10, 2.11)
```

Scripts run in numeric order (`01` → `31`). Note that the file numbers run **across**
directories, not within them: `12` and `13` (in `04_transcriptQuant/`) must run before `14`
(in `03_trinityQC/`), because Ex90N50 needs the expression matrix.

### 01_preProcessing — read quality control

| Script | Does |
| --- | --- |
| `01_rcorrector.sh` | Rcorrector v1.0.4, default settings — correct erroneous k-mers |
| `02_filterUncorrectableJob.sh` + `02a_filterUncorrectable.py` | Drop read pairs flagged uncorrectable |
| `03_trimgaloreJob.sh` | TrimGalore! v0.6.6 — trim phred < 5 and adapters, min length 36 |
| `04_rRNAalignment.sh` | bowtie2 v2.3.4.3 (`--nofw --very-sensitive-local`) against a SILVA rRNA database (SSUParc + LSUParc v138.1); keeps pairs where **neither** read mapped |
| `05_fastqc.sh` | FastQC reports |
| `06_removeFastqcOverrepJob.sh` + `06a_removeFastqcOverrep.py` | **Not used in the published analysis.** Written to strip FastQC-flagged overrepresented sequences, but the paper's results use the reads from step 04 — hence the `bf` ("before removal") tag on every downstream file. |

The SILVA database is not included here; download SSUParc and LSUParc v138.1 from
https://www.arb-silva.de (the paper used files accessed December 2022) and build a bowtie2 index.

### 02_trinityAssembly

`07a_trinityAssembly_BA.sh`, `07b_trinityAssembly_BN.sh` — Trinity v2.15.0 via Singularity,
`--SS_lib_type RF` (dUTP libraries), one assembly per species. The sample FASTQ paths are listed
explicitly in each script.

Each script also retains an earlier, commented-out version of the same command that bound the read
directory into the container as `/data`; the live command uses absolute scratch paths instead.

### 03_trinityQC — assembly assessment

| Script | Does |
| --- | --- |
| `08_bowtieBuild.sh` | bowtie2 index + realign reads to their own assembly for RNA-seq read representation |
| `09_NxExN50.sh` | `TrinityStats.pl` — contig Nx statistics, transcript/gene counts, median contig length |
| `10_busco.sh` | BUSCO v5.2.2, `vertebrata_odb10`, transcriptome mode → paper Table 2 |
| `11_rnaQUAST.sh` | rnaQUAST/QUAST v5.0.2. Run during development; **not reported in the paper.** |
| `14_Ex90N50.sh` | Ex90N50 and TPM filtering thresholds. **Run after quantification** (steps 12–13). |

BUSCO is run `--offline`, so download `vertebrata_odb10` first (the paper used the set accessed
February 2023).

### 04_transcriptQuant — quantification and differential expression

| Script | Does |
| --- | --- |
| `12_salmon_alignment.sh` | Salmon via Trinity's `align_and_estimate_abundance.pl` (the paper reports Salmon v1.10.1, the version bundled in the Trinity v2.15.0 container) |
| `13_expMatrix.sh` | `abundance_estimates_to_matrix.pl` — gene and isoform count / TMM-EXPR matrices using the Trinity gene-transcript map |
| `15_qualityCheckReplicates.sh` | `PtR` — replicate comparison, sample correlation matrix, PCA |
| `16_DEGanalysis.sh` | `run_DE_analysis.pl --method edgeR` with `--min_reps_min_cpm '7,1'` (CPM ≥ 1 in ≥ 7 individuals — 7 being the smallest sampling group) |
| `17_clusterDEG.sh` | `analyze_diff_expr.pl -P 0.01 -C 2` — FDR < 0.01 and ≥ 4-fold change, the paper's DEG criteria |
| `EXTRA_filter_rerun.sh` | Exploratory rerun after filtering low-expression transcripts. **Not part of the published analysis**; kept for transparency. Almost entirely commented out. |

### 05_annotation

| Script | Does |
| --- | --- |
| `18_transDecoder.sh` | TransDecoder v5.7.0 — `LongOrfs` then `Predict` |
| `19_prepTrinotate.sh` | Trinotate v4.0.0 — Swiss-Prot and PFam (accessed May 2023), Infernal v1.1.4 vs Rfam v14.9, SignalP v6.0, tmHMM v2, eggNOG-mapper v2 |

`19_prepTrinotate.sh` is written as a sequence of commented blocks meant to be uncommented and run
one at a time (database creation → initialization → sequence analyses → report).

### 06_diffAnalysis — R analyses

Run from the repository root (open `af_fish.Rproj`); all paths are relative to it.

| Script | Does |
| --- | --- |
| `20_trinotateStats.R` | Annotation summary — unique Swiss-Prot hits, and hits at ≥ 80% identity |
| `21_statAnalysis.R` | PCA on all genes and on DEGs; hierarchical clustering (Euclidean distance, Ward linkage, `cluster` v2.1.4); Mfuzz v2.60.0 soft clustering (c = 9 for EN, c = 7 for EA, membership ≥ 0.7); `pheatmap` v1.0.12 heatmaps; holding-pond effect checks (Figure S1) |
| `22_compareEvP.R` | Builds the evolved / plastic / shared gene sets; Pearson correlation of log2FC in shared genes against a 10,000-iteration permutation null; chi-squared tests with Yates correction |
| `23_plotREVIGO.R` | Redraws REVIGO semantic-space plots. REVIGO is a web tool (http://revigo.irb.hr), so its output is pasted inline in this script as data — this is intentional, not stray output. |
| `24_goSeq.R` | goseq GO enrichment (length-bias corrected, Wallenius approximation) on the Mfuzz clusters of interest, background restricted to genes retained after expression filtering; `qvalue` v2.32.0, FDR < 0.05 |
| `25_ManU.R` | Mann–Whitney U tests on the magnitude of log2FC between species, plus `rstatix` effect sizes. Written as a standalone script to run on the cluster when the analysis would not fit in memory locally, so it reads `./plastManU.csv` and `./evolManU.csv` from the working directory rather than from `data/`. |
| `26_subset.R` | Exploratory subsetting / percentage-overlap checks. **Not part of the published analysis.** |

### 07_snpAnalysis

| Script | Does |
| --- | --- |
| `27_constructSuperTranscript.sh` | `Trinity_gene_splice_modeler.py` — collapse each assembly into SuperTranscripts, giving a genome-like reference for variant calling |
| `28_callSNPs.sh` | Trinity's `run_variant_calling.py` — STAR v2.7.11b alignment, Picard v2.26.3 duplicate marking |
| `29_filter_combine.sh`, `29_filter_combine_BN.sh` | bcftools v1.16 — sort, index and merge per-sample VCFs (the `_BN` variant differs only in requested resources and modules) |
| `30_freebayes.sh` | FreeBayes v1.3.6 on the merged BAM, then bcftools filtering: QUAL > 30, biallelic SNPs, MAF ≥ 5%, no missing data, HWE *q* < 0.05, and LD pruning in 10 kb windows keeping one SNP per pair with *r*² > 0.8 |
| `31_snp_analysis.R` | `hierfstat` Weir & Cockerham *F*<sub>ST</sub> between DO populations; `pcadapt` v4.4.0 outlier scan (5 PCs retained for EA, 8 for EN); `topGO` v2.56.0 GO enrichment on genes containing outlier SNPs (Fisher's exact test, `weight` algorithm, no multiple-testing correction per the algorithm's guidance) |
| `31_bayescan.sh` | BayeScan v2.1 outlier scan. **Not used in the paper** — the published outlier analysis is pcadapt. Kept for transparency. |
| `BA_files.txt`, `BN_files.txt` | Per-sample read-pair paths and output prefixes consumed by `28_callSNPs.sh` |
| `var_calling.py` | A copy of Trinity's own `run_variant_calling.py`, vendored for reference. It expects Trinity's `PyLib` on the path and will not run standalone from this repository. |

`31_snp_analysis.R` also contains commented-out OutFLANK and PANTHER sections from earlier
exploration. `library(OutFLANK)` is commented out to match, since OutFLANK installs from GitHub
rather than CRAN; uncomment it if you want to run those sections.

---

## Data layout for the R scripts

The small metadata files that define the experimental design **are** included under `data/`:

```
data/
├── samples_BA_sal_bf.txt     group assignment per individual, read by scripts 12, 15, 16, 17, 28
├── samples_BN_sal_bf.txt
├── pondData/                 pond_data_BA.csv, pond_data_BN.csv — holding-pond assignment per
│                             individual (Field / Pond1 / Pond2), read by 21_statAnalysis.R
└── raw/
    ├── collections.xlsx      master field record for all 76 fish: ID, species, site, total and
    │                         standard length, pond, sample type, whether the fish was sequenced
    ├── pond_data.xlsx        source spreadsheet for the pondData CSVs
    └── PondTempDOData.ods    daily holding-pond temperature and dissolved oxygen readings,
                              14 dates from 2017-06-06 — the source data for Table S1
```

> **On the pond readings.** `PondTempDOData.ods` logs four ponds, but only ponds 1 and 2 held fish
> used in this study — pond 1 held 20 individuals and pond 2 held 16 (36 acclimated fish in total,
> of which 33 were sequenced). Ponds 3 and 4 were monitored alongside them but contributed no
> individuals here, which is why Methods §2.3 describes two ponds. The per-individual pond
> assignment used in the analysis is `data/pondData/pond_data_*.csv`.

Everything else the R scripts read is **not** here — those are outputs of steps 01–19 and 27–30,
and the gene-level results tables are published as supplementary material with the paper. The
scripts expect them in this tree under the repository root:

```
data/
├── DEG/
│   ├── BN/salmon/     diffExpr.P0.01_C2.matrix
│   │                  BN_bf_new_sal.gene.counts.matrix.<contrast>.edgeR.DE_results
│   │                  BN_bf_new_sal.gene.counts.matrix.<contrast>.edgeR.DE_results.P0.01_C2.DE.subset
│   ├── BA/salmon/     (same, BA_ prefix)
│   ├── filtered/      BN/, BA/ — outputs of EXTRA_filter_rerun.sh (exploratory)
│   └── BN_evolved.tsv, BN_plastic.tsv, BA_evolved.txt, BA_plastic.tsv, *_list.tsv
├── totalExpr/         BN_bf_new_sal.gene.TMM.EXPR.matrix, .gene.counts.matrix (and BA_)
├── trinotate/         BN_uniprot_sprot.diamond.blastx.outfmt6 (and BA_), Trinotate reports
├── goSeqData/         <species>.Trinity.gene_lengths.txt, per-cluster and per-contrast
│                      up/down gene factor lists
└── snp/               pop_EA.txt, pop_EN.txt (individual → population; derivable from
                       data/raw/collections.xlsx),
                       BA.freebayes.results.filtered.hwe.ld.norm.edited.vcf (and BN)
```

`<contrast>` is one of `stream_field_vs_stream_pond`, `stream_field_vs_swamp_field`,
`stream_pond_vs_swamp_pond`, `swamp_field_vs_swamp_pond`.

Figures and tables are written to the working directory (repository root), not to an output
folder.

### Sample-to-group definitions

Scripts `12`, `15`, `16`, `17` and `28` take a Trinity samples file — `data/samples_BA_sal_bf.txt`
and `data/samples_BN_sal_bf.txt` — which assigns each individual to one of the four groups and is
therefore what defines the experimental design for edgeR. Each is tab-delimited, one line per
individual:

```
<group>	<group>_rep<n>_sal	<left.fq>	<right.fq>
```

where `<group>` is `stream_field`, `stream_pond`, `swamp_field` or `swamp_pond`. The read paths are
the cluster's container-internal `/data/...` paths and will need editing for your own system.

> **One exclusion.** `samples_BN_sal_bf.txt` lists 33 individuals, but the published
> *E. neumayeri* analysis used 32. Fish **0616g** (`swamp_field_BN_rep7_sal`) was dropped for poor
> sample quality after this file was written — this is the removed low-DO immediate sample noted in
> Methods §2.5. It is absent from `data/pondData/pond_data_BN.csv`, which has the final 32. Exclude
> it to reproduce the paper's group sizes (Table 1). `samples_BA_sal_bf.txt` needs no such
> adjustment: all 37 *E. apleurogramma* individuals were retained.

The same group assignment underlies `data/snp/pop_EA.txt` and `pop_EN.txt`, and the whole design can
be rebuilt from `data/raw/collections.xlsx` if needed.

---

## Software

**Compute environment.** Steps 01–19 and 27–31 are SLURM batch scripts written for the Digital
Research Alliance of Canada (formerly Compute Canada) clusters. They carry `#SBATCH` headers with
a specific allocation (`def-barrett`), `module load` calls for that software stack, absolute
scratch paths under `/home/janayfox/scratch/` or `/lustre04/scratch/janayfox/`, and references to
Singularity/Apptainer images (`trinityrnaseq.v2.15.0.simg`, `trinotate.v4.0.0.simg`) that are not
included here. **To reuse these scripts you will need to edit the account, module names and all
paths for your own system.** Two spellings of the project directory appear (`afFishRNA` and
`af_fish_RNA`) because it was renamed partway through the project.

Many scripts hold several alternative or successive commands with all but the active one commented
out — that was the working style during the project. Read each script before running it rather
than assuming every line executes.

**Command-line tools:** Rcorrector 1.0.4, TrimGalore! 0.6.6, FastQC, bowtie2 2.3.4.3, Trinity
2.15.0 (and 2.15.1 for SuperTranscripts), BUSCO 5.2.2, QUAST 5.0.2, Salmon 1.10.1, TransDecoder
5.7.0, Trinotate 4.0.0, Infernal 1.1.4, SignalP 6.0, tmHMM 2, eggNOG-mapper 2, STAR 2.7.11b,
Picard 2.26.3, FreeBayes 1.3.6, BCFtools 1.16, BayeScan 2.1 (exploratory only).

**Reference databases:** SILVA SSUParc + LSUParc 138.1, BUSCO `vertebrata_odb10`, Swiss-Prot,
PFam, Rfam 14.9. All are external downloads and are not stored in this repository.

**R 4.2.2** with: edgeR 3.40.2, Mfuzz 2.60.0, cluster 2.1.4, pheatmap 1.0.12, goseq, qvalue
2.32.0, topGO 2.56.0, pcadapt 4.4.0, hierfstat, vcfR, radiator, mmod, rstatix, plus tidyverse,
ggplot2, ggpubr, patchwork, viridis, RColorBrewer, ggrepel, scales, data.table, reshape2, tibble,
stringr, dendextend, factoextra, GGally, ggupset, VennDiagram, GO.db. `OutFLANK` is loaded in
`31_snp_analysis.R` but used only by commented-out code.

---

## License

The code in this repository is released under the [MIT License](LICENSE). The paper itself is open
access under CC BY-NC-ND; the sequence data are available from SRA under NCBI's terms. Please cite
the paper if you reuse this code.

---

## Contact

Janay A. Fox — janay.fox@mail.mcgill.ca — [ORCID 0000-0002-7895-2357](https://orcid.org/0000-0002-7895-2357)
