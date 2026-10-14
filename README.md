# An African-specific k-mer database for improved human read removal in viral and bacterial metagenomic data

Code, benchmark data and results for the paper of the same name (citation to
be added on publication).

The NCBI SRA Human Scrubber removes human reads from metagenomic data by
matching them against a k-mer database built from the human reference
(the "default" database below). Human sequence that is common in African
genomes but missing from that reference is not removed. This repository builds
a supplementary k-mer database from three African genome assemblies, with all
k-mers shared with bacteria, archaea or viruses removed, and benchmarks it
against the default on synthetic bacterial and viral metagenomes and on real
metagenomes from the AWI-Gen 2 cohort.

## Contents

- [Repository layout](#repository-layout)
- [Reproducing the figures and tables](#reproducing-the-figures-and-tables)
- [Using the database](#using-the-database)
- [Full pipeline](#full-pipeline)
- [Human background panels](#human-background-panels)
- [Tagging convention](#tagging-convention)
- [External dependencies](#external-dependencies)
- [Results](#results)
- [Limitations](#limitations)

## Repository layout

```
├── config.sh              Paths for the cluster steps (edit before running stages 1-3 and 5)
├── config/                CAMISIM settings for the simulated bacterial community
├── scripts/
│   ├── db_build/              1. Build the African k-mer database           (01-11)
│   ├── titration_gen/         2. Build the synthetic benchmark data          (12-22)
│   ├── benchmarks/            3. Classify reads with each database           (23-27)
│   ├── analysis/              4. Statistics and figures                      (28-31, 35)
│   └── real_data_validation/  5. AWI-Gen 2 real-data validation              (32-34)
├── data/                  Input tables: 1000 Genomes accessions, gut species ranking
├── metadata/              Provenance of every genome and read set (start with PROVENANCE.md)
├── results/               Manifests for the genome panel, read pools and titrations
├── benchmark_results/     Classifier output, one folder per benchmark arm
├── stats/                 Statistical tables
├── figures/               Figures (PDF and PNG)
└── logs/                  SLURM logs are written here (empty in the repository)
```

Scripts are numbered in the order they were run.

## Reproducing the figures and tables

Stage 4 needs only R and the files already in this repository. From the
repository root:

```bash
Rscript scripts/analysis/28_supplementary_table_S1.R
Rscript scripts/analysis/30_paper_figures.R
Rscript scripts/analysis/31_per_individual_analysis.R
Rscript scripts/analysis/35_figure_pipeline.R
Rscript scripts/analysis/36_display_tables_S5_S7.R
```

| Figure file | Script | Shows |
|---|---|---|
| `figures/pipeline_schematic` | 35 | How the database is built |
| `figures/sensitivity` | 30 | Sensitivity of each method across 10-90% human content |
| `figures/specificity` | 30 | Specificity of each method across 10-90% human content |
| `figures/per_individual_sensitivity` | 31 | Sensitivity per individual, Global panel |
| `figures/ancestry_gain` | 31 | Per-individual gain by ancestry group, Global panel |
| `figures/awigen_validation` | 34 | Extra human reads removed in 600 AWI-Gen 2 metagenomes |

| Table file | Script | Contents |
|---|---|---|
| `stats/supplementary_table_S1_full_stats.csv` | 28 | Sensitivity and specificity with Wilson 95% CIs and McNemar tests, every arm and titration |
| `stats/population_comparison_wilson_ci.csv` | 28 | Global minus Continental-African difference for each cell |
| `stats/supp_tab_per_individual_wilcoxon.csv` | 31 | Per-individual paired Wilcoxon tests, both panels |
| `stats/wilcoxon_per_individual_forbes27.csv` | 31 | The same tests, Global panel only |
| `stats/per_individual_summary_annotated.csv` | 31 | Per-individual sensitivities with ancestry annotation |
| `stats/ancestry_gain_summary.csv`, `stats/ancestry_mannwhitney.csv` | 31 | Gain by ancestry group within the Global panel |
| `stats/stats_display.csv` | 36 | Supplementary Table S5: sensitivity and specificity with 95% CIs, formatted for LaTeX |
| `stats/gap_display.csv` | 36 | Supplementary Table S7: sensitivity gain over the default with 95% CIs and McNemar p-values, formatted for LaTeX |
| `stats/awigen_validation_wilcoxon.csv` | 34 | AWI-Gen 2 validation, overall and by country |

Script 34 needs the AWI-Gen 2 participant-level files, which are not
distributed (see [stage 5](#5-awi-gen-2-real-data-validation-scriptsreal_data_validation));
its outputs are included.

## Using the database

The database is built from ethically protected human genome data and is not
distributed. It can be used through our web service,
<https://human-scrubber-dev.afrigen-d.dev/>, without needing access to the
underlying genomes.

In the benchmarks it was run two ways with the NCBI SRA Human Scrubber:
alone, and as a second stage on the reads left after the default database.
In the titration benchmarks the two gave sensitivities within 0.04 percentage
points of each other in every arm (see [Results](#results)).

## Full pipeline

Stages 1, 2, 3 and 5 ran on a SLURM cluster and need large downloads and
several hundred GB to a few TB of intermediate data. They are included so that
every step can be inspected and rerun. To rerun them, set the paths in
`config.sh`, run `mkdir -p logs`, and submit each script from the repository
root (for example `sbatch scripts/db_build/01_build_index.sh`). The `#SBATCH`
lines give the resources each step used.

### 1. Database construction (`scripts/db_build/`)

```
African assembly k-mers − bacterial and archaeal k-mers − viral k-mers = african_human_filtered.db
```

The African k-mers come from three diploid African genome assemblies (from
South Africa, Kenya and Burkina Faso; six haplotypes). Every genome is indexed
with the NCBI STAT tools as 32-mers, keeping every position in each 4-position
window (a dense index, not a MinHash sketch). The African set is then
subtracted against 60,294 bacterial, 707 archaeal and 14,626 viral complete
RefSeq genomes. The final database contains 1,012,561,615 k-mers.

| Scripts | Step |
|---|---|
| 01-02 | Index and merge the six African haplotypes |
| 03-06 | Download, decompress, index and merge bacterial and archaeal genomes |
| 07-10 | Download, decompress, index and merge viral genomes |
| 11 | Subtract the microbial, then the viral, k-mers |

Steps run by hand and not scripted here:

- Downloading NCBI RefSeq `assembly_summary.txt` for bacteria, archaea and
  viruses into `${WORK_DIR}/taxonomy/` (read by scripts 03 and 07).
- Writing the genome list files: `files.list` (the six African haplotype
  FASTAs; scripts 01-02), `microbial_files.list` (script 06), and the 16-way
  split lists read by the index arrays (scripts 05 and 09). The merge lists
  use the format `file size <tab> path` that script 10 writes for the viral
  genomes.

The African assemblies and the database built from them are ethically protected human data and are not distributed.

### 2. Benchmark data construction (`scripts/titration_gen/`)

| Scripts | Step |
|---|---|
| 12-14 | Resolve and download a 29-species gut bacterial panel ranked by abundance in the AWI-Gen 2 cohort (Maghini et al. 2025), then simulate a community with CAMISIM |
| 15-16 | Download 1000 Genomes human reads and build tagged read pools |
| 17-18 | Mix human and non-human reads at 10-90% human content and verify the ground truth |
| 19 | Compile the provenance tables in `metadata/` |
| 20 | Extract viral reads from Forbes et al. (2025) clinical samples |
| 21 | Resolve the 25 Continental-African individuals |
| 22 | Check that the viral and bacterial Global arms use the same 27 individuals |

Every genome, run accession and publication behind the benchmark data is
listed in [`metadata/PROVENANCE.md`](metadata/PROVENANCE.md). Script 19 also
reads two sample sheets from Forbes et al. (2025) that are not included;
its output, `metadata/viral_reads.tsv`, is.

### 3. Benchmark runs (`scripts/benchmarks/`)

Each titration file is classified three ways: default database, African
database alone, and the two-stage combination, using the Human Scrubber's
`aligns_to` program. A read can match more than one database k-mer, so hits
are deduplicated by read ID before counting true and false positives.

| Script | Arm |
|---|---|
| 23 | Bacterial community × Global panel |
| 24 | Bacterial community and viral × Continental-African panel |
| 25 | Viral × Global panel |
| 26 | Per-individual sensitivity, Global panel (30,000 reads per person) |
| 27 | Per-individual sensitivity, Continental-African panel (30,000 reads per person) |

### 4. Analysis (`scripts/analysis/`)

Wilson score confidence intervals, McNemar's test on paired read-level
calls, and, because reads from the same person are not independent, a paired
Wilcoxon signed-rank test with one value per individual. See
[Reproducing the figures and tables](#reproducing-the-figures-and-tables).

### 5. AWI-Gen 2 real-data validation (`scripts/real_data_validation/`)

Script 32 selects 600 individuals from the AWI-Gen 2 microbiome cohort (150
per country: South Africa, Kenya, Burkina Faso and Ghana), excluding the one
participant whose genome was also used to build the database and one with no
sequencing data. Script 33 counts
the reads each database removes from every individual's raw sequencing data,
and script 34 summarises the results. The cohort data, the list of selected
participants and the per-participant results are controlled-access and are not
distributed (they are listed in `.gitignore`); access is described in the
paper's Data availability statement.

## Human background panels

| Name in the paper | Name in file and script names | Individuals |
|---|---|---|
| Global | Forbes-27, mixed ancestry | 27 individuals from all 26 1000 Genomes populations (CEU twice), the background used by Forbes et al. (2025) |
| Continental-African | african25 | 25 individuals, five from each of YRI, ESN, GWD, MSL and LWK |

The Continental-African panel leaves out ASW and ACB, which share the 1000
Genomes AFR label but were sampled outside Africa. The bacterial background is
labelled "Microbial" in the figures and tables.

## Tagging convention

Every simulated read carries its source in field 2 of the FASTA header: an
ENA run accession (`ERR...`) for human reads, and a RefSeq genome accession
(`GCF_...`) or viral sample type for non-human reads. All ground-truth counts
are computed from this field.

## External dependencies

| Software | Version |
|---|---|
| [NCBI SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) | 2.2.1 (Singularity image) |
| [NCBI STAT / ngs-tools](https://github.com/ncbi/ngs-tools) | branch `tax` |
| [CAMISIM](https://github.com/CAMI-challenge/CAMISIM) | tag 2.0.1, commit `7ce6013c6d5a0fac8ba8a80e52a03560d3546fa6` |
| NCBI taxonomy (`taxdump.tar.gz`, for CAMISIM) | downloaded 25 June 2026, placed in `data/ncbi_taxdump/` |
| R | analysis verified with R 4.6.1: ggplot2 4.0.3, dplyr 1.2.1, readr 2.2.0, tidyr 1.3.2, purrr 1.2.2, patchwork 1.3.2, scales 1.4.0 (ggrepel optional) |
| Python | 3.9; standard library only |

## Results

**Per-individual sensitivity** (paired Wilcoxon signed-rank test, 30,000 reads
per individual):

| Panel | n | Median, default | Median, African DB | Median gain | Missed reads recovered | p |
|---|---|---|---|---|---|---|
| Global | 27 | 95.81% | 99.03% | +3.2 points | 77.7% | 5.9×10⁻⁶ |
| Continental-African | 25 | 90.88% | 98.33% | +7.4 points | 81.1% | 1.3×10⁻⁵ |

"Missed reads recovered" is the share of the human reads missed by the
default that the African database removes, pooled across individuals. Every
individual in both panels gained (smallest gains 2.4 and 2.9 points), so each
Wilcoxon statistic is at its maximum. The gain on the Continental-African
panel is more than twice that on the Global panel.

**Specificity**, summed over every background, panel and titration:

| Method | False positives | Non-human reads | Specificity |
|---|---|---|---|
| Default | 2 | 3,600,000 | 99.99994% |
| African DB only | 20 | 3,600,000 | 99.99944% |
| Combined (two-stage) | 20 | 3,600,000 | 99.99944% |

## Limitations

The Continental-African panel comes from 1000 Genomes populations in West
Africa and Kenya. Southern African populations (for example San and Khoe) and
Central African rainforest populations, which are likely to differ most from
existing reference sequence, are not represented, so the gains reported here
may not hold for them.

## Citation

Citation to be added on publication.

## Licence

MIT; see [LICENSE](LICENSE).

## Author

Aaron Y. Berkman
Division of Human Genetics, Sydney Brenner Institute for Molecular Bioscience,
University of the Witwatersrand, Johannesburg, South Africa

Supervisors: Prof. Scott Hazelhurst, Dr. Luicer Olubayo

## References

1. Katz KS, Shutov O, Lapoint R, Kimelman M, Brister JR, O'Sullivan C. STAT: a fast, scalable, MinHash-based k-mer tool to assess Sequence Read Archive next-generation sequence submissions. *Genome Biology* 2021;22:270.
2. Forbes M, Ng DYK, Boggan RM, et al. Benchmarking of human read removal strategies for viral and microbial metagenomics. *Cell Reports Methods* 2025;5:101218.
3. 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature* 2015;526:68–74.
4. Maghini DG, Oduaran OH, Olubayo LAI, et al. Expanding the human gut microbiome atlas of Africa. *Nature* 2025;638:718–728.
