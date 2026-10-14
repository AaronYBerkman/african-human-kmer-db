# Provenance: bacterial, human, and viral titration components

This is the single reference for exactly where every read in this project's titration benchmarks comes from — accession numbers, names, source database, and the publication backing each choice. It covers all components used across the HRR benchmarking effort: the bacterial genome panel, human reads (both the Global and Continental-African panels), and viral reads. The Global panel is called Forbes-27 or mixed-ancestry in file and script names, and the Continental-African panel is called african25.

Compiled by `scripts/titration_gen/19_compile_provenance_metadata.py`. Full per-row data in the sibling `.tsv` files in this directory.

---

## 1. Bacterial genomes — `bacterial_genomes.tsv` (29 rows)

**What:** One representative RefSeq genome assembly per bacterial species in the simulated gut community.

**Names/accessions:** Species binomial + RefSeq assembly accession (`GCF_...`) + NCBI taxid + species_taxid. See `bacterial_genomes.tsv` for all 29; full resolution detail (including the 25 AWI-Gen candidates that did *not* resolve, and why) is in `../results/species_resolution_log.tsv`.

**Source database:** NCBI RefSeq (`https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/`) — the actual genome sequence/assembly data.

**Selection rationale (why these 29 species specifically):**
> Maghini DG, Oduaran OH, Olubayo LAI, et al. Expanding the human gut microbiome atlas of Africa. *Nature* 2025;638:718–728. https://doi.org/10.1038/s41586-024-08485-8

The species panel is not a hand-picked "known gut taxa" list — it's resolved directly from this paper's Supplementary Data 3 (real measured species-abundance ranking across 6 AWI-Gen 2 African cohort study sites: Agincourt, Nanoro, DIMAMO, Navrongo, Nairobi, Soweto). See `../scripts/titration_gen/12_download_bacterial_genomes.py` for the deduplication and resolution method.

---

## 2. Human reads — `human_reads.tsv` (27 rows, 26 populations)

**What:** One representative ENA sequencing run per 1000 Genomes population, replicating Forbes et al.'s mixed-ancestry human background design.

**Names/accessions:** Population code + full population name + superpopulation (EUR/AFR/EAS/SAS/AMR) + ENA run accession (`ERR...`) + ENA study accession (`ERP...`). See `human_reads.tsv` for all 27.

**Source database:** ENA, under the 1000 Genomes Project 30× WGS umbrella study `PRJEB31736`. (Correction: the table's `ena_study_accession` column, `ERP114329`–`ERP114355`, does *not* actually resolve to per-population sub-studies as originally assumed — confirmed by directly querying ENA, where the real `study_accession` for these runs is just `PRJEB31736` itself, and one of the values, `ERP114331`, is an unrelated sea cucumber sequencing study on ENA. Kept as a descriptive field only; nothing in this pipeline relies on it.)

**Data publication (the underlying sequencing):**
> 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature* 2015;526:68–74. https://doi.org/10.1038/nature15393

**Design rationale (why these specific 27 runs / mixed-ancestry composition):**
> Forbes M, Ng DYK, Boggan RM, et al. Benchmarking of human read removal strategies for viral and microbial metagenomics. *Cell Reports Methods* 2025;5:101218. https://doi.org/10.1016/j.crmeth.2025.101218

**Cross-tier identity guarantee, verified:** the microbial × Forbes-27 arm (`results/titrations/`, built from this pool) and the viral × Forbes-27 arm (the ERP169281 accessions, used as-is) must share the same 27 human individuals, or any microbial-vs-viral difference observed later could be an artifact of different human backgrounds rather than the content type. `scripts/titration_gen/22_verify_forbes27_cross_tier_identity.py` checks this directly: every one of Forbes' own mixed ERP169281 samples carries its embedded human reads' ENA run accession in field 2 of the read header (same tagging convention used throughout this repo). Scanning all 27 of the viral benchmark's own source accessions (`reads/viral_raw/`, from `20_extract_viral_pool.py`) shows **every single one embeds exactly this same set of 27 accessions** — confirmed, not assumed (`results/forbes27_cross_tier_identity_check.tsv`). This isn't a coincidence: Forbes et al. evidently used one fixed 27-individual mixed-ancestry pool uniformly across their entire titration dataset, and the table in this row is that same pool.

---

## 3. Continental-African human reads, 5 individuals/site — `human_reads_african25.tsv` (25 rows, 5 populations)

**What:** Five independent individuals per continental-African (sub-Saharan) 1000 Genomes population (25 total) — the ablation arm restricting Forbes et al.'s 27-population background to the 5 populations actually sampled within Africa, with 5 individuals per site rather than 1 to reduce the risk of a single idiosyncratic genome skewing the African background.

**Names/accessions:** YRI (Nigeria), LWK (Kenya), MSL (Sierra Leone), ESN (Nigeria), GWD (Gambia) — 5 individuals each. Excludes ASW (US) and ACB (Barbados) — diaspora populations, not continental-Africa samples, despite their "AFR" 1000 Genomes superpopulation label. See `human_reads_african25.tsv` for all 25 rows (same columns as `human_reads.tsv`); the originating individual IDs (`NA.../HG...`) are in `../data/continental_african_5x5_accessions.tsv`.

**Resolution method:** joined the canonical [1000 Genomes/IGSR pedigree file](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped) (individual → population) against `PRJEB31736`'s own run list (`sample_title` consistently `Coriell GM/HG<digits>` → individual ID, confirmed for all 37,090 runs) via `../scripts/titration_gen/21_resolve_african_multi_individual_accessions.py`, which documents the method in full. This process also surfaced a correction: the `ena_study_accession` values in the original per-population table don't actually resolve to per-population sub-studies as assumed (see §2 above).

**Source database / publications:** same as §2 above (ENA / 1000 Genomes `PRJEB31736`; 1000 Genomes Consortium 2015; Forbes et al. 2025 for the underlying mixed-ancestry design this ablation restricts).

---

## 4. Viral reads — `viral_reads.tsv` (27 rows, *only the accessions actually used*)

**What:** The exact 27 ERP169281 accessions used by this repo's viral benchmark arms — 3 sample types × human-titre deciles 1-9 × iteration 1: `clinical_rsv_b`, `clinical_sars_cov_2_high_ct`, `h3n2_control`.

**This is a deliberate subset, not the full catalog:** the source titration/accession sheets actually catalog 351 accessions (13 sample types × 9 titre deciles × 3 iterations). `viral_reads.tsv` is filtered down to the 27 that are genuinely used anywhere in this project (`_ena_utils.BENCHMARK_VIRUS_TYPES` / `BENCHMARK_ITERATION`, the single source of truth both `19` and `20` filter against), so this sheet reflects actual usage rather than the broader universe of what's available.

**Source database:** ENA project `ERP169281` (the complete deposited dataset, already mixed with human reads by the original authors — `20_extract_viral_pool.py` strips that human background back out and remixes with our own human pool; see that script's docstring).

**Publication:**
> Forbes M, Ng DYK, Boggan RM, et al. Benchmarking of human read removal strategies for viral and microbial metagenomics. *Cell Reports Methods* 2025;5:101218. https://doi.org/10.1016/j.crmeth.2025.101218

---

## Summary table

| Component | # entries | Source database | Selection/design publication |
|---|---|---|---|
| Bacterial genomes | 29 | NCBI RefSeq | Maghini et al. 2025, *Nature* |
| Human reads (Forbes-27 arm) | 27 (26 populations) | ENA / 1000 Genomes (`PRJEB31736`) | 1000 Genomes Consortium 2015, *Nature*; Forbes et al. 2025, *Cell Reports Methods* — **verified identical to the viral × Forbes-27 arm's embedded human accessions** |
| Human reads (continental-African, 5/site) | 25 (5 populations) | ENA / 1000 Genomes (`PRJEB31736`) | same as above |
| Viral reads (used accessions only) | 27 | ENA (`ERP169281`) | Forbes et al. 2025, *Cell Reports Methods* |
