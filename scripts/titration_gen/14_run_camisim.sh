#!/usr/bin/env bash
# Simulate one bacterial community sample (29-genome AWI-Gen-ranked panel,
# log-normal abundances) with CAMISIM 2.0's metagenomic Nextflow pipeline.
#
# CAMISIM itself is an external dependency, not vendored in this repo - see
# README.md "External dependencies" for the exact pinned version (tag,
# commit). Set CAMISIM_DIR and CONDA_PROFILE_SCRIPT to where CAMISIM and
# conda are installed; nextflow must be on PATH.
#
# Fixed design parameters (abundance distribution, size, seed, gsa/binning
# disabled) are applied additively via config/camisim_overrides.config (-c),
# so CAMISIM's own nested configs (art.config, conda.config,
# distribution.config) stay intact - see that file for what's overridden
# and why. The four paths below are passed as CLI flags instead, computed
# from this repo's own location, so the whole thing is portable to wherever
# this repo and CAMISIM are actually cloned.
#
# genomes_total/genomes_real and min/max_strains_per_otu are left at their
# CAMISIM defaults (2/2, 1/2): genomes_total==genomes_real disables
# CAMISIM's synthetic extra-strain generation, so all 29 genomes we
# explicitly listed are used directly and deterministically - see README
# for why that's the right call for a reproducible designed community.
#
# Run with ~4 CPUs and 8 GB RAM, e.g.:
#   CAMISIM_DIR=/path/to/CAMISIM bash scripts/titration_gen/14_run_camisim.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CAMISIM_DIR="${CAMISIM_DIR:?set CAMISIM_DIR to your CAMISIM 2.0.1 clone}"
CONDA_PROFILE_SCRIPT="${CONDA_PROFILE_SCRIPT:?set CONDA_PROFILE_SCRIPT to your conda.sh}"

source "$CONDA_PROFILE_SCRIPT"
command -v nextflow >/dev/null || { echo "nextflow not found on PATH" >&2; exit 1; }

cd "$CAMISIM_DIR"

nextflow run main.nf \
    --pipeline metagenomic \
    -c "${REPO_ROOT}/config/camisim_overrides.config" \
    --genome_locations_file "${REPO_ROOT}/camisim_inputs/genome_locations.tsv" \
    --metadata_file "${REPO_ROOT}/camisim_inputs/metadata.tsv" \
    --outdir "${REPO_ROOT}/camisim_out" \
    --ncbi_taxdump_file "${REPO_ROOT}/data/ncbi_taxdump/taxdump.tar.gz" \
    -resume
