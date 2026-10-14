#!/usr/bin/env python3
"""
Select individuals from the AWI-Gen 2 microbiome cohort for the real-data
validation arm: does the African-specific database catch more human reads
than the default on genuine African metagenomic sequencing data, not just
the synthetic titration benchmark?

Input: the cohort metadata for the individuals included in the AWI-Gen 2
microbiome publication (Maghini et al. 2025), with columns metaG_id
(sequencing sample ID), study_id and site. These data are available under
controlled access and are not distributed with this repository.

Exclusions are read from a separate file (--exclusions; columns study_id,
reason) that is likewise not distributed, because it names cohort
participants. Two individuals were excluded:
  - the one cohort participant whose genome assembly was also used to build
    the African database, so that the database is never tested on a sample
    it was partly built from;
  - one participant with metadata but no sequencing files.

Selection: 6 sites map to 4 countries. South Africa has 3 sites (Agincourt,
Dimamo, Soweto); Kenya, Burkina Faso, and Ghana have 1 each. 150 individuals
per country, with South Africa split evenly (not proportionally) across its
3 sites - 50 each - so no single South African site dominates the country's
representation the way Agincourt's much larger cohort (534 of 959) would if
sampled proportionally.

Output: results/awigen_validation_samples.tsv (600 rows: 50+50+50+150+150+150;
not distributed, as it lists cohort sample IDs)
"""
import argparse
import csv
import random
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_OUTPUT_PATH = REPO_ROOT / "results" / "awigen_validation_samples.tsv"

BASE_SEED = 20260707  # date this arm was designed, for a stable/citable seed

# site -> (country, n_to_select)
SITE_PLAN = {
    "Agincourt": ("South Africa", 50),
    "Dimamo":    ("South Africa", 50),
    "Soweto":    ("South Africa", 50),
    "Nairobi":   ("Kenya", 150),
    "Nanoro":    ("Burkina Faso", 150),
    "Navrongo":  ("Ghana", 150),
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata-path", type=Path, required=True,
                   help="AWI-Gen 2 microbiome cohort metadata TSV (columns metaG_id, study_id, site)")
    p.add_argument("--exclusions", type=Path, required=True,
                   help="TSV with columns study_id, reason: participants to exclude")
    p.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH,
                   help="Output manifest TSV path")
    p.add_argument("--seed", type=int, default=BASE_SEED)
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.exclusions, newline="") as f:
        exclusion_reasons = {r["study_id"]: r["reason"] for r in csv.DictReader(f, delimiter="\t")}
    excluded_study_ids = set(exclusion_reasons)

    with open(args.metadata_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    excluded = [r for r in rows if r["study_id"] in excluded_study_ids]
    for r in excluded:
        reason = exclusion_reasons[r["study_id"]]
        print(f"Excluding {r['metaG_id']} (study_id={r['study_id']}, site={r['site']}) - {reason}.")

    rows = [r for r in rows if r["study_id"] not in excluded_study_ids]

    by_site = {}
    for r in rows:
        by_site.setdefault(r["site"], []).append(r)

    selected = []
    rng = random.Random(args.seed)
    for site, (country, n) in SITE_PLAN.items():
        pool = by_site.get(site, [])
        if len(pool) < n:
            raise ValueError(f"{site}: need {n}, only {len(pool)} eligible individuals available")
        draw = rng.sample(pool, n)
        for r in draw:
            selected.append({
                "metaG_id": r["metaG_id"],
                "study_id": r["study_id"],
                "site": site,
                "country": country,
            })
        print(f"{site} ({country}): selected {n} of {len(pool)} eligible")

    with open(args.output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["metaG_id", "study_id", "site", "country"], delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    print(f"\nWrote {len(selected)} selected individuals to {args.output_path}")
    print(f"Seed: {args.seed}")


if __name__ == "__main__":
    main()
