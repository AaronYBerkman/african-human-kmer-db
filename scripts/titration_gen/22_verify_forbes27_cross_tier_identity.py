#!/usr/bin/env python3
"""
Verify the design guarantee the Forbes-27 tier depends on: the microbial x
Forbes-27 arm's human pool (results/pools/human_pool.fasta, built from
data/forbes_human_ancestry_accessions.tsv) and the viral x Forbes-27 arm
(the ERP169281 accessions, used as-is - see metadata/viral_reads.tsv) must
draw from the SAME 27 human individuals, or any microbial-vs-viral
difference observed later could be an artifact of different human
backgrounds rather than the content type.

This is checkable directly: Forbes' own mixed ERP169281 samples carry their
embedded human reads' ENA run accession in field 2 (same convention this
repo uses throughout - see README "Tagging convention"). Scans the raw
viral accessions already downloaded by 20_extract_viral_pool.py
(reads/viral_raw/, gitignored - re-run 20 first if not present) and
extracts every embedded human ERR accession, then diffs that set against
our own Forbes-27 accession list.

Exits non-zero if the sets differ for any of the 27 source accessions.
Writes results/forbes27_cross_tier_identity_check.tsv (one row per source
accession: embedded human accession count, whether it matches our list).
"""
import csv
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUR_ACCESSIONS_PATH = REPO_ROOT / "data" / "forbes_human_ancestry_accessions.tsv"
VIRAL_RAW_DIR = REPO_ROOT / "reads" / "viral_raw"
VIRAL_POOL_MANIFEST_PATH = REPO_ROOT / "results" / "viral_pool_manifest.tsv"
REPORT_PATH = REPO_ROOT / "results" / "forbes27_cross_tier_identity_check.tsv"

ERR_PREFIX_PATTERN = re.compile(r"^ERR(\d+)")


def load_our_accessions():
    with open(OUR_ACCESSIONS_PATH, newline="") as f:
        return {row["ena_run_accession"] for row in csv.DictReader(f, delimiter="\t")}


def embedded_human_accessions(fastq_path):
    found = set()
    with open(fastq_path) as f:
        for i, line in enumerate(f):
            if i % 4 != 0:
                continue
            fields = line[1:].rstrip("\n").split(" ")
            if len(fields) < 2:
                continue
            m = ERR_PREFIX_PATTERN.match(fields[1])
            if m:
                found.add(f"ERR{m.group(1)}")
    return found


def main():
    our_accessions = load_our_accessions()
    print(f"Our Forbes-27 human pool accessions ({len(our_accessions)}): {sorted(our_accessions)}")

    with open(VIRAL_POOL_MANIFEST_PATH, newline="") as f:
        source_accessions = [row["ena_run_accession"] for row in csv.DictReader(f, delimiter="\t")]

    if not VIRAL_RAW_DIR.exists():
        print(f"ERROR: {VIRAL_RAW_DIR} not found - run scripts/titration_gen/20_extract_viral_pool.py first.", file=sys.stderr)
        sys.exit(2)

    report_rows = []
    all_match = True
    for run_acc in source_accessions:
        embedded = set()
        for mate in (1, 2):
            path = VIRAL_RAW_DIR / f"{run_acc}_{mate}.fastq"
            if not path.exists():
                print(f"ERROR: {path} missing.", file=sys.stderr)
                sys.exit(2)
            embedded |= embedded_human_accessions(path)

        missing_from_ours = embedded - our_accessions
        extra_in_ours = our_accessions - embedded
        matches = not missing_from_ours and not extra_in_ours
        all_match &= matches

        report_rows.append(
            {
                "viral_source_accession": run_acc,
                "n_embedded_human_accessions": len(embedded),
                "n_our_accessions": len(our_accessions),
                "identical_set": matches,
                "missing_from_our_pool": ";".join(sorted(missing_from_ours)) or "",
                "extra_in_our_pool_not_embedded": ";".join(sorted(extra_in_ours)) or "",
            }
        )
        print(f"  {run_acc}: {len(embedded)} embedded human accessions, identical_set={matches}")

    with open(REPORT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(report_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows)
    print(f"Wrote {REPORT_PATH}")

    if not all_match:
        print("\nFAILED: microbial x Forbes-27 and viral x Forbes-27 do NOT share the same 27 individuals.", file=sys.stderr)
        sys.exit(1)
    print(f"\nOK: all {len(source_accessions)} viral-benchmark source accessions embed exactly our {len(our_accessions)} "
          f"Forbes-27 human accessions - microbial x Forbes-27 and viral x Forbes-27 share an identical human pool.")


if __name__ == "__main__":
    main()
