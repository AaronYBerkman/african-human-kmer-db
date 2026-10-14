#!/usr/bin/env python3
"""
Spot-check ground-truth read tagging across a set of titration FASTAs:

  - every read whose header field 2 matches ^ERR\\d+ must be one of the
    known human ENA run accessions, and only those reads
  - every other read's field 2 must be one of the known content IDs
    (bacterial RefSeq genome accessions, or virus sample-type tags), and
    only those reads
  - the achieved human fraction per file must match its target (manifest.tsv)

Exits non-zero on any mismatch.

CLI args let this verify any arm (the default Forbes-27 microbial arm, or
the continental-African ablation's new arms) without duplicating this
script - see README "Human background panels" section.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_MANIFEST_PATH = REPO_ROOT / "results" / "manifest.tsv"
DEFAULT_CONTENT_MANIFEST_PATH = REPO_ROOT / "results" / "genome_manifest.tsv"
DEFAULT_HUMAN_ACCESSIONS_PATH = REPO_ROOT / "results" / "human_accessions.tsv"
DEFAULT_TITRATIONS_DIR = REPO_ROOT / "results" / "titrations"

ERR_PATTERN = re.compile(r"^ERR\d+$")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--titrations-dir", type=Path, default=DEFAULT_TITRATIONS_DIR)
    p.add_argument("--manifest-path", type=Path, default=DEFAULT_MANIFEST_PATH)
    p.add_argument("--content-manifest", type=Path, default=DEFAULT_CONTENT_MANIFEST_PATH,
                    help="TSV with the known content IDs (default: genome_manifest.tsv)")
    p.add_argument("--content-id-column", type=str, default="genome_ID",
                    help="Column in --content-manifest holding the content tag (default: genome_ID)")
    p.add_argument("--content-label", type=str, default="bacterial",
                    help="Matches 17_generate_titrations.py's --content-label, to find the right manifest.tsv column (default: bacterial)")
    p.add_argument("--human-accessions-path", type=Path, default=DEFAULT_HUMAN_ACCESSIONS_PATH)
    return p.parse_args()


def load_known_ids(content_manifest_path, content_id_column, human_accessions_path):
    with open(content_manifest_path, newline="") as f:
        content_ids = {row[content_id_column] for row in csv.DictReader(f, delimiter="\t")}
    with open(human_accessions_path, newline="") as f:
        human_ids = {row["run_accession"] for row in csv.DictReader(f, delimiter="\t")}
    return content_ids, human_ids


def check_file(path, content_ids, human_ids):
    n_human = 0
    n_content = 0
    bad_lines = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            fields = line[1:].rstrip("\n").split(" ")
            if len(fields) < 2:
                bad_lines.append(line)
                continue
            tag = fields[1]
            is_err_shaped = bool(ERR_PATTERN.match(tag))
            if is_err_shaped:
                if tag not in human_ids:
                    bad_lines.append(line)
                else:
                    n_human += 1
            else:
                if tag not in content_ids:
                    bad_lines.append(line)
                else:
                    n_content += 1
    return n_human, n_content, bad_lines


def main():
    args = parse_args()
    content_count_col = f"n_{args.content_label}_reads"
    content_ids, human_ids = load_known_ids(args.content_manifest, args.content_id_column, args.human_accessions_path)

    with open(args.manifest_path, newline="") as f:
        manifest_rows = list(csv.DictReader(f, delimiter="\t"))

    all_ok = True
    print(f"{'file':<45} {'target%':>8} {'achieved%':>10} {'n_human':>9} {'n_content':>10} {'bad':>5}")
    for row in manifest_rows:
        path = args.titrations_dir / row["filename"]
        n_human, n_content, bad_lines = check_file(path, content_ids, human_ids)
        total = n_human + n_content
        achieved_pct = 100 * n_human / total if total else 0
        target_pct = float(row["human_pct_target"])
        ok = (
            not bad_lines
            and n_human == int(row["n_human_reads"])
            and n_content == int(row[content_count_col])
            and abs(achieved_pct - target_pct) < 0.5
        )
        all_ok &= ok
        flag = "" if ok else "  <-- MISMATCH"
        print(
            f"{row['filename']:<45} {target_pct:>8.1f} {achieved_pct:>10.2f} "
            f"{n_human:>9} {n_content:>10} {len(bad_lines):>5}{flag}"
        )
        if bad_lines:
            for line in bad_lines[:5]:
                print(f"    bad header: {line.rstrip()}", file=sys.stderr)

    if not all_ok:
        print("\nFAILED: ground-truth tagging verification found mismatches.", file=sys.stderr)
        sys.exit(1)
    print(f"\nOK: all {len(manifest_rows)} titration files have exact, unambiguous ground-truth tagging.")


if __name__ == "__main__":
    main()
