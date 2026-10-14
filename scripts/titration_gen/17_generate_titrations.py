#!/usr/bin/env python3
"""
Subsample the bacterial and human read pools (16_build_tagged_pools.py) into
9 titration-point FASTA files at human content 10%-90% (10% steps), each a
fixed total_reads draw mixing the SAME bacterial community pool with the
SAME human pool - only the mixing ratio changes between files, isolating
that as the one variable under test.

Output: results/titrations/gut_complex_titre<NN>pct_human.fasta
        results/manifest.tsv

--content-pool/--human-pool/--output-dir/--filename-prefix/--manifest-path
let this mix any content pool (bacterial, or a viral pool from
20_extract_viral_pool.py) against any human pool (the default Forbes-27
pool, or the continental-African ablation's restricted pool) without
duplicating this script - see README "Human background panels"
section.
"""
import argparse
import csv
import random
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_BACTERIAL_POOL_PATH = REPO_ROOT / "results" / "pools" / "bacterial_pool.fasta"
DEFAULT_HUMAN_POOL_PATH = REPO_ROOT / "results" / "pools" / "human_pool.fasta"
DEFAULT_TITRATIONS_DIR = REPO_ROOT / "results" / "titrations"
DEFAULT_MANIFEST_PATH = REPO_ROOT / "results" / "manifest.tsv"
DEFAULT_FILENAME_PREFIX = "gut_complex"

TOTAL_READS = 200_000
HUMAN_PERCENTAGES = range(10, 100, 10)  # 10, 20, ..., 90
BASE_SEED = 632741178  # same repo-wide seed as the CAMISIM run; +pct per file below


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--content-pool", type=Path, default=DEFAULT_BACTERIAL_POOL_PATH,
                    help="Non-human content pool FASTA (bacterial or viral) (default: bacterial_pool.fasta)")
    p.add_argument("--human-pool", type=Path, default=DEFAULT_HUMAN_POOL_PATH)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_TITRATIONS_DIR)
    p.add_argument("--filename-prefix", type=str, default=DEFAULT_FILENAME_PREFIX)
    p.add_argument("--manifest-path", type=Path, default=DEFAULT_MANIFEST_PATH)
    p.add_argument("--content-label", type=str, default="bacterial",
                    help="Used to name the manifest's read-count column, e.g. 'viral' -> n_viral_reads (default: bacterial)")
    return p.parse_args()


def load_pool(path):
    records = []
    with open(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip("\n")
            tag = header[1:].rstrip("\n").split(" ", 1)[1]
            records.append((tag, seq))
    return records


def write_titration(path, records, rng):
    rng.shuffle(records)
    with open(path, "w") as f:
        for i, (tag, seq) in enumerate(records):
            f.write(f">read_{i} {tag}\n{seq}\n")


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    content_count_col = f"n_{args.content_label}_reads"

    print(f"Loading content pool from {args.content_pool} ...")
    content_records = load_pool(args.content_pool)
    print(f"Loading human pool from {args.human_pool} ...")
    human_records = load_pool(args.human_pool)
    print(f"{len(content_records)} {args.content_label} reads, {len(human_records)} human reads available.")

    manifest_rows = []
    for pct in HUMAN_PERCENTAGES:
        n_human = round(TOTAL_READS * pct / 100)
        n_content = TOTAL_READS - n_human
        if n_human > len(human_records) or n_content > len(content_records):
            raise ValueError(
                f"Pool too small for titre {pct}%: need {n_human} human / {n_content} {args.content_label}, "
                f"have {len(human_records)} / {len(content_records)}"
            )

        seed = BASE_SEED + pct
        rng = random.Random(seed)
        human_draw = rng.sample(human_records, n_human)
        content_draw = rng.sample(content_records, n_content)
        combined = human_draw + content_draw

        filename = f"{args.filename_prefix}_titre{pct}pct_human.fasta"
        write_titration(args.output_dir / filename, combined, rng)
        print(f"{filename}: {n_human} human + {n_content} {args.content_label} = {len(combined)} reads (seed={seed})")

        manifest_rows.append(
            {
                "filename": filename,
                "human_pct_target": pct,
                "total_reads": len(combined),
                "n_human_reads": n_human,
                content_count_col: n_content,
                "seed": seed,
            }
        )

    with open(args.manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote {args.manifest_path}")


if __name__ == "__main__":
    main()
