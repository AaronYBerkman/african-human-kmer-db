#!/usr/bin/env python3
"""
Pull a reproducible pool of real human WGS reads replicating the
mixed-ancestry composition Forbes et al. (2025) used as the human
background in their viral titration benchmark: one representative ENA run
per 1000 Genomes population (data/forbes_human_ancestry_accessions.tsv, 27
rows covering the 26 standard populations - CEU appears twice in the
source table and is kept as given).

Only a partial read pool is needed (not full ~30x WGS runs, which are tens
of GB each per individual): for each run, the remote fastq.gz is streamed
through `zcat | head` so the download stops as soon as enough reads have
been read, without ever fetching the full file.

CLI args support restricting to a population subset at greater per-run
depth (e.g. the continental-African ablation: 5 of the 27 populations,
sampled deeper so the smaller pool still supports the same total-reads
titration design) without duplicating this script - see README "Human
background panels" section.
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _ena_utils import fetch_fastq_urls, stream_partial_fastq  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_ACCESSIONS_TABLE_PATH = REPO_ROOT / "data" / "forbes_human_ancestry_accessions.tsv"
DEFAULT_READS_PER_MATE_PER_RUN = 15_000  # 27 runs x 2 mates x 15,000 = 810,000 pool entries
DEFAULT_READS_DIR = REPO_ROOT / "reads" / "human"
DEFAULT_ACCESSIONS_PATH = REPO_ROOT / "results" / "human_accessions.tsv"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--accessions-table", type=Path, default=DEFAULT_ACCESSIONS_TABLE_PATH,
                    help="Source population/run accession TSV (default: data/forbes_human_ancestry_accessions.tsv)")
    p.add_argument("--population-filter", type=str, default=None,
                    help="Comma-separated population codes to keep (default: all rows in the accessions table)")
    p.add_argument("--reads-per-mate", type=int, default=DEFAULT_READS_PER_MATE_PER_RUN,
                    help=f"Reads to stream per mate per run (default: {DEFAULT_READS_PER_MATE_PER_RUN})")
    p.add_argument("--output-dir", type=Path, default=DEFAULT_READS_DIR,
                    help="Directory to write downloaded FASTQ files (default: reads/human)")
    p.add_argument("--manifest-path", type=Path, default=DEFAULT_ACCESSIONS_PATH,
                    help="Output manifest TSV path (default: results/human_accessions.tsv)")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with open(args.accessions_table, newline="") as f:
        populations = list(csv.DictReader(f, delimiter="\t"))

    if args.population_filter:
        keep = {code.strip() for code in args.population_filter.split(",")}
        populations = [p for p in populations if p["population_code"] in keep]
        missing = keep - {p["population_code"] for p in populations}
        if missing:
            raise ValueError(f"Population code(s) not found in {args.accessions_table.name}: {sorted(missing)}")

    print(f"{len(populations)} population/run rows loaded from {args.accessions_table.name}"
          f"{f' (filtered to {args.population_filter})' if args.population_filter else ''}.")

    manifest_rows = []
    for pop in populations:
        run_acc = pop["ena_run_accession"]
        url_1, url_2 = fetch_fastq_urls(run_acc)
        dest_1 = args.output_dir / f"{run_acc}_1.fastq"
        dest_2 = args.output_dir / f"{run_acc}_2.fastq"

        for url, dest in [(url_1, dest_1), (url_2, dest_2)]:
            if dest.exists():
                print(f"Already have {dest.name}, skipping.")
                continue
            print(f"[{pop['population_code']}] Streaming {args.reads_per_mate} reads from {url} ...")
            stream_partial_fastq(url, args.reads_per_mate, dest)

        manifest_rows.append(
            {
                "population_code": pop["population_code"],
                "population_name": pop["population_name"],
                "run_accession": run_acc,
                "ena_study_accession": pop["ena_study_accession"],
                "fastq_1_url": url_1,
                "fastq_2_url": url_2,
                "reads_taken_per_mate": args.reads_per_mate,
            }
        )

    with open(args.manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    n_pops = len({p["population_code"] for p in populations})
    print(f"Wrote {args.manifest_path} ({len(manifest_rows)} runs across {n_pops} populations)")


if __name__ == "__main__":
    main()
