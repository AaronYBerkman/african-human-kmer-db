#!/usr/bin/env python3
"""
Extract pure viral reads from 27 ERP169281 accessions (Forbes et al. 2025):
3 sample types x titre deciles 1-9 x iteration 1 - clinical_rsv_b,
clinical_sars_cov_2_high_ct, h3n2_control. Confirmed against
metadata/viral_reads.tsv, not re-transcribed by hand.

Each accession is downloaded in full (these are small - tens to hundreds of
MB, not the 30x WGS scale of the human reads; checked via ENA's filereport
API first) and filtered by this project's ground-truth convention: a read is
human if field 2 of its header STARTS WITH an ERR accession (a prefix match
against Forbes' own read IDs like "ERR3239334.100342045/1", not the
exact-equality match this repo's own simplified tags use elsewhere). Human
reads are discarded; viral reads are kept and retagged with their
virus_sample_type (discarding Forbes' own chromosome-coordinate-style tag,
not needed here). The resulting pure-viral pool is then remixed with this
project's own human read pool (see 24_african25_benchmark.sh and
25_viral_mixed_benchmark.sh), rather than using Forbes' pre-mixed files
directly, which embed a human background of different, uncontrolled
composition.

Output: results/pools/viral_pool.fasta (all 3 types pooled, tagged by type)
        results/viral_pool_manifest.tsv (per-source-accession breakdown)
        results/viral_content_manifest.tsv (content_tag column, for
            18_verify_ground_truth.py --content-manifest)
"""
import csv
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _ena_utils import BENCHMARK_ITERATION, BENCHMARK_VIRUS_TYPES, fetch_fastq_urls  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
VIRAL_READS_MANIFEST_PATH = REPO_ROOT / "metadata" / "viral_reads.tsv"

RAW_DOWNLOAD_DIR = REPO_ROOT / "reads" / "viral_raw"
POOLS_DIR = REPO_ROOT / "results" / "pools"
VIRAL_POOL_PATH = POOLS_DIR / "viral_pool.fasta"
POOL_MANIFEST_PATH = REPO_ROOT / "results" / "viral_pool_manifest.tsv"
CONTENT_MANIFEST_PATH = REPO_ROOT / "results" / "viral_content_manifest.tsv"

# Forbes' own convention (not this repo's simplified exact-match tags): a
# read is human if field 2 STARTS WITH an ERR accession, e.g.
# "ERR3239334.100342045/1" - confirmed by inspecting actual downloaded reads.
ERR_PREFIX_PATTERN = re.compile(r"^ERR\d+")


def load_benchmark_accessions():
    with open(VIRAL_READS_MANIFEST_PATH, newline="") as f:
        rows = [
            row for row in csv.DictReader(f, delimiter="\t")
            if row["virus_sample_type"] in BENCHMARK_VIRUS_TYPES
            and int(row["iteration"]) == BENCHMARK_ITERATION
        ]
    if len(rows) != len(BENCHMARK_VIRUS_TYPES) * 9:
        raise ValueError(f"Expected {len(BENCHMARK_VIRUS_TYPES) * 9} accessions, found {len(rows)}")
    return rows


def download_full_fastq(url, dest_path):
    cmd = f"curl -s --max-time 1800 '{url}' | zcat"
    with open(dest_path, "w") as out:
        subprocess.run(["bash", "-c", cmd], stdout=out, check=True)


def iter_fastq_records(path):
    with open(path) as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline().rstrip("\n")
            f.readline()  # '+' line
            f.readline()  # quality line
            yield header, seq


def main():
    RAW_DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    POOLS_DIR.mkdir(parents=True, exist_ok=True)

    rows = load_benchmark_accessions()
    print(f"{len(rows)} source accessions loaded from {VIRAL_READS_MANIFEST_PATH.name} "
          f"({BENCHMARK_VIRUS_TYPES}, iteration={BENCHMARK_ITERATION}).")

    pool_manifest_rows = []
    counter = 0
    with open(VIRAL_POOL_PATH, "w") as out_f:
        for row in rows:
            virus_type = row["virus_sample_type"]
            titre = row["human_titre_decile"]
            run_acc = row["ena_run_accession"]
            url_1, url_2 = fetch_fastq_urls(run_acc)

            n_viral, n_human = 0, 0
            for mate_idx, url in enumerate([url_1, url_2], start=1):
                raw_path = RAW_DOWNLOAD_DIR / f"{run_acc}_{mate_idx}.fastq"
                if not raw_path.exists():
                    print(f"[{virus_type} titre={titre}] Downloading {run_acc} mate {mate_idx} ...")
                    download_full_fastq(url, raw_path)
                for header, seq in iter_fastq_records(raw_path):
                    fields = header[1:].rstrip("\n").split(" ")
                    tag = fields[1] if len(fields) > 1 else ""
                    if ERR_PREFIX_PATTERN.match(tag):
                        n_human += 1
                        continue
                    out_f.write(f">read_{counter} {virus_type}\n{seq}\n")
                    counter += 1
                    n_viral += 1

            total = n_viral + n_human
            viral_fraction = n_viral / total if total else 0
            print(f"  {run_acc}: {n_viral} viral kept, {n_human} human discarded ({viral_fraction:.1%} viral)")

            pool_manifest_rows.append(
                {
                    "virus_sample_type": virus_type,
                    "human_titre_decile": titre,
                    "iteration": BENCHMARK_ITERATION,
                    "ena_run_accession": run_acc,
                    "n_viral_reads_kept": n_viral,
                    "n_human_reads_discarded": n_human,
                    "viral_fraction": round(viral_fraction, 4),
                }
            )

    with open(POOL_MANIFEST_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(pool_manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(pool_manifest_rows)
    print(f"Wrote {POOL_MANIFEST_PATH}")

    with open(CONTENT_MANIFEST_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["content_tag"], delimiter="\t")
        writer.writeheader()
        writer.writerows({"content_tag": vt} for vt in sorted(BENCHMARK_VIRUS_TYPES))
    print(f"Wrote {CONTENT_MANIFEST_PATH}")

    total_viral = sum(r["n_viral_reads_kept"] for r in pool_manifest_rows)
    print(f"Viral pool: {total_viral} reads from {len(rows)} accessions -> {VIRAL_POOL_PATH}")


if __name__ == "__main__":
    main()
