"""Shared ENA fetch helpers and constants used by 15_download_human_reads.py,
19_compile_provenance_metadata.py, and 20_extract_viral_pool.py - not a
numbered pipeline step itself."""
import csv
import subprocess
import urllib.request

# Exact viral accessions actually used: 3 sample types x human-titre deciles
# 1-9 x iteration 1. Single source of truth for "which viral accessions did
# we use" - both 19_compile_provenance_metadata.py and 20_extract_viral_pool.py
# filter metadata/viral_reads.tsv against these.
BENCHMARK_VIRUS_TYPES = {"clinical_rsv_b", "clinical_sars_cov_2_high_ct", "h3n2_control"}
BENCHMARK_ITERATION = 1


def filereport_url(run_accession):
    return (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={run_accession}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"
    )


def fetch_fastq_urls(run_accession):
    with urllib.request.urlopen(filereport_url(run_accession), timeout=60) as resp:
        text = resp.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    fastq_ftp = rows[0]["fastq_ftp"]
    urls = fastq_ftp.split(";")
    if len(urls) != 2:
        raise ValueError(f"{run_accession}: expected 2 paired-end fastq files, got {urls}")
    return [normalize_url(u) for u in urls]


def normalize_url(u):
    return u if u.startswith(("http://", "https://", "ftp://")) else f"https://{u}"


def stream_partial_fastq(url, n_reads, dest_path, max_attempts=4):
    n_lines = n_reads * 4
    # Note: a successful run normally exits 141 (SIGPIPE) here, not 0 - once
    # `head` has read enough lines it closes the pipe, and curl/zcat upstream
    # get SIGPIPE. That's the expected/correct outcome of this pattern, not a
    # failure, so success is judged purely on actual line count.
    cmd = f"curl -s --max-time 600 '{url}' | zcat | head -n {n_lines}"
    for attempt in range(1, max_attempts + 1):
        with open(dest_path, "w") as out:
            subprocess.run(["bash", "-c", cmd], stdout=out)
        actual_lines = sum(1 for _ in open(dest_path))
        if actual_lines == n_lines:
            return
        print(f"  attempt {attempt}/{max_attempts} got {actual_lines}/{n_lines} lines, retrying ...")
    raise RuntimeError(f"Failed to fetch {n_reads} reads from {url} after {max_attempts} attempts")
