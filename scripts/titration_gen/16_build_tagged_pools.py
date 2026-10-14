#!/usr/bin/env python3
"""
Build tagged read pools that 17_generate_titrations.py samples from:

  results/pools/bacterial_pool.fasta - every CAMISIM-simulated bacterial
      read (both mates, independently), header field 2 = the genome's
      RefSeq assembly accession (e.g. GCF_000025985.1).
  results/pools/human_pool.fasta - every downloaded 1000 Genomes read
      (both mates, independently), header field 2 = its ENA run accession
      (e.g. ERR3239334), matching the existing project's human-read
      ground-truth tagging convention.

CAMISIM's `simulate_reads_art` process (read_simulators/read_simulator_art.nf)
publishes one FASTQ pair per genome, per sample, independently of the final
merged per-sample FASTQ: sample<id>_<genome_id>1.fq.gz / ...2.fq.gz. Reading
those directly (rather than the merged file) means read-to-genome
attribution comes for free from the filename - no need to parse ART's own
internal read headers.

--skip-bacterial / configurable human paths let this build just a human
pool against a different accessions table (e.g. the continental-African
ablation's restricted population set) without touching the unchanged
bacterial pool or re-running CAMISIM - see README "Human background
panels" section.
"""
import argparse
import csv
import gzip
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_GENOME_MANIFEST_PATH = REPO_ROOT / "results" / "genome_manifest.tsv"
DEFAULT_HUMAN_ACCESSIONS_PATH = REPO_ROOT / "results" / "human_accessions.tsv"
DEFAULT_CAMISIM_SAMPLE_DIR = REPO_ROOT / "camisim_out" / "sample_0" / "reads" / "fastq"
DEFAULT_HUMAN_READS_DIR = REPO_ROOT / "reads" / "human"
DEFAULT_BACTERIAL_POOL_PATH = REPO_ROOT / "results" / "pools" / "bacterial_pool.fasta"
DEFAULT_HUMAN_POOL_PATH = REPO_ROOT / "results" / "pools" / "human_pool.fasta"
DEFAULT_POOL_SUMMARY_PATH = REPO_ROOT / "results" / "pool_summary.tsv"

SAMPLE_ID = 0  # number_of_samples=1 -> CAMISIM always emits "sample_0"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genome-manifest", type=Path, default=DEFAULT_GENOME_MANIFEST_PATH)
    p.add_argument("--camisim-sample-dir", type=Path, default=DEFAULT_CAMISIM_SAMPLE_DIR)
    p.add_argument("--bacterial-pool-output", type=Path, default=DEFAULT_BACTERIAL_POOL_PATH)
    p.add_argument("--human-accessions-path", type=Path, default=DEFAULT_HUMAN_ACCESSIONS_PATH)
    p.add_argument("--human-reads-dir", type=Path, default=DEFAULT_HUMAN_READS_DIR)
    p.add_argument("--human-pool-output", type=Path, default=DEFAULT_HUMAN_POOL_PATH)
    p.add_argument("--pool-summary-path", type=Path, default=DEFAULT_POOL_SUMMARY_PATH)
    p.add_argument("--skip-bacterial", action="store_true", help="Don't rebuild the bacterial pool (e.g. it's unchanged for a new human-only ablation arm)")
    p.add_argument("--skip-human", action="store_true", help="Don't rebuild the human pool")
    return p.parse_args()


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def iter_fastq_sequences(path):
    with open_maybe_gzip(path) as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline().rstrip("\n")
            f.readline()  # '+' line
            f.readline()  # quality line
            yield seq


def write_pool_entry(out_f, counter, tag, seq):
    out_f.write(f">read_{counter} {tag}\n{seq}\n")
    return counter + 1


def build_bacterial_pool(genome_manifest_path, camisim_sample_dir, bacterial_pool_path):
    with open(genome_manifest_path, newline="") as f:
        genome_ids = [row["genome_ID"] for row in csv.DictReader(f, delimiter="\t")]

    counter = 0
    per_genome_counts = {}
    with open(bacterial_pool_path, "w") as out_f:
        for genome_id in genome_ids:
            mate_paths = [
                camisim_sample_dir / f"sample{SAMPLE_ID}_{genome_id}1.fq.gz",
                camisim_sample_dir / f"sample{SAMPLE_ID}_{genome_id}2.fq.gz",
            ]
            n_reads = 0
            for mate_path in mate_paths:
                if not mate_path.exists():
                    raise FileNotFoundError(
                        f"Expected CAMISIM output not found: {mate_path}\n"
                        f"Contents of {camisim_sample_dir}: "
                        f"{sorted(p.name for p in camisim_sample_dir.glob('*'))[:20]}"
                    )
                for seq in iter_fastq_sequences(mate_path):
                    counter = write_pool_entry(out_f, counter, genome_id, seq)
                    n_reads += 1
            per_genome_counts[genome_id] = n_reads
            print(f"{genome_id}: {n_reads} reads")

    print(f"Bacterial pool: {counter} reads from {len(genome_ids)} genomes -> {bacterial_pool_path}")
    return counter, per_genome_counts


def build_human_pool(human_accessions_path, human_reads_dir, human_pool_path):
    with open(human_accessions_path, newline="") as f:
        run_accessions = [row["run_accession"] for row in csv.DictReader(f, delimiter="\t")]

    counter = 0
    per_run_counts = {}
    with open(human_pool_path, "w") as out_f:
        for run_acc in run_accessions:
            mate_paths = [
                human_reads_dir / f"{run_acc}_1.fastq",
                human_reads_dir / f"{run_acc}_2.fastq",
            ]
            n_reads = 0
            for mate_path in mate_paths:
                if not mate_path.exists():
                    raise FileNotFoundError(f"Expected human read file not found: {mate_path}")
                for seq in iter_fastq_sequences(mate_path):
                    counter = write_pool_entry(out_f, counter, run_acc, seq)
                    n_reads += 1
            per_run_counts[run_acc] = n_reads
            print(f"{run_acc}: {n_reads} reads")

    print(f"Human pool: {counter} reads from {len(run_accessions)} runs -> {human_pool_path}")
    return counter, per_run_counts


def main():
    args = parse_args()
    args.bacterial_pool_output.parent.mkdir(parents=True, exist_ok=True)
    args.human_pool_output.parent.mkdir(parents=True, exist_ok=True)

    per_genome_counts, per_run_counts = {}, {}
    n_bacterial, n_human = 0, 0

    if not args.skip_bacterial:
        n_bacterial, per_genome_counts = build_bacterial_pool(
            args.genome_manifest, args.camisim_sample_dir, args.bacterial_pool_output
        )
    if not args.skip_human:
        n_human, per_run_counts = build_human_pool(
            args.human_accessions_path, args.human_reads_dir, args.human_pool_output
        )

    with open(args.pool_summary_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["pool", "source_id", "n_reads"])
        for genome_id, n in per_genome_counts.items():
            writer.writerow(["bacterial", genome_id, n])
        for run_acc, n in per_run_counts.items():
            writer.writerow(["human", run_acc, n])

    print(f"Total: {n_bacterial} bacterial + {n_human} human reads pooled.")
    print(f"Wrote per-source breakdown to {args.pool_summary_path}")


if __name__ == "__main__":
    main()
