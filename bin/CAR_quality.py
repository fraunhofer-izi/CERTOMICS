#!/usr/bin/env python3
"""
CAR-specific quality control of scRNA-seq data
Author: Christina Kuhn
Date: 2025-13-03
"""
import argparse
import csv
import logging
import os
from typing import List, Tuple
import pyfaidx  # to handle fasta files https://pythonhosted.org/pyfaidx/
import pysam  # to handle bam files https://pysam.readthedocs.io/en/latest/api.html

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def parse_arguments() -> argparse.Namespace:
    """
    Parses command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Perform quality control on scRNA-seq CellRanger \
                                     output specific for assesment of CAR therapy efficacy.")
    parser.add_argument("--sample_names",
                        nargs="+",
                        type=str,
                        required=True,
                        help="List of samples names to be analyzed. Same order as Bam files.")
    parser.add_argument("--bam_files",
                        nargs="+",
                        type=str,
                        required=True,
                        help="Path(s) to the BAM file(s) of the mapping.")
    parser.add_argument("--CAR_fasta_file",
                        type=str, required=True,
                        help="Path to the FASTA file of the CAR construct.")
    parser.add_argument("--CAR_gtf_file",
                        type=str, required=True,
                        help="Path to the GTF file of the CAR domain annotations.")

    args = parser.parse_args()

    #log each argument passed
    logging.info("Sample names: %s", args.sample_names)
    logging.info("BAM files: %d", len(args.bam_files))
    logging.info("CAR FASTA path: %s", args.CAR_fasta_file)
    logging.info("CAR GTF path: %s", args.CAR_gtf_file)

    return args

def load_data(
        bam_paths: List[str],
        fasta_path: str,
        gtf_path_car: str) -> Tuple[List[pysam.AlignmentFile], List[str], pyfaidx.Fasta, str]:
    """
    Load data from BAM, FASTA, and GTF files.

    Args:
        bam_paths (List[str]): Paths to the BAM files.
        fasta_path (List[str]): Paths to the FASTA file.
        gtf_path_car (List[str]): Paths to the GTF file for the CAR construct.

    Returns:
        Tuple[List[pysam.AlignmentFile], List[pyfaidx.Fasta], List[str], List[str]]:
            - List of pysam.AlignmentFile objects for each BAM file.
            - List of pyfaidx.Fasta objects for each FASTA file.
            - List of paths for each CAR GTF file verified to exist.
            - List of paths for each concatenated annotation GTF file verified to exist.
    """
    bam_files = []
    for path in bam_paths:
        try:
            bam_file = pysam.AlignmentFile(path, "rb")
            bam_files.append(bam_file)
            logging.info("Loaded BAM file: %s", bam_file)
        except FileNotFoundError:
            logging.error("Error: The BAM file at %s was not found.", path)
            raise

    fasta_path = verify_file_exists(fasta_path, "FASTA file")
    fasta_file = pyfaidx.Fasta(fasta_path)
    verified_gtf_path_car = verify_file_exists(gtf_path_car, "CAR GTF file")

    return bam_files, fasta_file, verified_gtf_path_car

def verify_file_exists(file_path: str, file_description: str) -> str:
    """Verifies that a file exists and logs the result."""
    if os.path.exists(file_path):
        logging.info("%s verified at %s", file_description, file_path)
        return file_path

    logging.error("Error: %s at %s was not found.", file_description, file_path)
    raise FileNotFoundError(f"Error: {file_description} at {file_path} was not found.")

def get_car_start_end_from_gtfs(
        gtf_path_car: str
) -> Tuple[str, int, int]:
    """Get the CAR construct from the GTF file.

    Args:
        gtf_path_car (str): Path to the GTF file.

    Raises:
        ValueError: If the CAR region is not identified in the GTF file.

    Returns:
        Tuple[str, int, int]: CAR construct name, minimum start, and maximum end positions.
    """

    start = float("inf")
    end = float("-inf")
    car = ""

    with open(gtf_path_car, encoding="utf-8") as f:
        lines = f.readlines()

        if not lines:
            raise ValueError(f"Error: GTF file {gtf_path_car} is empty.")

        first_line = lines[0].strip().split("\t")
        last_line = lines[-1].strip().split("\t")

        car = str(first_line[0])
        start = int(first_line[3])
        end = int(last_line[4])

    logging.info("CAR construct: %s, Start: %d, End: %d", car, start, end)
    
    return car, start, end  # ✅ Return a tuple instead of a list


def generate_quality_metrics(
    bams: List[pysam.AlignmentFile],
    sample_names: List[str],
    gtf_path_car: str
):
    """
    Analyze multiple BAM files to perform quality control specific to CAR-T constructs.
    
    This function calculates the absolute number of reads mapping to the CAR region,
    as well as the coverage across this region for each sample.

    Args:
        bams (List[pysam.AlignmentFile]): List of BAM file objects for reads.
        sample_names (List[str]): List of sample names corresponding to each BAM file.
        gtf_path_car (str): Path to the GTF file of the CAR domain annotations.
    """
    logging.info("Calculating metrics and saving as CSVs...")
    logging.info("GTF Path: %s", gtf_path_car)
    logging.info("Samples: %s", sample_names)

    car_used, min_start, max_end  = get_car_start_end_from_gtfs(gtf_path_car)
    results_reads = {}
    results_coverage = {}
    results_coverage_unique = {}

    for bam, sample in zip(bams, sample_names):
        print(f"Processing sample {sample}")
        logging.info("Processing sample %s", sample)
        logging.info("Contig '%s' found in BAM file.", car_used)

        # Count reads in the CAR region
        total_reads, forward_reads, reverse_reads = count_reads_in_region(bam, car_used)
        results_reads[sample] = {
            "absolute_reads": total_reads,
            "forward_reads": forward_reads,
            "reverse_reads": reverse_reads
        }

        # Compute coverage metrics
        results_coverage[sample] = calculate_coverage(bam, car_used, min_start, max_end)
        results_coverage_unique[sample] = calculate_unique_coverage(
            bam, car_used, min_start, max_end
        )
    print(results_reads)
    # Save results as CSV files
    write_dict_to_csv("results_metrics_reads_CAR.csv", results_reads)

    with open("results_coverage_against_CAR.csv", "w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=results_coverage.keys())
        writer.writeheader()
        writer.writerow(results_coverage)

    with open("results_coverage_against_CAR_unique.csv", "w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=results_coverage_unique.keys())
        writer.writeheader()
        writer.writerow(results_coverage_unique)

def count_reads_in_region(
        bam: pysam.AlignmentFile,
        car_used: str
) -> Tuple[int, int, int]:
    """Counts total, forward, and reverse reads mapping to a specific CAR region.
    If region is not found, returns default values (0, 0, 0).
    """
    try:
        reads = list(bam.fetch(car_used))
        total_reads = len(reads)
        forward_reads = sum(1 for read in reads if not read.is_reverse)
        reverse_reads = total_reads - forward_reads
    except (ValueError, KeyError):  # e.g., region not found in BAM file
        total_reads = forward_reads = reverse_reads = 0

    return total_reads, forward_reads, reverse_reads

def calculate_coverage(
        bam: pysam.AlignmentFile,
        car_used: str,
        min_start: int,
        max_end: int
) -> List[int]:
    """Calculates coverage across the CAR region.
    If region is not found, returns coverage as all zeros.
    """
    coverage_counts = [0] * (max_end - min_start + 1)

    try:
        for read in bam.fetch(car_used, min_start, max_end):
            overlap_start = max(read.reference_start, min_start)
            overlap_end = min(read.reference_end, max_end)

            for j in range(overlap_start, overlap_end):
                coverage_counts[j - min_start] += 1
    except (ValueError, KeyError):
        print(f"Warning: CAR region '{car_used}' not found in BAM. Coverage set to zeros.")

    return coverage_counts

def calculate_unique_coverage(
        bam: pysam.AlignmentFile,
        car_used: str,
        min_start: int,
        max_end: int
) -> List[int]:
    """Calculates coverage of uniquely mapped reads in the CAR region.
    If region is not found, returns all-zero coverage.
    """
    coverage_counts = [0] * (max_end - min_start + 1)

    try:
        for read in bam.fetch(car_used, min_start, max_end):
            if read.has_tag("NH") and read.get_tag("NH") == 1:  # Only uniquely mapped reads
                overlap_start = max(read.reference_start, min_start)
                overlap_end = min(read.reference_end, max_end)

                for j in range(overlap_start, overlap_end):
                    coverage_counts[j - min_start] += 1
    except (ValueError, KeyError):
        print(f"Warning: CAR region '{car_used}' not found in BAM. Unique coverage set to zeros.")

    return coverage_counts

def write_dict_to_csv(
        filename: str,
        data: dict[str, dict[str, int]]):
    """Writes dictionary data to a CSV file."""
    with open(filename, "w", encoding="utf-8") as f:
        if not data:
            # If data is empty, write default header and one line with zeros
            fieldnames = ["sample","absolute_reads","forward_reads","reverse_reads"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({"sample": "NA", "absolute_reads": 0, "forward_reads": 0,"reverse_reads":0})  
        else:
            first_sample_key = next(iter(data))  # Get the first key
            fieldnames = ["sample"] + list(data[first_sample_key].keys())

            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for sample, metrics in data.items():
                writer.writerow({"sample": sample, **metrics})

def main():
    """
    Main function
    Example usage:
    python CAR_quality.py
        --sample_names "sample1" "sample2"
        --bam_files "/sample1/count/sample_alignments.bam" "/sample2/count/sample_alignments.bam"
        --CAR_fasta_file Idecel.fasta --CAR_gtf_file Anno_Idecel.gtf
    """
    args = parse_arguments()
    bams, CAR_fasta, CAR_anno = load_data(args.bam_files, args.CAR_fasta_file, args.CAR_gtf_file)
    print("Data loaded.")
    print("Samples:", args.sample_names)
    print("CAR fasta:", CAR_fasta)
    print("CAR anno:", CAR_anno)
    generate_quality_metrics(bams, args.sample_names, CAR_anno)

if __name__ == "__main__":
    main()
