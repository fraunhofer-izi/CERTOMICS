#!/usr/bin/env python3
"""
CAR-specific quality control plotting functions
Author: Christina Kuhn
Date: 2025-13-03
"""
import argparse
import csv
import logging
from typing import Tuple
import ast
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pygenomeviz import GenomeViz
from matplotlib.patches import Patch

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def parse_arguments() -> argparse.Namespace:
    """
    Parses command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Perform quality control on scRNA-seq CellRanger output \
        specific for assesment of CAR therapy efficacy.")
    parser.add_argument('--read_metrics',
                        type=str,
                        required=True,
                        help='The CSV file with read metrics against CAR')
    parser.add_argument('--coverage_all',
                        type=str,
                        required=True,
                        help='The CSV file with coverage metrics against CAR')
    parser.add_argument('--coverage_unique',
                        type=str,
                        required=True,
                        help='The CSV file with coverage metrics against CAR')
    parser.add_argument("--CAR_gtf_file",
                        type=str,
                        required=True,
                        help="Path to the GTF file of the CAR domain annotations.")

    args = parser.parse_args()

    # Log each argument passed
    logging.info("Path to read metrics: %s", args.read_metrics)
    logging.info("Path to coverage all: %s", args.coverage_all)
    logging.info("Path to coverage unique: %s", args.coverage_unique)
    logging.info("CAR GTF file: %s", args.CAR_gtf_file)

    return args

def read_data(read_metrics, coverage_all, coverage_unique):
    """_summary_

    Args:
        read_metrics (_type_): _description_
        coverage_all (_type_): _description_
        coverage_unique (_type_): _description_
    """

    # Open the CSV file
    with open(read_metrics, "r", encoding="utf-8") as f:
        # Create a CSV reader
        reader = csv.DictReader(f)
        # Read the CSV file into a list of dictionaries
        data = list(reader)

    # Convert the list of dictionaries to a dictionary of dictionaries
    # The keys are the 'sample' values and the values are dictionaries of the other columns
    read_metrics = {row["sample"]: {k: v for k, v in row.items() if k != "sample"} for row in data}

    # Read the CSV file
    cov_all = pd.read_csv(coverage_all, header=0)

    # Convert string arrays to actual arrays
    for column in cov_all.columns:
        cov_all[column] = cov_all[column].apply(ast.literal_eval)

    # Read the CSV file
    cov_unique = pd.read_csv(coverage_unique)

    # Convert string arrays to actual arrays
    for column in cov_unique.columns:
        cov_unique[column] = cov_unique[column].apply(ast.literal_eval)


    return read_metrics, cov_all, cov_unique

def plot_coverage(cov_all, cov_unique, gtf_path_car):
    """_summary_

    Args:
        cov_all (_type_): _description_
        cov_unique (_type_): _description_
        gtf_path_car (_type_): _description_
    """
    samples = cov_all.columns.tolist()

    coverage_plot(samples, gtf_path_car, cov_all, cov_unique)


def coverage_plot(samples, gtf_path_car, cov_all, cov_unique):
    """
    Generate a coverage plot for the given samples and CAR construct.
    """
    plt.rcParams.update({'font.size': 14})

    #get max_end from gtf_path_car
    car_start_end = get_car_start_end_from_gtfs(gtf_path_car)
    car = car_start_end[0]

    max_end = car_start_end[2]

    cov_all_dict = {col: cov_all[col][0] for col in cov_all.columns}
    cov_unique_dict = {col: cov_unique[col][0] for col in cov_unique.columns}

    # Calculate the global maximum coverage across all samples
    global_max_coverage = max(
        max(max(values) for values in cov_all_dict.values()),
        max(max(values) for values in cov_unique_dict.values())
    )

    #Use pyGenomeWiz to plot the coverage profile
    gv = GenomeViz(tick_style="axis", fig_track_height=0.8)
    track = gv.add_feature_track(car, max_end)

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

    color_index = 0
    #get start stop from each line in gtf_path_car
    #print(f"GTF of CAR construct used at: {gtf_path_car}")
    with open(gtf_path_car, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 9:
                logging.error(
                    "Error: The GTF file at %s is not formatted correctly.",
                    gtf_path_car)
            if car in line:
                if line.strip().split("\t")[2] =="transcript":
                    continue
                start = int(line.strip().split("\t")[3])
                end = int(line.strip().split("\t")[4])

                attributes = line.strip().split("\t")[-1].split(";")
                parsed_attrs = {
                    k.strip(): v.strip().strip('"')
                    for attr in attributes if attr
                    for k, v in [attr.strip().split(" ", 1)]
                }
                exon_id = parsed_attrs.get("exon_id", "Exon ID not found")
                track.add_feature(
                    start, end,
                    label=exon_id,
                    facecolor=colors[color_index % len(colors)],
                    edgecolor="black",
                    linewidth=0.5
                )
                color_index += 1

    for sample in samples:
        gv.get_track(car).add_subtrack(name=sample+"all_reads")
        gv.get_track(car).add_subtrack(name=sample+"unique_reads")


    fig = gv.plotfig()

    for sample in samples:
        coverage = cov_all_dict[sample]
        coverage_unique = cov_unique_dict[sample]

        pos_list = np.arange(1, max_end + 1)
        subtrack = gv.get_track(car).get_subtrack(sample+"all_reads")
        subtrack.ax.fill_between(pos_list, coverage, alpha=0.6, color="#7e4794")

        subtrack_2 = gv.get_track(car).get_subtrack(sample+"unique_reads")
        subtrack_2.ax.fill_between(pos_list, coverage_unique, alpha=0.6, color="#59a89c")

        subtrack.ax.set_ylim(0, global_max_coverage + 10)
        subtrack_2.ax.set_ylim(0, global_max_coverage + 10)

        subtrack.ax.text(
            gv.top_track.offset,
            global_max_coverage / 2,
            sample,
            ha="right",
            va="center"
        )

    handles = [
        Patch(color="#7e4794", alpha=0.6, label="All reads"),
        Patch(color="#59a89c",alpha=0.6, label="Uniquely mapped reads"),
    ]
    fig.legend(
        handles=handles,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05),
        ncol=2,
        fontsize=14
    )
    # fig.savefig("coverage_plot.png", bbox_inches='tight')

def get_car_start_end_from_gtfs(gtf_path_car: str) -> Tuple[str,int,int]:
    """Get the CAR constructs from the GTF files.

    Args:
        gtf_paths_car (List[str]): List of paths to the GTF files.

    Raises:
        ValueError: If the CAR region is not identified in the GTF files.

    Returns:
        List[str,int,int]: List of CAR constructs with the minimum start and maximum end positions.
    """

    start = float("inf")
    end = float("-inf")
    car = ""

    with open(gtf_path_car, encoding="utf-8") as f:
        lines = f.readlines()
        first_line = lines[0].strip().split("\t")
        last_line = lines[-1].strip().split("\t")

        car = str(first_line[0])
        start = int(first_line[3])
        end = int(last_line[4])

    return car, start, end

def plot_read_metrices(read_metrics):
    """
    Generate a bar plot for absolute reads and annotate it with the number of absolute reads.

    Args:
        read_metrics (dict): A dictionary where keys are sample names and values contain 
            read count metrics such as 'absolute_reads', 'forward_reads', and 'reverse_reads'.
    """
    # Check if all values in read_metrics are zero
    if all(int(metrics['absolute_reads']) == 0 for metrics in read_metrics.values()):
        logging.info("All absolute read values are zero. Skipping plot generation.")
        
    # Extract unique samples
    samples = list(read_metrics.keys())
    bar_plot(read_metrics, samples)

def bar_plot(read_metrics, samples):
    """
    Generate a bar plot for absolute reads and annotate the plot with the number of absolute reads.

    Args:
        read_metrics (dict): Read metrics dictionary (absolute, forward, and reverse reads).
        samples (list): List of sample names to include in the plot.
    """
    # Filter the read metrics for the given samples
    read_metrics = {
        sample: metrics
        for sample, metrics in read_metrics.items()
        if sample in samples
    }

    # Extract absolute reads
    absolute_reads = [int(value['absolute_reads']) for value in read_metrics.values()]

    # Create the plot
    _, ax = plt.subplots(figsize=(20, 10))

    # Bar plot for absolute reads
    bars = ax.bar(samples, absolute_reads, color='#7e4794', label='Absolute Reads')

    # Annotate the bars with the number of absolute reads
    for bar_patch, value in zip(bars, absolute_reads):
        ax.text(
            bar_patch.get_x() + bar_patch.get_width() / 2,  # X position
            bar_patch.get_height() + 5,  # Y position (slightly above the bar)
            str(value),  # Text to display
            ha='center', va='bottom', fontsize=20  # Alignment and font size
        )

    # Formatting the plot
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, fontsize=22)
    ax.set_xlabel('Sample', fontsize=22)
    ax.set_ylabel('Reads (abs.)', fontsize=22)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Save the plot
    plt.tight_layout()
    # plt.savefig(f"Read_metrics_{car}.png", bbox_inches='tight')

def main():
    """
    Main function
    """
    args = parse_arguments()
    read_metrics, cov_all, cov_unique = read_data(\
        args.read_metrics, args.coverage_all, args.coverage_unique)
    plot_coverage(cov_all, cov_unique, args.CAR_gtf_file)
    plot_read_metrices(read_metrics)

if __name__ == "__main__":
    main()
