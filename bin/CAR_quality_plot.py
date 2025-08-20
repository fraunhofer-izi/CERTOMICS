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
import math
import plotly.graph_objs as go
from plotly.subplots import make_subplots

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

    interactive_coverage_plot(samples, cov_all, cov_unique, gtf_path_car)
    

def interactive_coverage_plot(samples, cov_all, cov_unique, gtf_path_car, output_html="interactive_coverage.html"):
    # === Parse GTF for exon features ===
    exon_features = []
    max_end = 0
    with open(gtf_path_car, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 9 or parts[2] == "transcript":
                continue
            start = int(parts[3])
            end = int(parts[4])
            max_end = max(max_end, end)
            attributes = {
                k.strip(): v.strip().strip('"')
                for attr in parts[8].split(";") if attr
                for k, v in [attr.strip().split(" ", 1)]
            }
            exon_id = attributes.get("exon_id", "Exon")
            exon_features.append((start, end, exon_id))

    # === Build coverage data ===
    cov_all_dict = {col: cov_all[col][0] for col in cov_all.columns}
    cov_unique_dict = {col: cov_unique[col][0] for col in cov_unique.columns}

    row_heights = [0.15]  # GTF feature track at top
    row_names = ["Features"]

    for sample in samples:
        max_val = max(max(cov_all_dict[sample]), max(cov_unique_dict[sample]))
        row_heights.append(max_val)
        row_names.append(sample)

    # Normalize row heights
    log_heights = [math.log10(h + 1) for h in row_heights]  # +1 to avoid log(0)
    norm_heights = [max(h / max(log_heights), 0.1) for h in log_heights]

    fig = make_subplots(
        rows=len(row_names),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.01,
        row_heights=norm_heights
    )

    pos_list = np.arange(1, max_end + 1)

    fig.update_yaxes(
        range=[0, 1],              # restrict y-axis to exact arrow height
        showticklabels=False,
        showgrid=False,
        fixedrange=True,           # prevent zoom from stretching it again
        row=1, col=1
    )

    y0 = 0
    y1 = 1  
    y_center = 0.5
    arrow_tip_length = 10

    colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
        "#bcbd22", "#17becf"
    ]

    for i, (start, end, exon_id) in enumerate(exon_features):
        # Avoid negative arrow base if exon too short
        arrow_base = max(start, end - arrow_tip_length)

        path = (
            f"M {start},{y0} "
            f"L {arrow_base},{y0} "
            f"L {end},{y_center} "
            f"L {arrow_base},{y1} "
            f"L {start},{y1} "
            f"Z"
        )

        fig.add_shape(
            type="path",
            path=path,
            xref="x",
            yref="y1",
            fillcolor=colors[i % len(colors)],
            line=dict(color="black"),
            layer="below"
        )

        fig.add_annotation(
            x=(start + end) // 2,
            y=0.99,
            xref='x',
            yref='y1',
            text=exon_id,
            showarrow=True,
            arrowhead=2,
            ax=0,
            ay=-40,
            font=dict(size=13),
            textangle=-45
        )

    # === Add coverage plots ===
    for i, sample in enumerate(samples):
        all_y = cov_all_dict[sample]
        unique_y = cov_unique_dict[sample]
        row = i + 2

        # Add all reads first
        fig.add_trace(go.Scatter(
            x=pos_list,
            y=all_y,
            name="All reads",
            fill='tozeroy',
            mode='lines',
            line=dict(color="rgba(126, 71, 148, 0.8)", width=2.5),  # purple, transparent
            showlegend=(i == 0),
        ), row=row, col=1)

        # Then uniquely mapped reads — this one will appear on top
        fig.add_trace(go.Scatter(
            x=pos_list,
            y=unique_y,
            name="Unique reads",
            fill='tozeroy',
            mode='lines',
            line=dict(color="rgba(89,168,156,0.3)"),
            showlegend=(i == 0),
        ), row=row, col=1)

        # Add sample name as annotation (left of tracks)
        fig.add_annotation(
            xref="paper",
            yref=f'y{row} domain',
            x=-0.055,
            y=0.5,
            text=sample,
            textangle=-90,
            showarrow=False,
            font=dict(size=13, color="black"),
            xanchor="center"
        )

    fig.update_layout(
        height= 100 + 120 * len(samples),
        yaxis_title=None,
        legend=dict(
            orientation="h",
            yanchor="top",
            y=1.30,        # aligned with bottom edge of plot
            xanchor="center",
            x=0.5
        ),
        margin=dict(t=120, b=20, l=50, r=10)
    )
    fig.update_yaxes(showgrid=False)
    fig.update_xaxes(showgrid=False)
    # Add x-axis title only to bottom subplot
    num_rows = len(samples) + 1
    fig.update_xaxes(title_text="Genomic Position", row=num_rows, col=1)
    fig.write_html("coverage_plot.html")  # saved to current working directory
    return fig


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
            ha='center', va='bottom', fontsize=22  # Alignment and font size
        )

    # Formatting the plot
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, fontsize=24)
    ax.set_xlabel('Sample', fontsize=24)
    ax.set_ylabel('Reads (abs.)', fontsize=24)
    ax.tick_params(axis='y', labelsize=24)
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
