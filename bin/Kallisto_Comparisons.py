#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import glob

def extract_counts(quant_dirs):
    count_matrix = {}
    for qdir in quant_dirs:
        sample_id = os.path.basename(qdir).replace("quant_", "")
        abundance_path = os.path.join(qdir, "abundance.tsv")

        if os.path.exists(abundance_path):
            df = pd.read_csv(abundance_path, sep='\t')
            count_matrix[sample_id] = df.set_index("target_id")["est_counts"]
        else:
            print(f"{abundance_path} not found.")

    return pd.DataFrame(count_matrix).fillna(0).astype(int).T

def main():
    parser = argparse.ArgumentParser(
        description="Generate est_counts matrix from kallisto quant_* directories."
    )
    parser.add_argument(
        "--input-dir", nargs='+', required=True, help="One or more quant_* directories"
    )  
    parser.add_argument(
        "--output", default="matrix.csv", help="Output CSV file"
    )

    args = parser.parse_args()
    quant_dirs = args.input_dir

    if not quant_dirs:
        print("No quant_* directories found.")
        return

    print(f"Found {len(quant_dirs)} quant directories. Building count matrix...")

    count_df = extract_counts(quant_dirs)
    count_df.to_csv(args.output)

    print(f"Matrix saved to {args.output}")

if __name__ == "__main__":
    main()
