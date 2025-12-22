#!/usr/bin/env python3
import sys
import os
from collections import defaultdict

def read_fasta(path):
    """Simple FASTA reader that returns {name: sequence}."""
    seqs = {}
    name = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq_chunks).upper()
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name:
            seqs[name] = "".join(seq_chunks).upper()
    return seqs


def get_kmers(seq, k):
    """Return list of all k-mers in a sequence."""
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]


def main():
    if len(sys.argv) < 4:
        print("Usage: python unique_kmers.py transgenes.fasta [k] output_dir")
        sys.exit(1)

    fasta = sys.argv[1]
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 31
    outdir = sys.argv[3]
    os.makedirs(outdir, exist_ok=True)

    print(f"Loading FASTA: {fasta}")
    seqs = read_fasta(fasta)
    names = list(seqs.keys())
    
    # store kmers per transgene
    kmers = {name: set(get_kmers(seqs[name], k)) for name in names}

    # compute unique kmers per transgene
    unique = {}
    for name in names:
        others = set()
        for other in names:
            if other != name:
                others |= kmers[other]
        unique[name] = kmers[name] - others

    # compute shared counts matrix
    shared_matrix = defaultdict(dict)
    for a in names:
        for b in names:
            shared_matrix[a][b] = len(kmers[a] & kmers[b])

    # print summary
    print("\n=== UNIQUE k-mers per transgene ===")
    for name in names:
        print(f"{name}: {len(unique[name])} unique {k}-mers")

    print("\n=== SHARED k-mer matrix ===")
    header = "\t" + "\t".join(names)
    print(header)
    for a in names:
        row = [a] + [str(shared_matrix[a][b]) for b in names]
        print("\t".join(row))

    # optional: write unique kmers to files
    for name in names:
        out = os.path.join(outdir, f"{name}_unique_{k}mers.txt")
        with open(out, "w") as f:
            for km in unique[name]:
                f.write(km + "\n")
        print(f"Written unique kmers for {name} → {out}")

if __name__ == "__main__":
    main()
