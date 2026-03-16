#!/usr/bin/env python3



# This tool tries to find sequence motifs that are molecularly the same
# regarding the cyclic rotation or reverse complement of the DNA sequence.
# 
# By Shahryar Alavi
# University of Isfahan, Isfahan, Iran
# UCL Institute of Neurology, London, UK
# December 2025



import argparse
from collections import defaultdict

# -----------------------------
# Basic sequence utilities
# -----------------------------

def reverse(seq):
    return seq[::-1]

def complement(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)

def reverse_complement(seq):
    return reverse(complement(seq))

def rotations(seq):
    """All cyclic rotations of a sequence"""
    return {seq[i:] + seq[:i] for i in range(len(seq))}

# -----------------------------
# Motif equivalence class
# -----------------------------

def motif_equivalence_class(motif):
    """
    Generate the full equivalence class of a motif:
    rotations + reverse + reverse complement (+ their rotations)
    """
    eq = set()

    for s in (
        motif,
        reverse(motif),
        complement(motif),
        reverse_complement(motif)
    ):
        eq.update(rotations(s))

    return eq

def canonical_motif(motif):
    """
    Choose a canonical representative:
    lexicographically smallest string in the equivalence class
    """
    return min(motif_equivalence_class(motif))

# -----------------------------
# Main logic
# -----------------------------

def main(bed_file, out_file):

    # motif_length -> canonical_motif -> list of observed motifs
    motif_groups = defaultdict(lambda: defaultdict(list))

    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            motif = fields[3].upper()
            mlen = len(motif)

            canon = canonical_motif(motif)
            motif_groups[mlen][canon].append(motif)

    # Prepare output
    rows = []

    for mlen, groups in motif_groups.items():
        for canon, motifs in groups.items():
            rows.append((
                mlen,
                ",".join(sorted(set(motifs))),
                len(motifs)
            ))

    # Sort by count (descending)
    rows.sort(key=lambda x: x[2], reverse=True)

    # Write TSV
    with open(out_file, "w") as out:
        out.write("Motif_length\tEquivalent_motifs\tSTR_count\n")
        for mlen, motif_list, count in rows:
            out.write(f"{mlen}\t{motif_list}\t{count}\n")

# -----------------------------
# CLI
# -----------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Group permutation STR motifs - by their cyclic rotation, reverse, and reverse-complement equivalence",
        epilog="By Shahryar Alavi - December 2025"
    )
    parser.add_argument(
        "--bed",
        required=True,
        help="Input BED file (with motif sequence in column 4)"
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV file"
    )

    args = parser.parse_args()
    main(args.bed, args.out)
