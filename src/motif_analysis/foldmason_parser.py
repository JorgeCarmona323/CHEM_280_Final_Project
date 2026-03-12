"""
Parse FoldMason structural MSA output to extract conserved structural motifs.

FoldMason outputs a FASTA-like structural alignment where each position
encodes secondary structure and residue identity.

Usage:
    python foldmason_parser.py \
        --msa results/foldmason/all_tau/msa.fasta \
        --pdb_dir results/esmfold/ \
        --output results/motifs/foldmason_motifs.csv \
        --conservation_threshold 0.7
"""

import argparse
import os
import csv
from collections import Counter


SS_MAP = {
    "H": "helix",
    "E": "strand",
    "C": "coil",
    "T": "turn",
    "-": "gap",
}


def read_structural_msa(msa_path: str) -> dict[str, str]:
    """Read FoldMason structural alignment FASTA."""
    records = {}
    header, lines = None, []
    with open(msa_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(lines)
                header = line[1:].split()[0]
                lines = []
            else:
                lines.append(line)
    if header is not None:
        records[header] = "".join(lines)
    return records


def compute_conservation(msa: dict[str, str]) -> list[dict]:
    """
    Compute per-column conservation statistics from structural MSA.

    Returns list of dicts with column index, dominant character, conservation score.
    """
    sequences = list(msa.values())
    if not sequences:
        return []

    n_cols = len(sequences[0])
    columns = []

    for col_idx in range(n_cols):
        col_chars = [seq[col_idx] for seq in sequences if col_idx < len(seq)]
        non_gap = [c for c in col_chars if c != "-"]

        if not non_gap:
            columns.append({
                "col": col_idx,
                "dominant": "-",
                "conservation": 0.0,
                "gap_fraction": 1.0,
                "ss_type": "gap",
            })
            continue

        counter = Counter(non_gap)
        dominant_char, dominant_count = counter.most_common(1)[0]
        conservation = dominant_count / len(col_chars)
        gap_fraction = (len(col_chars) - len(non_gap)) / len(col_chars)

        columns.append({
            "col": col_idx,
            "dominant": dominant_char,
            "conservation": conservation,
            "gap_fraction": gap_fraction,
            "ss_type": SS_MAP.get(dominant_char.upper(), "unknown"),
        })

    return columns


def find_conserved_motifs(
    columns: list[dict],
    conservation_threshold: float = 0.7,
    min_motif_length: int = 4,
) -> list[dict]:
    """
    Identify stretches of conserved columns as structural motifs.

    Returns list of motifs with start, end, consensus pattern, avg conservation.
    """
    motifs = []
    in_motif = False
    motif_start = 0
    motif_cols = []

    for col in columns:
        if col["conservation"] >= conservation_threshold and col["ss_type"] != "gap":
            if not in_motif:
                in_motif = True
                motif_start = col["col"]
                motif_cols = []
            motif_cols.append(col)
        else:
            if in_motif and len(motif_cols) >= min_motif_length:
                motifs.append({
                    "start_col": motif_start,
                    "end_col": motif_cols[-1]["col"],
                    "length": len(motif_cols),
                    "consensus": "".join(c["dominant"] for c in motif_cols),
                    "avg_conservation": sum(c["conservation"] for c in motif_cols) / len(motif_cols),
                    "dominant_ss": Counter(c["ss_type"] for c in motif_cols).most_common(1)[0][0],
                })
            in_motif = False
            motif_cols = []

    # Handle motif running to end
    if in_motif and len(motif_cols) >= min_motif_length:
        motifs.append({
            "start_col": motif_start,
            "end_col": motif_cols[-1]["col"],
            "length": len(motif_cols),
            "consensus": "".join(c["dominant"] for c in motif_cols),
            "avg_conservation": sum(c["conservation"] for c in motif_cols) / len(motif_cols),
            "dominant_ss": Counter(c["ss_type"] for c in motif_cols).most_common(1)[0][0],
        })

    return motifs


def main():
    parser = argparse.ArgumentParser(description="Parse FoldMason structural MSA for conserved motifs")
    parser.add_argument("--msa", required=True, help="FoldMason structural alignment FASTA")
    parser.add_argument("--pdb_dir", default=None, help="Directory of source PDBs (for annotation)")
    parser.add_argument("--output", required=True, help="Output CSV for motifs")
    parser.add_argument("--conservation_threshold", type=float, default=0.7)
    parser.add_argument("--min_motif_length", type=int, default=4)
    args = parser.parse_args()

    print(f"Reading structural MSA from {args.msa}")
    msa = read_structural_msa(args.msa)
    print(f"  {len(msa)} sequences, alignment length {len(next(iter(msa.values())))}")

    columns = compute_conservation(msa)
    motifs = find_conserved_motifs(
        columns,
        conservation_threshold=args.conservation_threshold,
        min_motif_length=args.min_motif_length,
    )

    print(f"Found {len(motifs)} conserved structural motifs "
          f"(conservation >= {args.conservation_threshold}, length >= {args.min_motif_length})")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "start_col", "end_col", "length", "consensus",
            "avg_conservation", "dominant_ss"
        ])
        writer.writeheader()
        writer.writerows(motifs)

    print(f"Motifs written to {args.output}")

    # Print top motifs
    print("\nTop 10 conserved motifs:")
    for m in sorted(motifs, key=lambda x: -x["avg_conservation"])[:10]:
        print(f"  cols {m['start_col']}-{m['end_col']} | "
              f"len={m['length']} | "
              f"ss={m['dominant_ss']} | "
              f"cons={m['avg_conservation']:.2f} | "
              f"pattern={m['consensus'][:20]}{'...' if len(m['consensus']) > 20 else ''}")


if __name__ == "__main__":
    main()
