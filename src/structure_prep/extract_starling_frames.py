"""
Extract individual PDB frames from Starling's .xtc trajectory output
and optionally pre-filter by predicted disorder score.

Starling outputs:
    {name}_STARLING.pdb   - topology file
    {name}_STARLING.xtc   - trajectory (N conformers)

This script extracts every Nth frame as a PDB file for downstream
FoldMason analysis.

Usage:
    python extract_starling_frames.py \
        --topology tau_monomer_STARLING.pdb \
        --trajectory tau_monomer_STARLING.xtc \
        --output data/structures/tau_monomer_ensemble/ \
        --n_frames 100 \
        --stride 4
"""

import argparse
import os


def extract_frames(
    topology: str,
    trajectory: str,
    output_dir: str,
    n_frames: int = 100,
    stride: int = None,
):
    """
    Extract frames from Starling XTC trajectory to individual PDB files.

    topology:    _STARLING.pdb file (Starling topology)
    trajectory:  _STARLING.xtc file (Starling trajectory)
    n_frames:    how many frames to extract total
    stride:      take every Nth frame (computed from n_frames if not given)
    """
    try:
        import MDAnalysis as mda
    except ImportError:
        raise ImportError(
            "MDAnalysis required: pip install MDAnalysis"
        )

    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading trajectory: {trajectory}")
    u = mda.Universe(topology, trajectory)
    total_frames = len(u.trajectory)
    print(f"  Total frames in trajectory: {total_frames}")

    if stride is None:
        stride = max(1, total_frames // n_frames)

    frames_to_write = range(0, total_frames, stride)
    print(f"  Extracting every {stride}th frame → {len(list(frames_to_write))} PDBs")

    protein = u.select_atoms("protein")
    name_stem = os.path.splitext(os.path.basename(topology))[0].replace("_STARLING", "")

    written = []
    for i, ts in enumerate(u.trajectory[::stride]):
        if len(written) >= n_frames:
            break
        out_path = os.path.join(output_dir, f"{name_stem}_frame{i:04d}.pdb")
        protein.write(out_path)
        written.append(out_path)

    print(f"  Wrote {len(written)} PDB frames to {output_dir}")
    return written


def predict_disorder_boundaries(sequence: str, window: int = 9) -> list[tuple[int, int]]:
    """
    Simple hydrophobicity + charge-based disorder prediction to suggest trim points.
    Uses Uversky's charge-hydropathy plot heuristic.

    Returns list of (start, end) tuples for predicted disordered regions.
    This is a rough heuristic — confirm with IUPred2A or ESMFold pLDDT.
    """
    # Kyte-Doolittle hydrophobicity scale
    kd = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    }

    n = len(sequence)
    half = window // 2
    scores = []

    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        window_seq = sequence[start:end]
        h = sum(kd.get(aa, 0) for aa in window_seq) / len(window_seq)
        # Normalize to [0, 1] (roughly: -4.5 to 4.5 range)
        h_norm = (h + 4.5) / 9.0
        scores.append(h_norm)

    # Disordered = low hydrophobicity (score < 0.45 roughly)
    threshold = 0.45
    disordered_regions = []
    in_disorder = False
    start_idx = 0

    for i, score in enumerate(scores):
        if score < threshold and not in_disorder:
            in_disorder = True
            start_idx = i
        elif score >= threshold and in_disorder:
            in_disorder = False
            if i - start_idx >= 10:  # min 10 aa disordered stretch
                disordered_regions.append((start_idx + 1, i))  # 1-indexed

    if in_disorder and n - start_idx >= 10:
        disordered_regions.append((start_idx + 1, n))

    return disordered_regions


def main():
    parser = argparse.ArgumentParser(
        description="Extract PDB frames from Starling XTC trajectory"
    )
    parser.add_argument("--topology", required=True,
                        help="Starling topology PDB (_STARLING.pdb)")
    parser.add_argument("--trajectory", required=True,
                        help="Starling trajectory XTC (_STARLING.xtc)")
    parser.add_argument("--output", required=True,
                        help="Output directory for extracted PDB frames")
    parser.add_argument("--n_frames", type=int, default=100,
                        help="Number of frames to extract (default: 100)")
    parser.add_argument("--stride", type=int, default=None,
                        help="Extract every Nth frame (default: auto from n_frames)")
    parser.add_argument("--sequence", default=None,
                        help="Optional: amino acid sequence for disorder boundary prediction")
    args = parser.parse_args()

    # Optionally show disorder boundary suggestion
    if args.sequence:
        regions = predict_disorder_boundaries(args.sequence)
        print("Predicted disordered regions (rough heuristic — verify with IUPred2A):")
        for start, end in regions:
            print(f"  Residues {start}–{end} ({end - start + 1} aa)")
        print("  Recommendation: trim your input to the largest disordered region before running Starling")
        print()

    extract_frames(
        topology=args.topology,
        trajectory=args.trajectory,
        output_dir=args.output,
        n_frames=args.n_frames,
        stride=args.stride,
    )


if __name__ == "__main__":
    main()
