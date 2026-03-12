"""
Stitch IDP ensemble structures onto stable domain scaffolds.

Usage:
    python stitch_constructs.py \
        --stable data/structures/tau_af2/tau_full.pdb \
        --idp_ensemble data/structures/tau_monomer_ensemble/ \
        --idp_residues 244-368 \
        --output data/structures/tau_stitched/
"""

import argparse
import os
import glob
import numpy as np
from pathlib import Path


def parse_pdb_ca(pdb_path: str) -> dict:
    """Extract Cα coordinates and residue info from a PDB file."""
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                chain = line[21]
                resnum = int(line[22:26].strip())
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                residues[(chain, resnum)] = np.array([x, y, z])
    return residues


def superpose_junction(
    stable_ca: dict,
    idp_ca: dict,
    junction_residues: list,
    chain_stable: str = "A",
    chain_idp: str = "A",
) -> np.ndarray:
    """
    Compute rotation+translation to superpose IDP onto stable domain at junction.

    junction_residues: list of residue numbers that are shared (±5 of the IDP start).
    Returns (R, t) such that idp_coords_aligned = R @ idp_coords + t
    """
    stable_pts = np.array([stable_ca[(chain_stable, r)] for r in junction_residues
                           if (chain_stable, r) in stable_ca])
    idp_pts = np.array([idp_ca[(chain_idp, r)] for r in junction_residues
                        if (chain_idp, r) in idp_ca])

    if len(stable_pts) < 3 or len(idp_pts) < 3:
        raise ValueError(f"Not enough junction residues found for superposition. "
                         f"Stable: {len(stable_pts)}, IDP: {len(idp_pts)}")

    n = min(len(stable_pts), len(idp_pts))
    stable_pts, idp_pts = stable_pts[:n], idp_pts[:n]

    # Center
    stable_center = stable_pts.mean(axis=0)
    idp_center = idp_pts.mean(axis=0)
    A = stable_pts - stable_center
    B = idp_pts - idp_center

    # SVD-based rotation
    H = B.T @ A
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = stable_center - R @ idp_center

    return R, t


def stitch_pdb(
    stable_pdb: str,
    idp_pdb: str,
    idp_start: int,
    idp_end: int,
    output_pdb: str,
    junction_window: int = 5,
):
    """
    Replace residues idp_start..idp_end in stable_pdb with the IDP structure,
    superposing at the junction.
    """
    stable_ca = parse_pdb_ca(stable_pdb)
    idp_ca = parse_pdb_ca(idp_pdb)

    # Junction: ±junction_window residues around IDP start
    junction = list(range(idp_start - junction_window, idp_start + junction_window + 1))

    R, t = superpose_junction(stable_ca, idp_ca, junction)

    # Read all atoms, transform IDP region, write stitched PDB
    stable_atoms = []
    idp_atoms = []

    with open(stable_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resnum = int(line[22:26].strip())
                if resnum < idp_start or resnum > idp_end:
                    stable_atoms.append(line)

    with open(idp_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resnum = int(line[22:26].strip())
                if idp_start <= resnum <= idp_end:
                    # Transform coordinates
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    coords = R @ np.array([x, y, z]) + t
                    new_line = (line[:30] +
                                f"{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}" +
                                line[54:])
                    idp_atoms.append(new_line)

    os.makedirs(os.path.dirname(output_pdb), exist_ok=True)
    with open(output_pdb, "w") as f:
        # Write stable (non-IDP) atoms
        for line in sorted(stable_atoms, key=lambda l: int(l[22:26].strip())):
            f.write(line)
        # Write IDP atoms (transformed)
        for line in sorted(idp_atoms, key=lambda l: int(l[22:26].strip())):
            f.write(line)
        f.write("END\n")

    print(f"Stitched structure written to {output_pdb}")


def main():
    parser = argparse.ArgumentParser(description="Stitch IDP ensemble onto stable domain")
    parser.add_argument("--stable", required=True, help="Stable domain PDB (e.g. AF2 output)")
    parser.add_argument("--idp_ensemble", required=True,
                        help="Directory of IDP Starling ensemble PDBs, or single PDB")
    parser.add_argument("--idp_residues", required=True,
                        help="Residue range of IDP region, e.g. 244-368")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--max_structures", type=int, default=5,
                        help="Max ensemble members to stitch (default: 5)")
    args = parser.parse_args()

    idp_start, idp_end = map(int, args.idp_residues.split("-"))

    # Collect IDP structures
    if os.path.isdir(args.idp_ensemble):
        idp_pdbs = sorted(glob.glob(os.path.join(args.idp_ensemble, "*.pdb")))[:args.max_structures]
    else:
        idp_pdbs = [args.idp_ensemble]

    if not idp_pdbs:
        raise FileNotFoundError(f"No PDB files found in {args.idp_ensemble}")

    os.makedirs(args.output, exist_ok=True)

    for i, idp_pdb in enumerate(idp_pdbs):
        name = Path(idp_pdb).stem
        output_pdb = os.path.join(args.output, f"stitched_{name}.pdb")
        try:
            stitch_pdb(
                stable_pdb=args.stable,
                idp_pdb=idp_pdb,
                idp_start=idp_start,
                idp_end=idp_end,
                output_pdb=output_pdb,
            )
        except Exception as e:
            print(f"Warning: failed to stitch {idp_pdb}: {e}")

    print(f"\nStitching complete. {len(idp_pdbs)} structures written to {args.output}")


if __name__ == "__main__":
    main()
