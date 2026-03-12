"""
Use FoldMason structural alignment to filter and rank an IDP ensemble,
selecting the most structurally representative conformers before stitching.

Workflow:
    1. Run FoldMason on the full Starling ensemble
    2. Parse the structural MSA to compute per-conformer deviation from consensus
    3. Rank conformers by structural consistency score
    4. Output top-N representative conformers for downstream stitching

Usage:
    python foldmason_refine.py \
        --ensemble data/structures/tau_monomer_ensemble/ \
        --output data/structures/tau_monomer_refined/ \
        --top_n 5 \
        --tmp_dir tmp/foldmason_tau/
"""

import argparse
import os
import glob
import shutil
import subprocess
import csv
from pathlib import Path
import numpy as np


def run_foldmason(
    pdb_dir: str,
    output_prefix: str,
    tmp_dir: str,
) -> tuple[str, str]:
    """
    Run FoldMason easy-msa on all PDBs in pdb_dir.

    Returns (msa_fasta_path, lddt_report_path)
    """
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(os.path.dirname(output_prefix) or ".", exist_ok=True)

    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in {pdb_dir}")

    print(f"Running FoldMason on {len(pdb_files)} structures...")

    # Write a file list for FoldMason
    filelist = os.path.join(tmp_dir, "filelist.txt")
    with open(filelist, "w") as f:
        f.write("\n".join(pdb_files) + "\n")

    cmd = [
        "foldmason", "easy-msa",
        *pdb_files,
        output_prefix,
        tmp_dir,
        "--report-mode", "1",    # produces HTML with per-structure lDDT
        "--threads", "4",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"FoldMason failed:\n{result.stderr}")

    # Real FoldMason outputs: {prefix}_aa.fa, {prefix}_3di.fa, {prefix}.nw, {prefix}.html
    msa_path = output_prefix + "_aa.fa"
    html_path = output_prefix + ".html"   # contains per-structure lDDT (no separate TSV)

    if not os.path.exists(msa_path):
        raise FileNotFoundError(
            f"FoldMason MSA output not found: {msa_path}\n"
            f"Expected outputs: {output_prefix}_aa.fa, {output_prefix}_3di.fa, "
            f"{output_prefix}.nw, {output_prefix}.html"
        )

    print(f"FoldMason MSA written to {msa_path}")
    if os.path.exists(html_path):
        print(f"FoldMason HTML report (with lDDT): {html_path}")

    return msa_path, html_path


def parse_foldmason_lddt(html_path: str) -> dict[str, float]:
    """
    Parse per-structure lDDT scores from FoldMason HTML report.

    FoldMason does not produce a separate TSV — lDDT data is embedded in
    the HTML report (--report-mode 1). This function extracts it via regex.

    Returns dict: {structure_name -> mean_lddt}
    Falls back to empty dict (MSA-only scoring) if HTML is not parseable.
    """
    import re

    scores = {}
    if not os.path.exists(html_path):
        print(f"  Warning: FoldMason HTML report not found at {html_path}, "
              f"falling back to MSA-consistency-only scoring")
        return scores

    with open(html_path) as f:
        content = f.read()

    # FoldMason embeds lDDT as JSON data in the HTML — pattern may vary by version
    # Try to extract: "name": "...", "lddt": 0.xx patterns
    pattern = r'"name"\s*:\s*"([^"]+)"[^}]*?"lddt"\s*:\s*([\d.]+)'
    for m in re.finditer(pattern, content):
        name, lddt = m.group(1), float(m.group(2))
        scores[name] = lddt

    if scores:
        print(f"  Parsed lDDT scores for {len(scores)} structures from HTML report")
    else:
        print(f"  Could not parse lDDT from HTML (format may differ by FoldMason version); "
              f"using MSA-consistency-only scoring")

    return scores


def compute_msa_consistency(msa_path: str) -> dict[str, float]:
    """
    Compute per-conformer structural consistency from FoldMason MSA.

    Consistency score = fraction of non-gap positions where the conformer
    matches the column-wise dominant character (structural consensus).

    Higher score = more representative of the ensemble consensus.
    """
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

    if not records:
        return {}

    sequences = list(records.values())
    names = list(records.keys())
    n_cols = len(sequences[0])

    # Build consensus per column
    from collections import Counter
    consensus = []
    for col in range(n_cols):
        chars = [s[col] for s in sequences if col < len(s)]
        non_gap = [c for c in chars if c != "-"]
        if non_gap:
            consensus.append(Counter(non_gap).most_common(1)[0][0])
        else:
            consensus.append("-")

    # Score each conformer: fraction of non-gap positions matching consensus
    consistency = {}
    for name, seq in zip(names, sequences):
        matches = 0
        total = 0
        for col, (char, cons) in enumerate(zip(seq, consensus)):
            if char != "-" and cons != "-":
                total += 1
                if char == cons:
                    matches += 1
        consistency[name] = matches / total if total > 0 else 0.0

    return consistency


def rank_conformers(
    consistency_scores: dict[str, float],
    lddt_scores: dict[str, float],
    lddt_weight: float = 0.4,
    consistency_weight: float = 0.6,
) -> list[tuple[str, float]]:
    """
    Rank conformers by combined score: consistency + lDDT.

    If lDDT scores are unavailable, rank by consistency alone.
    """
    all_names = list(consistency_scores.keys())

    if lddt_scores:
        # Normalize lDDT to [0, 1]
        lddt_vals = np.array([lddt_scores.get(n, 0.0) for n in all_names])
        lddt_min, lddt_max = lddt_vals.min(), lddt_vals.max()
        if lddt_max > lddt_min:
            lddt_norm = (lddt_vals - lddt_min) / (lddt_max - lddt_min)
        else:
            lddt_norm = np.ones(len(all_names))
    else:
        lddt_norm = np.ones(len(all_names))
        lddt_weight = 0.0
        consistency_weight = 1.0

    cons_vals = np.array([consistency_scores[n] for n in all_names])

    combined = consistency_weight * cons_vals + lddt_weight * lddt_norm
    ranked = sorted(zip(all_names, combined.tolist()), key=lambda x: -x[1])

    return ranked


def copy_top_conformers(
    ranked: list[tuple[str, float]],
    ensemble_dir: str,
    output_dir: str,
    top_n: int = 5,
) -> list[str]:
    """
    Copy top-N ranked conformers to output directory.
    Returns list of output PDB paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    copied = []

    for rank, (name, score) in enumerate(ranked[:top_n]):
        # Find the source PDB — name may be filename stem
        candidates = glob.glob(os.path.join(ensemble_dir, f"{name}.pdb"))
        if not candidates:
            # Try partial match
            candidates = glob.glob(os.path.join(ensemble_dir, f"*{name}*.pdb"))
        if not candidates:
            print(f"  Warning: could not find PDB for {name}, skipping")
            continue

        src = candidates[0]
        dst = os.path.join(output_dir, f"rank{rank+1:02d}_{Path(src).name}")
        shutil.copy2(src, dst)
        copied.append(dst)
        print(f"  Rank {rank+1}: {name} (score={score:.3f}) → {dst}")

    return copied


def main():
    parser = argparse.ArgumentParser(
        description="FoldMason-based ensemble filtering for IDP conformers"
    )
    parser.add_argument("--ensemble", required=True,
                        help="Directory of Starling ensemble PDBs")
    parser.add_argument("--output", required=True,
                        help="Output directory for refined/ranked conformers")
    parser.add_argument("--top_n", type=int, default=5,
                        help="Number of representative conformers to keep (default: 5)")
    parser.add_argument("--tmp_dir", default="tmp/foldmason_refine",
                        help="Temporary directory for FoldMason files")
    parser.add_argument("--lddt_weight", type=float, default=0.4,
                        help="Weight for lDDT in combined score (default: 0.4)")
    parser.add_argument("--consistency_weight", type=float, default=0.6,
                        help="Weight for MSA consistency in combined score (default: 0.6)")
    parser.add_argument("--report", default=None,
                        help="Optional CSV path to write full ranking report")
    args = parser.parse_args()

    output_prefix = os.path.join(args.tmp_dir, "msa")

    # Step 1: Run FoldMason
    msa_path, lddt_path = run_foldmason(
        pdb_dir=args.ensemble,
        output_prefix=output_prefix,
        tmp_dir=args.tmp_dir,
    )

    # Step 2: Parse scores
    lddt_scores = parse_foldmason_lddt(lddt_path)
    consistency_scores = compute_msa_consistency(msa_path)

    print(f"\nScored {len(consistency_scores)} conformers")
    if lddt_scores:
        print(f"lDDT scores available for {len(lddt_scores)} conformers")
    else:
        print("lDDT scores not available — ranking by MSA consistency only")

    # Step 3: Rank
    ranked = rank_conformers(
        consistency_scores,
        lddt_scores,
        lddt_weight=args.lddt_weight,
        consistency_weight=args.consistency_weight,
    )

    print(f"\nTop {args.top_n} conformers (by consistency + lDDT):")
    copied = copy_top_conformers(ranked, args.ensemble, args.output, top_n=args.top_n)

    # Step 4: Write report
    report_path = args.report or os.path.join(args.output, "conformer_ranking.csv")
    with open(report_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "name", "combined_score", "consistency", "lddt"])
        for rank, (name, score) in enumerate(ranked):
            writer.writerow([
                rank + 1,
                name,
                f"{score:.4f}",
                f"{consistency_scores.get(name, 0):.4f}",
                f"{lddt_scores.get(name, 'N/A')}",
            ])
    print(f"\nFull ranking report written to {report_path}")
    print(f"\nRefined ensemble ({len(copied)} conformers) ready in: {args.output}")
    print("Next step: pass these to stitch_constructs.py --idp_ensemble")


if __name__ == "__main__":
    main()
