"""
Validate stitched "Frankenstein" constructs using FoldSeek identity confirmation.

Workflow:
    1. Starling generates IDP conformers
    2. FoldMason finds the most structurally consistent IDP motifs across the ensemble
    3. Best motifs are stitched to stable/solved regions → Frankenstein construct
    4. MD relaxation (run_md_relaxation.py)
    5. FoldSeek searches the relaxed construct against PDB/reference DB
       → confirms the structure is recognizable as the intended protein

The FoldSeek step answers: "Does the Frankenstein protein still look like
the protein it's supposed to be?" High TM-score hit to the reference protein
= the IDP conformer is structurally plausible in context.

Usage:
    python validate_construct.py \
        --construct data/structures/tau_relaxed/tau_monomer_relaxed.pdb \
        --reference_name "Tau" \
        --foldseek_db pdb \
        --output results/validation/construct_validation/ \
        --min_tm 0.5
"""

import argparse
import os
import subprocess
import csv
from pathlib import Path


def run_foldseek_search(
    query_pdb: str,
    db: str,
    output_tsv: str,
    tmp_dir: str,
    sensitivity: float = 9.5,
) -> str:
    """
    Run FoldSeek easy-search on a single query structure.

    Returns path to results TSV.
    """
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)

    cmd = [
        "foldseek", "easy-search",
        query_pdb,
        db,
        output_tsv,
        tmp_dir,
        "--format-output", "query,target,evalue,alntmscore,rmsd,prob,taxname,taxid,tlen,qlen",
        "-s", str(sensitivity),
        "--exhaustive-search", "0",
    ]

    print(f"Running FoldSeek: {query_pdb} vs {db}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"FoldSeek failed:\n{result.stderr}")

    print(f"FoldSeek results written to {output_tsv}")
    return output_tsv


def parse_foldseek_hits(tsv_path: str, min_tm: float = 0.5) -> list[dict]:
    """Parse FoldSeek TSV, filter by TM-score threshold."""
    hits = []
    if not os.path.exists(tsv_path):
        return hits

    with open(tsv_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 4:
                continue
            try:
                tm_score = float(row[3])
            except (ValueError, IndexError):
                continue
            if tm_score >= min_tm:
                hits.append({
                    "query": row[0] if len(row) > 0 else "",
                    "target": row[1] if len(row) > 1 else "",
                    "evalue": row[2] if len(row) > 2 else "",
                    "tm_score": tm_score,
                    "rmsd": row[4] if len(row) > 4 else "",
                    "prob": row[5] if len(row) > 5 else "",
                    "taxname": row[6] if len(row) > 6 else "",
                    "taxid": row[7] if len(row) > 7 else "",
                })

    return sorted(hits, key=lambda x: -x["tm_score"])


def confirm_identity(
    hits: list[dict],
    reference_name: str,
    top_k: int = 10,
) -> dict:
    """
    Check if the top FoldSeek hits match the expected protein identity.

    reference_name: string to look for in hit target names / taxonomy
                    (e.g. "Tau", "MAPT", "Homo sapiens")

    Returns a dict with pass/fail + supporting evidence.
    """
    if not hits:
        return {
            "confirmed": False,
            "reason": "No FoldSeek hits above TM-score threshold",
            "top_hit": None,
            "identity_hits": [],
        }

    top_hits = hits[:top_k]
    top_hit = top_hits[0]

    # Check if reference name appears in any top hit
    ref_lower = reference_name.lower()
    identity_hits = [
        h for h in top_hits
        if ref_lower in h["target"].lower() or ref_lower in h["taxname"].lower()
    ]

    confirmed = len(identity_hits) > 0 or top_hit["tm_score"] >= 0.7

    return {
        "confirmed": confirmed,
        "reason": (
            f"Top hit TM-score={top_hit['tm_score']:.3f} to {top_hit['target']} ({top_hit['taxname']})"
            + (f"; {len(identity_hits)}/{top_k} hits contain '{reference_name}'" if identity_hits else "")
        ),
        "top_hit": top_hit,
        "identity_hits": identity_hits,
    }


def write_report(
    construct_name: str,
    hits: list[dict],
    confirmation: dict,
    output_dir: str,
):
    """Write validation report CSV and summary."""
    os.makedirs(output_dir, exist_ok=True)

    # Full hits
    hits_path = os.path.join(output_dir, f"{construct_name}_foldseek_hits.csv")
    if hits:
        with open(hits_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(hits[0].keys()))
            writer.writeheader()
            writer.writerows(hits)

    # Summary
    summary_path = os.path.join(output_dir, f"{construct_name}_validation_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"Construct Validation Report: {construct_name}\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Identity confirmed: {confirmation['confirmed']}\n")
        f.write(f"Reason: {confirmation['reason']}\n\n")

        if confirmation["top_hit"]:
            f.write("Top FoldSeek Hit:\n")
            for k, v in confirmation["top_hit"].items():
                f.write(f"  {k}: {v}\n")

        f.write(f"\nTotal hits above threshold: {len(hits)}\n")
        if confirmation["identity_hits"]:
            f.write(f"Identity-confirmed hits: {len(confirmation['identity_hits'])}\n")
            for h in confirmation["identity_hits"][:5]:
                f.write(f"  TM={h['tm_score']:.3f} | {h['target']} | {h['taxname']}\n")

    print(f"\n{'✓' if confirmation['confirmed'] else '✗'} Construct identity: {confirmation['confirmed']}")
    print(f"  {confirmation['reason']}")
    print(f"  Report written to {summary_path}")

    return summary_path


def main():
    parser = argparse.ArgumentParser(
        description="Validate Frankenstein construct identity via FoldSeek"
    )
    parser.add_argument("--construct", required=True,
                        help="Relaxed stitched construct PDB (output of run_md_relaxation.py)")
    parser.add_argument("--reference_name", required=True,
                        help="Expected protein name to find in FoldSeek hits (e.g. 'Tau', 'MAPT')")
    parser.add_argument("--foldseek_db", default="pdb",
                        help="FoldSeek database to search (default: pdb). "
                             "Can also be a path to a custom DB.")
    parser.add_argument("--output", required=True,
                        help="Output directory for validation reports")
    parser.add_argument("--min_tm", type=float, default=0.5,
                        help="Minimum TM-score for a hit to be reported (default: 0.5)")
    parser.add_argument("--top_k", type=int, default=10,
                        help="Top-K hits to check for identity confirmation (default: 10)")
    parser.add_argument("--tmp_dir", default="tmp/foldseek_validate",
                        help="Temporary directory for FoldSeek")
    args = parser.parse_args()

    construct_name = Path(args.construct).stem

    # Run FoldSeek
    hits_tsv = os.path.join(args.output, f"{construct_name}_raw.tsv")
    run_foldseek_search(
        query_pdb=args.construct,
        db=args.foldseek_db,
        output_tsv=hits_tsv,
        tmp_dir=args.tmp_dir,
    )

    # Parse and filter
    hits = parse_foldseek_hits(hits_tsv, min_tm=args.min_tm)
    print(f"Found {len(hits)} hits with TM-score >= {args.min_tm}")

    # Confirm identity
    confirmation = confirm_identity(hits, args.reference_name, top_k=args.top_k)

    # Write report
    write_report(construct_name, hits, confirmation, args.output)

    # Exit code: 0 = confirmed, 1 = not confirmed
    return 0 if confirmation["confirmed"] else 1


if __name__ == "__main__":
    exit(main())
