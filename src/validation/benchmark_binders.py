"""
Benchmark generated binders against known binders and failed clinical compounds.

Loads the known binder validation set, computes recall at different thresholds,
and scores each generated binder against the validation set.

Usage:
    python benchmark_binders.py \
        --generated results/motifs/integrated_ranking.csv \
        --known_binders data/known_binders/metadata.csv \
        --embeddings results/embeddings/umap_clusters.csv \
        --output results/validation/benchmark_report.csv
"""

import argparse
import csv
import os
from collections import defaultdict


KNOWN_BINDER_DB = {
    # Tau antibodies (epitope residues, 0-indexed in repeat domain)
    "AT8_tau_pS202_pT205": {
        "target": "tau",
        "type": "antibody",
        "epitope_residues": [202, 205],
        "state": "monomer",
        "status": "clinical_diagnostic",
        "notes": "Phospho-specific; marks disease Tau",
    },
    "PHF1_tau_pS396_pS404": {
        "target": "tau",
        "type": "antibody",
        "epitope_residues": [396, 404],
        "state": "monomer",
        "status": "research",
        "notes": "C-terminal phospho epitope",
    },
    "MC1_tau_conformational": {
        "target": "tau",
        "type": "antibody",
        "epitope_residues": [7, 8, 9, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322],
        "state": "oligomer",
        "status": "research",
        "notes": "Conformational epitope; disease-specific",
    },
    "HJ8.5_tau_Nterminal": {
        "target": "tau",
        "type": "antibody",
        "epitope_residues": [25, 26, 27, 28, 29, 30],
        "state": "monomer",
        "status": "research",
        "pdb": "6NWP",
        "notes": "Crystal structure available",
    },
    # Failed clinical trial drugs
    "LMTM_TRx0237": {
        "target": "tau",
        "type": "small_molecule",
        "epitope_residues": [],  # binds repeat domain, exact contacts unclear
        "state": "monomer",
        "status": "failed_phase3",
        "notes": "Tau aggregation inhibitor; Phase 3 failed 2016",
    },
    "Semorinemab": {
        "target": "tau",
        "type": "antibody",
        "epitope_residues": list(range(1, 50)),  # N-terminal
        "state": "monomer",
        "status": "failed_phase2",
        "notes": "Anti-Tau N-terminal; Phase 2 failed 2021",
    },
    "Gosuranemab": {
        "target": "tau",
        "type": "antibody",
        "epitope_residues": list(range(6, 24)),
        "state": "monomer",
        "status": "failed_phase2",
        "notes": "N-terminal aa 6-23; Phase 2 failed 2021",
    },
    # Baker lab IDP binder controls (fill in after reading papers)
    "Baker2025_Nature_IDP_binder_1": {
        "target": "IDP_general",
        "type": "de_novo_protein",
        "epitope_residues": [],
        "state": "monomer",
        "status": "published_positive_control",
        "notes": "Baker lab 2025 Nature; ADD PDB when available",
    },
    "Baker2025_Science_IDP_binder_1": {
        "target": "IDP_general",
        "type": "de_novo_protein",
        "epitope_residues": [],
        "state": "monomer",
        "status": "published_positive_control",
        "notes": "Baker lab 2025 Science; ADD PDB when available",
    },
    # Structured protein controls
    "Cetuximab_EGFR": {
        "target": "EGFR",
        "type": "antibody",
        "epitope_residues": list(range(350, 420)),
        "state": "structured",
        "status": "approved_drug",
        "pdb": "1YY9",
        "notes": "EGFR domain III binder; positive control",
    },
}


def load_generated_binders(csv_path: str) -> list[dict]:
    with open(csv_path) as f:
        return list(csv.DictReader(f))


def load_embeddings(csv_path: str) -> dict[str, dict]:
    """Load UMAP cluster assignments keyed by header."""
    result = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            result[row["header"]] = row
    return result


def score_binders_against_db(
    generated: list[dict],
    known_db: dict,
    embeddings: dict,
) -> list[dict]:
    """
    For each generated binder, compute:
    - which known binders it clusters with (same HDBSCAN cluster as known binder)
    - summary flags: positive_control_hit, failed_drug_hit, novel
    """
    # Map: cluster -> list of known binder names that fall in cluster
    # (This requires known binders to be embedded + clustered too;
    # here we use a placeholder — populate after running motif_scanner on known binders)
    cluster_to_known = defaultdict(list)
    for name, meta in known_db.items():
        cluster = meta.get("cluster", None)
        if cluster is not None:
            cluster_to_known[cluster].append(name)

    scored = []
    for binder in generated:
        header = binder.get("header", binder.get("name", "unknown"))
        emb_info = embeddings.get(header, {})
        cluster = emb_info.get("cluster", "-1")

        known_hits = cluster_to_known.get(cluster, [])
        positive_hits = [k for k in known_hits
                         if known_db[k]["status"] in ("published_positive_control", "approved_drug", "research")]
        failed_drug_hits = [k for k in known_hits
                            if known_db[k]["status"].startswith("failed")]

        scored.append({
            **binder,
            "cluster": cluster,
            "known_binder_hits": "|".join(known_hits) if known_hits else "",
            "positive_control_hit": len(positive_hits) > 0,
            "failed_drug_hit": len(failed_drug_hits) > 0,
            "novel": len(known_hits) == 0,
            "n_known_hits": len(known_hits),
        })

    return scored


def compute_recall(
    scored: list[dict],
    known_db: dict,
    top_k: int = 50,
) -> dict:
    """
    Compute recall of known binder cluster coverage in top-K candidates.
    """
    top_k_binders = scored[:top_k]
    hit_clusters = set(b["cluster"] for b in top_k_binders if b["cluster"] != "-1")

    # All known binder clusters
    known_clusters = set(str(v.get("cluster")) for v in known_db.values()
                         if v.get("cluster") is not None)

    if not known_clusters:
        return {"recall_at_k": None, "note": "No cluster assignments for known binders yet"}

    recall = len(hit_clusters & known_clusters) / len(known_clusters) if known_clusters else 0
    return {
        f"recall_at_{top_k}": recall,
        "known_clusters_covered": len(hit_clusters & known_clusters),
        "total_known_clusters": len(known_clusters),
    }


def main():
    parser = argparse.ArgumentParser(description="Benchmark generated binders vs. known binders")
    parser.add_argument("--generated", required=True,
                        help="CSV of generated binders (from motif analysis ranking)")
    parser.add_argument("--embeddings", required=True,
                        help="CSV of UMAP cluster assignments (from motif_scanner.py)")
    parser.add_argument("--output", required=True, help="Output benchmark report CSV")
    parser.add_argument("--top_k", type=int, default=50,
                        help="Top-K candidates for recall computation")
    args = parser.parse_args()

    generated = load_generated_binders(args.generated)
    embeddings = load_embeddings(args.embeddings)

    print(f"Loaded {len(generated)} generated binders")
    print(f"Loaded {len(embeddings)} embedding assignments")
    print(f"Validation DB: {len(KNOWN_BINDER_DB)} known binders/drugs")

    scored = score_binders_against_db(generated, KNOWN_BINDER_DB, embeddings)

    recall = compute_recall(scored, KNOWN_BINDER_DB, top_k=args.top_k)
    print(f"\nRecall metrics: {recall}")

    # Summary
    n_pos = sum(1 for b in scored if b["positive_control_hit"])
    n_failed = sum(1 for b in scored if b["failed_drug_hit"])
    n_novel = sum(1 for b in scored if b["novel"])
    print(f"\nTop {len(scored)} binders:")
    print(f"  Positive control hits: {n_pos}")
    print(f"  Failed drug region hits: {n_failed}")
    print(f"  Novel (no known cluster): {n_novel}")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as f:
        if scored:
            writer = csv.DictWriter(f, fieldnames=list(scored[0].keys()))
            writer.writeheader()
            writer.writerows(scored)

    print(f"\nBenchmark report written to {args.output}")


if __name__ == "__main__":
    main()
