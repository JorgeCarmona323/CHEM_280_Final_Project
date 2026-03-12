"""
Single-residue ESM2 embedding extraction + UMAP + HDBSCAN clustering
for motif discovery across generated binder sequences.

Usage:
    python motif_scanner.py \
        --sequences results/mpnn/all_constructs/sequences.fasta \
        --output_dir results/embeddings/ \
        --esm_model esm2_t33_650M_UR50D
"""

import argparse
import os
import csv
from pathlib import Path
import numpy as np


def read_fasta(fasta_path: str) -> list[tuple[str, str]]:
    """Read FASTA file, return list of (header, sequence) tuples."""
    records = []
    header, seq_lines = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records


def extract_esm2_embeddings(
    sequences: list[tuple[str, str]],
    model_name: str = "esm2_t33_650M_UR50D",
    batch_size: int = 8,
    device: str = "auto",
) -> dict[str, np.ndarray]:
    """
    Extract per-residue ESM2 embeddings for each sequence.

    Returns dict: {header -> array of shape (seq_len, embedding_dim)}
    """
    try:
        import torch
        import esm
    except ImportError:
        raise ImportError("Install ESM: pip install fair-esm")

    if device == "auto":
        import torch
        device = "cuda" if torch.cuda.is_available() else "cpu"

    import torch
    print(f"Loading ESM2 model: {model_name} on {device}")
    model, alphabet = esm.pretrained.__dict__[model_name]()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    embeddings = {}
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        batch_labels, batch_strs, batch_tokens = batch_converter(batch)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=False)

        token_representations = results["representations"][33]

        for j, (header, seq) in enumerate(batch):
            # Trim off BOS/EOS tokens
            emb = token_representations[j, 1:len(seq) + 1].cpu().numpy()
            embeddings[header] = emb

        print(f"  Processed {min(i + batch_size, len(sequences))}/{len(sequences)}")

    return embeddings


def pool_embeddings(
    embeddings: dict[str, np.ndarray],
    method: str = "mean",
) -> tuple[np.ndarray, list[str]]:
    """
    Pool per-residue embeddings to per-sequence vectors.

    method: 'mean' | 'max' | 'cls' (first token)
    Returns (matrix, headers)
    """
    headers = list(embeddings.keys())
    vectors = []
    for h in headers:
        emb = embeddings[h]
        if method == "mean":
            vectors.append(emb.mean(axis=0))
        elif method == "max":
            vectors.append(emb.max(axis=0))
        elif method == "cls":
            vectors.append(emb[0])
        else:
            raise ValueError(f"Unknown pooling method: {method}")
    return np.array(vectors), headers


def run_umap_hdbscan(
    vectors: np.ndarray,
    headers: list[str],
    output_dir: str,
    n_neighbors: int = 30,
    min_dist: float = 0.1,
    min_cluster_size: int = 20,
):
    """Run UMAP dimensionality reduction + HDBSCAN clustering, save results."""
    try:
        import umap
    except ImportError:
        raise ImportError("Install UMAP: pip install umap-learn")
    try:
        import hdbscan
    except ImportError:
        raise ImportError("Install HDBSCAN: pip install hdbscan")

    print(f"Running UMAP on {len(vectors)} sequences...")
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, random_state=42)
    embedding_2d = reducer.fit_transform(vectors)

    print("Running HDBSCAN clustering...")
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        prediction_data=True,
    )
    labels = clusterer.fit_predict(embedding_2d)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = (labels == -1).sum()
    print(f"  Found {n_clusters} clusters, {n_noise} noise points")

    # Save UMAP + cluster assignments
    os.makedirs(output_dir, exist_ok=True)
    umap_path = os.path.join(output_dir, "umap_clusters.csv")
    with open(umap_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["header", "umap_x", "umap_y", "cluster"])
        for h, (x, y), label in zip(headers, embedding_2d, labels):
            writer.writerow([h, f"{x:.4f}", f"{y:.4f}", int(label)])
    print(f"UMAP clusters saved to {umap_path}")

    # Cluster summary
    summary_path = os.path.join(output_dir, "cluster_summary.csv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cluster", "size", "fraction"])
        total = len(labels)
        for c in sorted(set(labels)):
            size = (labels == c).sum()
            label_str = "noise" if c == -1 else str(c)
            writer.writerow([label_str, size, f"{size / total:.3f}"])
    print(f"Cluster summary saved to {summary_path}")

    return embedding_2d, labels


def main():
    parser = argparse.ArgumentParser(description="ESM2 embedding + UMAP + HDBSCAN motif scanner")
    parser.add_argument("--sequences", required=True, help="FASTA file of binder sequences")
    parser.add_argument("--output_dir", required=True, help="Output directory for embeddings/clusters")
    parser.add_argument("--esm_model", default="esm2_t33_650M_UR50D",
                        help="ESM2 model name (default: esm2_t33_650M_UR50D)")
    parser.add_argument("--pooling", default="mean", choices=["mean", "max", "cls"],
                        help="Embedding pooling method")
    parser.add_argument("--umap_n_neighbors", type=int, default=30)
    parser.add_argument("--umap_min_dist", type=float, default=0.1)
    parser.add_argument("--hdbscan_min_cluster_size", type=int, default=20)
    parser.add_argument("--batch_size", type=int, default=8)
    args = parser.parse_args()

    sequences = read_fasta(args.sequences)
    print(f"Loaded {len(sequences)} sequences from {args.sequences}")

    embeddings = extract_esm2_embeddings(
        sequences,
        model_name=args.esm_model,
        batch_size=args.batch_size,
    )

    # Save raw embeddings as numpy archive
    os.makedirs(args.output_dir, exist_ok=True)
    raw_path = os.path.join(args.output_dir, "embeddings_raw.npz")
    np.savez(raw_path, **{h.replace("/", "_"): v for h, v in embeddings.items()})
    print(f"Raw embeddings saved to {raw_path}")

    vectors, headers = pool_embeddings(embeddings, method=args.pooling)

    run_umap_hdbscan(
        vectors=vectors,
        headers=headers,
        output_dir=args.output_dir,
        n_neighbors=args.umap_n_neighbors,
        min_dist=args.umap_min_dist,
        min_cluster_size=args.hdbscan_min_cluster_size,
    )


if __name__ == "__main__":
    main()
