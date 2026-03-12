# Phase 3: Motif Analysis

## Goal
Discover whether generated binders share recurring structural or sequence motifs across different target constructs (monomer vs. oligomer vs. NFT). If they do, those motifs are candidate "binding signatures" worth validating.

---

## 3.1 FoldMason — Structural Multiple Sequence Alignment

FoldMason builds a structural MSA from a set of protein structures, identifying structurally conserved regions even when sequence identity is low.

### Input
All ESMFold-predicted binder structures (from Phase 2), organized by target construct.

```bash
# Install FoldMason
conda install -c conda-forge -c bioconda foldmason

# Build structural MSA for binders against tau_monomer
foldmason easy-msa \
    results/esmfold/tau_monomer/*.pdb \
    results/foldmason/tau_monomer/msa \
    tmp/ \
    --report-mode 1

# Compare across all constructs (look for cross-state motifs)
foldmason easy-msa \
    results/esmfold/tau_monomer/*.pdb \
    results/esmfold/tau_oligomer/*.pdb \
    results/esmfold/tau_nft/*.pdb \
    results/foldmason/all_tau/msa \
    tmp/
```

### What we look for
- Conserved structural elements (α-helices, β-strands, loops) across binders from different states
- Enriched secondary structure patterns at the binding interface
- If the same ~10-20 aa structural motif appears in binders for monomer AND NFT → state-agnostic binding motif

```bash
# Parse FoldMason output and extract conserved motifs
python src/motif_analysis/foldmason_parser.py \
    --msa results/foldmason/all_tau/msa.fasta \
    --pdb_dir results/esmfold/ \
    --output results/motifs/foldmason_motifs.csv \
    --conservation_threshold 0.7
```

---

## 3.2 ESM2 Single-Residue Embedding Cluster Analysis

Each amino acid in a protein sequence gets a 1280-dim embedding vector from ESM2. We use these as a residue-level "fingerprint" for binding interface residues.

### Why single-residue embeddings?
- Capture chemical environment (neighbors, secondary structure, evolutionary context)
- State-sensitive: the same residue in a disordered vs. ordered context gets a different embedding
- Allow comparison across diverse sequences without alignment

### Pipeline

```bash
python src/motif_analysis/motif_scanner.py \
    --sequences results/mpnn/all_constructs/sequences.fasta \
    --output_dir results/embeddings/ \
    --esm_model esm2_t33_650M_UR50D \
    --umap_n_neighbors 30 \
    --umap_min_dist 0.1 \
    --hdbscan_min_cluster_size 20
```

**Steps inside the script**:
1. **Load ESM2** (`esm2_t33_650M_UR50D`, 650M params)
2. **Extract per-residue embeddings** for all binder sequences
3. **Pool interface residues** (top-k residues by AlphaFold2/ESMFold pLDDT × contact score)
4. **UMAP** reduction to 2D (n_neighbors=30, min_dist=0.1)
5. **HDBSCAN** clustering (min_cluster_size=20)
6. **Color by**: construct targeted, design method (RFDiffusion/BindCraft/Forge), ipTM score

### What we look for
- Do binders against different Tau states cluster separately? (state-specific signatures)
- Do binders from Track A (structure-based) and Track B (sequence-based) overlap? (orthogonal validation)
- Are there clusters enriched for high ipTM? (quality-correlated motifs)
- Do known binders (Phase 4 validation set) fall within any cluster? (recall)

---

## 3.3 FoldSeek — Structural Similarity Search

FoldSeek searches a database of known structures for structural neighbors of our generated binders.

```bash
# Create database from validation set + PDB
foldseek createdb data/known_binders/*.pdb DB_known_binders/
foldseek createdb results/esmfold/all_binders/*.pdb DB_generated/

# Search generated binders against known binders
foldseek easy-search \
    results/esmfold/all_binders/*.pdb \
    DB_known_binders/ \
    results/foldseek/hits.tsv \
    tmp/ \
    --format-output "query,target,evalue,alntmscore,rmsd,qaln,taln"

# Also search against full PDB for context
foldseek easy-search \
    results/esmfold/all_binders/*.pdb \
    pdb \
    results/foldseek/pdb_hits.tsv \
    tmp/
```

```bash
# Parse hits and annotate binders with structural neighbors
python src/motif_analysis/foldseek_utils.py \
    --hits results/foldseek/hits.tsv \
    --known_binder_metadata data/known_binders/metadata.csv \
    --output results/foldseek/annotated_hits.csv
```

### What we look for
- Which generated binders have structural precedent in known IDP binders?
- Do FoldSeek hits corroborate FoldMason motifs? (converging evidence)
- Any generated binders hitting failed drugs structurally? → these may occupy the same binding mode and could be deprioritized, OR they are more likely to actually bind (two-sided interpretation)

---

## 3.4 Integrative Motif Score

Combine all three analyses into a single score per binder:

```
Motif Score = w1 × FoldMason_conservation
            + w2 × HDBSCAN_cluster_quality
            + w3 × FoldSeek_recall
            + w4 × ipTM
```

Rank all candidates. Top-scoring binders enter Phase 4 validation.

---

## Output

| File | Description |
|------|-------------|
| `results/motifs/foldmason_motifs.csv` | Conserved structural motifs + position |
| `results/embeddings/umap_clusters.csv` | Cluster assignments + UMAP coords |
| `results/embeddings/cluster_summary.csv` | Per-cluster stats (ipTM, method, target) |
| `results/foldseek/annotated_hits.csv` | Structural neighbors in known binder DB |
| `results/motifs/integrated_ranking.csv` | Final ranked candidate list |
