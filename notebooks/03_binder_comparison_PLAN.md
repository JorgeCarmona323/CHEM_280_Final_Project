# Notebook 03: Binder Comparison & Motif Analysis — PLAN

> **TODO**: Convert this into `03_binder_comparison.ipynb`

---

## Pipeline Position

```
[Notebook 01] Starling ensemble → filtered frames
[Notebook 02] FoldMason → top 3–5 conformers (_binder_ready/*.pdb)
    ↓
RFDiffusion / BindCraft  →  binder structures + ipTM scores
Forge                    →  binder sequences (need ESMFold)
    ↓
[THIS NOTEBOOK] Binder comparison, motif analysis, validation
```

---

## Inputs

| File | Source |
|------|--------|
| `results/esmfold/forge/*.pdb` | ESMFold predictions for Forge binder sequences |
| `results/esmfold/bindcraft/*.pdb` | ESMFold predictions for BindCraft outputs |
| `results/bindcraft/scores.csv` | ipTM, pTM, pLDDT per BindCraft binder |
| `data/known_binders/metadata.csv` | AT8, PHF1, MC1, HJ8.5, Baker 2025 control sequences |

---

## Step 1 — ESMFold all binders (GPU required)

Fold every binder sequence to get a structure + per-residue pLDDT.

- Forge outputs sequences only → ESMFold each one
- BindCraft already outputs structures → ESMFold again for consistency, or use as-is
- Use `esm.pretrained.esmfold_v1()` via the `esm` library
- Output: one PDB per binder in `results/esmfold/{method}/`
- Also fold known control binders (AT8 VHH, PHF1 Fab, Baker 2025 miniproteins)

---

## Step 2 — ESM2 Embeddings → UMAP → HDBSCAN

Per-binder sequence fingerprint using ESM2 (`esm2_t33_650M_UR50D`, 1280-dim).

**Steps**:
1. Load ESM2, extract per-residue embeddings for all binders
2. Mean-pool to get one 1280-dim vector per binder
3. UMAP (n_neighbors=30, min_dist=0.1) → 2D projection
4. HDBSCAN (min_cluster_size=15) → cluster assignments

**Color the UMAP by**:
- Method (Forge / BindCraft) → overlap = convergent solutions
- ipTM score → are high-confidence binders clustered?
- Target construct (monomer / oligomer / NFT) → state-specific signatures?
- Known control (yes/no) → where do AT8, PHF1, Baker 2025 land?

**Key question**: Do Forge and BindCraft binders overlap in embedding space?
- Overlap = both methods independently found the same binding mode → high confidence
- Separation = orthogonal solutions → explore both, different mechanisms

---

## Step 3 — Contact Map → Epitope Mapping

For each binder-target complex, identify which target residues are contacted.

**Method**:
- For each binder PDB + target PDB: compute all pairwise Cα distances
- Contact = any binder residue within 8Å of a target residue
- Map contacts back onto target sequence → contact frequency per target residue

**Output**: contact frequency heatmap across target sequence, per method
- Peaks = dominant binding epitopes
- Compare peaks between Forge and BindCraft → shared epitopes = most robust targets
- Annotate PHF6* (residues 2–7) and PHF6 (33–38) on the plot

---

## Step 4 — FoldMason Structural MSA on All Binders

Run FoldMason on the merged binder set (Forge + BindCraft + controls).

```bash
foldmason easy-msa \
    results/esmfold/forge/*.pdb \
    results/esmfold/bindcraft/*.pdb \
    data/known_binders/*.pdb \
    results/foldmason/all_binders/msa \
    tmp/foldmason_all_binders/ \
    --report-mode 2
```

**What to look for**:
- Conserved structural motif shared across Forge + BindCraft binders → binding signature
- Do control binders share this motif? → validates it's a real binding mode
- 3Di conservation plot: positions with non-D dominant state = structurally locked interface residues

---

## Step 5 — Chemical Validation

Per binder, compute:

| Property | Method | Flag if |
|----------|--------|---------|
| GRAVY score (hydrophobicity) | Biopython ProteinAnalysis | > 0 (too hydrophobic → aggregation risk) |
| Net charge at pH 7 | Biopython | Extreme (< −10 or > +10) |
| Isoelectric point | Biopython | < 4 or > 10 |
| Secondary structure fraction | ESMFold pLDDT + DSSP | < 20% structured → likely disordered binder |
| Self-aggregation risk | AGGRESCAN3D score or Rosetta fa_atr | High score → deprioritize |

Plot chemical properties alongside the UMAP — clusters enriched for good chemical properties vs. problematic ones.

---

## Step 6 — Known Control Binder Injection

Embed these in the same ESM2 UMAP space:

| Control | Type | Expected behavior |
|---------|------|-------------------|
| AT8 (S202/T205) | Phospho-Tau antibody | Should cluster with monomer-targeting binders |
| PHF1 (S396/S404) | Phospho-Tau antibody | Similar |
| MC1 (conformational) | NFT-selective antibody | Should cluster with NFT-targeting binders |
| HJ8.5 | MTBR binder | Should cluster near PHF6 epitope binders |
| Baker 2025 amylin (3.8 nM) | IDP miniprotein binder | Structural reference for what a good IDP binder looks like |
| Baker 2025 G3BP1 (11 nM) | IDP miniprotein binder | Same |
| LMTM (failed drug) | Small molecule (encode as pseudosequence or skip) | Flag if generated binders are nearby |

**Key interpretation**:
- Generated binder near AT8/PHF1/HJ8.5 → validated binding mode, known epitope
- Generated binder near Baker 2025 miniproteins → structurally precedented IDP binder
- Generated binder near LMTM → caution, same binding mode as failed drug
- Generated binder in empty space → novel binding mode, most interesting

---

## Step 7 — Integrated Binder Score

Combine all signals into one score per binder:

```
Integrated Score = w1 × ipTM
                 + w2 × ESMFold_pLDDT
                 + w3 × HDBSCAN_cluster_quality   (silhouette score)
                 + w4 × epitope_contact_score      (contacts PHF6/PHF6*?)
                 + w5 × chemical_score             (GRAVY, charge, aggregation)
                 − w6 × failed_drug_similarity     (penalty if near LMTM etc.)
```

Rank all candidates. Top binders enter experimental validation shortlist.

---

## Output Files

| File | Description |
|------|-------------|
| `results/embeddings/umap_clusters.csv` | UMAP coords + cluster ID + method + ipTM per binder |
| `results/embeddings/umap_plot.png` | Master UMAP figure |
| `results/contacts/epitope_map.csv` | Contact frequency per target residue per method |
| `results/contacts/epitope_map.png` | Contact heatmap |
| `results/foldmason/all_binders/msa_aa.fa` | Structural MSA of all binders |
| `results/chemical/properties.csv` | Per-binder chemical properties |
| `results/integrated_ranking.csv` | Final ranked candidate list |

---

## Notes

- ESMFold step needs GPU (Colab T4 is fine)
- UMAP + HDBSCAN + chemical analysis runs on CPU — can do locally
- FoldMason runs locally (binary install) or on Colab
- Known control binder sequences: source from PDB/UniProt, verify against literature before embedding
