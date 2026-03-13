# Notebook 03: Binder Comparison & Motif Analysis — PLAN

> **TODO**: Convert this into `03_binder_comparison.ipynb`
> Reminder: build this after Forge and BindCraft outputs are in hand and filtered through 03a

---

## Pipeline Position

```
[Notebook 01] Starling ensemble → filtered frames
[Notebook 02] FoldMason → top 3–5 conformers (_binder_ready/*.pdb)
    ↓
RFDiffusion / BindCraft  →  binder structures + ipTM scores
Forge                    →  binder sequences → ESMFold → structures
    ↓
[Notebook 03a] Binder contact filter (run SEPARATELY per method)
    → forge_binder_filter_passing.csv   + binders_passing/*.pdb
    → bindcraft_binder_filter_passing.csv + binders_passing/*.pdb
    ↓
[THIS NOTEBOOK] Binder comparison, motif analysis, validation
    ↓
integrated_ranking.csv  →  experimental shortlist
```

---

## Inputs

| File | Source |
|------|--------|
| `forge_binder_filter_passing.csv` | Notebook 03a (Forge run) |
| `bindcraft_binder_filter_passing.csv` | Notebook 03a (BindCraft run) |
| `results/esmfold/forge/*.pdb` | ESMFolded Forge binder structures |
| `results/esmfold/bindcraft/*.pdb` | BindCraft output structures |
| `data/known_binders/sequences.fasta` | Control binders pulled from PDB/literature (see below) |
| `data/known_binders/metadata.csv` | Kd, target epitope, source per control binder |

---

## Embedding Strategy — SaProt (recommended)

### Why SaProt over ESM2 alone

Running ESM2 (sequence only) and 3Di (structure only) separately and fusing them
is valid but introduces normalization choices and potential modality imbalance.

**SaProt** (`westlake-repl/SaProt_650M_AF2`, HuggingFace) is trained on
**(sequence token, 3Di token)** pairs simultaneously — it natively encodes both
modalities in a single 1280-dim embedding with no fusion step needed.

Input format: interleaved sequence + 3Di tokens per residue:
```
K v Q a I i I n K k ...
↑   ↑   ↑   ↑   ↑       (uppercase = amino acid, lowercase = 3Di state)
```
Output: 1280-dim embedding per residue → mean-pool → one vector per binder

### Fallback — Early fusion (ESM2 + 3Di concatenation)

If SaProt is unavailable or too slow on Colab free tier:
1. ESM2 (`esm2_t33_650M_UR50D`) → 1280-dim per binder (mean-pool)
2. 3Di from FoldMason JSON → one-hot encode 20 states × N residues → mean-pool → 20-dim
3. L2-normalize each independently, then concatenate → 1300-dim
4. PCA to 50-dim → UMAP

### Why merging is valid and powerful

| Embedding | Captures |
|-----------|----------|
| ESM2 / SaProt sequence | Chemical identity, evolutionary context, local aa environment |
| 3Di / SaProt structural | Backbone geometry, interface shape, fold topology |

A cluster present in **both** modalities = conserved in sequence AND structure →
strongest signal, most likely a real binding mode. Cluster in sequence only =
sequence convergence without structural convergence (weaker). SaProt captures both
in one pass.

---

## Step 1 — ESMFold all binders (GPU, Colab T4)

Fold every binder sequence to get structure + per-residue pLDDT.
- Forge: sequences only → ESMFold each
- BindCraft: structures already available → ESMFold for consistency or use as-is
- Known controls: fold any without an existing PDB structure

```python
import esm
model = esm.pretrained.esmfold_v1()
# output: PDB string with pLDDT in B-factor column
```

---

## Step 2 — SaProt Embeddings → UMAP → HDBSCAN

```python
from transformers import EsmTokenizer, EsmModel

tokenizer = EsmTokenizer.from_pretrained("westlake-repl/SaProt_650M_AF2")
model     = EsmModel.from_pretrained("westlake-repl/SaProt_650M_AF2")

# Input: interleaved (aa, 3Di) token string per binder
# Output: mean-pooled 1280-dim vector per binder
```

UMAP parameters: `n_neighbors=30, min_dist=0.1, metric="cosine"`
HDBSCAN: `min_cluster_size=10` (smaller than before — filtered set is smaller)

**Color the UMAP by**:
- Method (Forge / BindCraft) — overlap = convergent solutions = high confidence
- ipTM score — are high-confidence binders clustered?
- Target construct (monomer / oligomer / NFT) — state-specific signatures?
- Source (generated / known control) — where do AT8, HJ8.5, Baker 2025 land?
- Epitope contacted (PHF6 / PHF6* / both / neither)

**Key interpretations**:

| Observation | Meaning |
|-------------|---------|
| Forge + BindCraft overlap in same cluster | Both methods independently found same binding mode → top candidates |
| Cluster contains known binder (AT8, HJ8.5) | Validated binding mode, known epitope |
| Cluster contains Baker 2025 miniprotein | Structurally precedented IDP binder architecture |
| Cluster near failed drug (LMTM) | Caution — same binding mode as clinical failure |
| Cluster with NO known binder | Novel binding mode — most scientifically interesting |

---

## Step 3 — Contact Map → Epitope Heatmap

For each passing binder complex, compute which target residues are contacted (8Å Cα–Cα).
Aggregate contact frequency per target residue across all binders per method.

Plot: contact frequency heatmap across target sequence, Forge vs BindCraft side by side.
- Peaks = dominant binding epitopes
- Annotate PHF6* (2–7) and PHF6 (33–38) on x-axis
- Shared peaks between methods = most robust epitopes

---

## Step 4 — FoldMason Structural MSA on All Binders

Run FoldMason on merged binder set (Forge passing + BindCraft passing + known controls).

```bash
foldmason easy-msa \
    results/esmfold/forge/*.pdb \
    results/esmfold/bindcraft/*.pdb \
    data/known_binders/*.pdb \
    results/foldmason/all_binders/msa \
    tmp/foldmason_all_binders/ \
    --report-mode 2
```

Parse `foldmason.json`:
- lDDT per binder → rank by structural quality
- 3Di conservation plot → positions with non-D dominant state = structurally locked interface
- Do Forge + BindCraft binders share a conserved structural motif? → binding signature

---

## Step 5 — Chemical Validation

Per binder using Biopython `ProteinAnalysis`:

| Property | Flag if |
|----------|---------|
| GRAVY score | > 0.5 (too hydrophobic → aggregation/solubility risk) |
| Net charge pH 7 | < −10 or > +10 |
| Isoelectric point | < 4 or > 10 |
| Secondary structure (helix+sheet fraction) | < 20% (likely disordered binder) |
| Instability index | > 40 (unstable in solution) |

Overlay chemical properties on UMAP — flag clusters enriched for problematic binders.

---

## Step 6 — Known Control Binder Sources

### Tau binders (pull from PDB / literature)

| Control | PDB / Source | Type | Notes |
|---------|-------------|------|-------|
| AT8 Fab | 6NWP | Phospho-Tau antibody (S202/T205) | Fab structure + sequence |
| HJ8.5 nanobody | 6NWQ | MTBR binder | Nanobody, good structural reference |
| MC1 | Literature (Jicha et al.) | Conformational NFT antibody | Sequence from patent |
| PHF1 | 7PQ4 | Phospho-Tau (S396/S404) | |
| Baker 2025 miniproteins | Nature SI Table S2 / Science SI | IDP miniprotein binders | amylin 3.8 nM, G3BP1 11 nM |
| Semorinemab | WHO INN / patent WO2016112078 | Clinical antibody (failed Ph2) | N-terminal Tau epitope |
| Gosuranemab | WHO INN / patent WO2017059186 | Clinical antibody (failed Ph2) | N-terminal Tau epitope |

### FMRP binders (limited — use as negative control space)

| Control | Source | Notes |
|---------|--------|-------|
| CYFIP1 binding interface | PDB 4UMJ | Protein-protein interface residues |
| NUFIP1 | Literature | No structure — sequence only |

### How to pull sequences from PDB

```bash
# Download FASTA for a PDB entry
curl https://www.rcsb.org/fasta/entry/6NWQ > data/known_binders/hj8.5_6nwq.fasta

# Or use Biopython
from Bio import Entrez, SeqIO
# fetch by PDB ID via RCSB REST API
```

All known binder sequences go in `data/known_binders/sequences.fasta` with metadata
in `data/known_binders/metadata.csv` (name, source, type, Kd, epitope, PDB_id).

---

## Step 7 — Integrated Binder Score

```
Integrated Score = w1 × ipTM
                 + w2 × ESMFold_pLDDT (interface residues only)
                 + w3 × HDBSCAN cluster confidence (DBCV score)
                 + w4 × epitope_contact_score  (contacts PHF6 / PHF6*?)
                 + w5 × chemical_score         (GRAVY, charge, instability)
                 + w6 × known_binder_proximity (near AT8/HJ8.5/Baker2025 in UMAP)
                 − w7 × failed_drug_proximity  (penalty if near LMTM/Semorinemab)
```

Rank all candidates → top binders enter experimental shortlist.

---

## Output Files

| File | Description |
|------|-------------|
| `results/embeddings/umap_plot.png` | Master UMAP — all binders + controls, colored by method/ipTM/cluster |
| `results/embeddings/umap_clusters.csv` | UMAP coords, cluster ID, method, ipTM, source per binder |
| `results/contacts/epitope_heatmap.png` | Contact frequency across target sequence, Forge vs BindCraft |
| `results/contacts/epitope_map.csv` | Per-residue contact frequency per method |
| `results/foldmason/all_binders/msa_aa.fa` | Structural MSA of all binders |
| `results/chemical/properties.csv` | Per-binder chemical properties |
| `results/integrated_ranking.csv` | Final ranked candidate list → experimental shortlist |
