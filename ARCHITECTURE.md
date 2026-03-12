# System Architecture

## Overview

The IDP Binder Generation pipeline is organized into four sequential phases, each with its own source module. Data flows one-way through the pipeline: structures → binders → motifs → validation.

```
data/constructs/          (input sequences + metadata)
        │
        ▼
src/structure_prep/       Phase 1: ensemble generation + stitching + MD
        │
        ▼
src/binder_generation/    Phase 2: RFDiffusion / BindCraft / Forge
        │
        ▼
src/motif_analysis/       Phase 3: FoldMason + ESM2 embeddings + FoldSeek
        │
        ▼
src/validation/           Phase 4: recall scoring + known binder benchmarking
        │
        ▼
results/                  (all outputs, gitignored)
```

---

## Phase 1 — Structure Preparation (`src/structure_prep/`)

### Data flow
```
data/constructs/tau/sequences.md
        │
        ▼
[Starling — external]           → data/structures/{construct}_ensemble/*.pdb
        │
        ├── IDP-only constructs: pass directly to Phase 2
        │
        └── IDP + stable domain constructs:
                │
                ├── foldmason_refine.py       → data/structures/{construct}_refined/*.pdb
                │       - FoldMason easy-msa on full ensemble
                │       - Per-conformer consistency score (vs. MSA consensus)
                │       - Combined rank: 0.6×consistency + 0.4×lDDT
                │       - Keep top-N representative conformers
                │
                ├── [AlphaFold2 — external]  → data/structures/{construct}_af2.pdb
                │
                ├── stitch_constructs.py      → data/structures/{construct}_stitched/*.pdb
                │       - Cα superposition at junction (±5 residues)
                │       - SVD rotation+translation
                │       - Atom coordinate transform for IDP region
                │
                └── run_md_relaxation.py      → data/structures/{construct}_relaxed.pdb
                        - AMBER14-SB + TIP3P
                        - Cα restraints on stable domain
                        - Energy minimize → NVT → NPT
```

### Key design decisions
- **Restraint during MD**: Only the IDP region relaxes freely; stable domain Cα atoms are harmonically restrained (k = 1000 kJ/mol/nm²). This prevents the AF2 stable domain from drifting while letting the disordered junction adapt.
- **Multiple conformers**: Starling generates an ensemble; we take 3–5 representative members (clustered by RMSD) to capture IDP heterogeneity without combinatorial explosion downstream.
- **NFT exception**: Tau NFT uses published cryo-EM coordinates (PDB: 6QJH) directly — no Starling needed.

---

## Phase 2 — Binder Generation (`src/binder_generation/`)

Two parallel tracks that produce complementary candidate sets.

### Track A: Structure-based (RFDiffusion / BindCraft)
```
data/structures/{construct}_relaxed.pdb
        │
        ├── RFDiffusion (backbone diffusion)
        │       - Input: target PDB + hotspot residues
        │       - Output: backbone-only binder scaffolds (~100/construct)
        │       │
        │       └── ProteinMPNN (sequence design)
        │               - Input: binder backbone + target
        │               - Output: 8 sequences per scaffold
        │               │
        │               └── ESMFold (folding + confidence)
        │                       - Filter: ipTM > 0.7, pLDDT > 70
        │
        └── BindCraft (end-to-end hallucination)
                - Input: target PDB + chain
                - Output: all-atom binder designs (~50/construct)
```

### Track B: Sequence-based (Forge)
```
data/constructs/{construct}.fasta
        │
        └── Forge (flow matching on Raygun/ESM2 embeddings)
                - Input: target sequence
                - Output: 100 generated binder sequences/construct
                │
                └── ESMFold
                        - Filter: pLDDT > 70
```

### Construct × method matrix
| Construct | RFDiffusion | BindCraft | Forge |
|-----------|-------------|-----------|-------|
| tau_monomer | ✓ (PHF6/PHF6* hotspots) | ✓ | ✓ |
| tau_oligomer | — (structure too undefined) | ✓ | ✓ |
| tau_nft | ✓ (surface grooves on 6QJH) | — | ✓ |
| targetA_construct | ✓ | ✓ | ✓ |
| targetA_full | ✓ | — | ✓ |

### Output
All filtered candidates merged into: `results/binders/all_candidates.fasta` + individual ESMFold PDBs in `results/esmfold/{construct}/`

---

## Phase 3 — Motif Analysis (`src/motif_analysis/`)

Three independent analyses that converge on a ranked candidate list.

```
results/esmfold/{construct}/*.pdb          results/binders/all_candidates.fasta
        │                                          │
        ▼                                          ▼
foldmason_parser.py                        motif_scanner.py
  FoldMason structural MSA                   ESM2 (esm2_t33_650M_UR50D)
  → per-column conservation                  → per-residue embeddings (1280d)
  → conserved motif windows                  → mean pool per sequence
  → results/motifs/foldmason_motifs.csv      → UMAP (30 neighbors, 0.1 min_dist)
                                             → HDBSCAN (min_cluster=20)
        │                                    → results/embeddings/umap_clusters.csv
        │                                          │
        └──────────────────┬────────────────────────┘
                           │
                           ▼
                   foldseek_utils.py
                     FoldSeek search
                     → results/esmfold/ vs data/known_binders/
                     → results/foldseek/annotated_hits.csv
                           │
                           ▼
                   Integrated Motif Score
                   = w1×FoldMason_conservation
                   + w2×HDBSCAN_cluster_quality
                   + w3×FoldSeek_recall
                   + w4×ipTM
                   → results/motifs/integrated_ranking.csv
```

### ESM2 embedding rationale
Single-residue embeddings from ESM2 encode:
- Local sequence context (neighbors in sequence)
- Secondary structure propensity
- Evolutionary co-variation signal

By pooling interface residues and clustering, we ask: **do binders for different Tau states share a common "chemical language"?** Cluster overlap = state-agnostic motif. Cluster separation = state-specific binder.

---

## Phase 4 — Validation (`src/validation/`)

```
results/motifs/integrated_ranking.csv      data/known_binders/metadata.csv
        │                                          │
        └──────────────────┬────────────────────────┘
                           │
                           ▼
                   benchmark_binders.py
                     - Map clusters → known binder DB entries
                     - Flag: positive_control_hit | failed_drug_hit | novel
                     - Compute recall@K
                     → results/validation/benchmark_report.csv
```

### Validation tiers
```
Tier 1: Structured control (EGFR / Hemoglobin)
    → Confirms pipeline produces sensible binders for well-defined targets
    → Pass: FoldSeek recall ≥ 10% of known binders in top-50

Tier 2: Baker 2025 IDP binder controls
    → Confirms embedding space is consistent with published solutions
    → Pass: any cluster overlap between our binders and Baker 2025 binders

Tier 3: Tau-specific validation
    → AT8, PHF1, MC1, HJ8.5 epitope overlap
    → Failed drug region flagging (LMTM, Semorinemab, Gosuranemab)
```

---

## Data Directory Structure

```
data/
├── constructs/
│   ├── tau/
│   │   ├── sequences.md          ← K18 / NFT sequences + UniProt references
│   │   ├── tau_monomer.fasta     ← (to add)
│   │   ├── tau_oligomer.fasta    ← (to add)
│   │   └── tau_nft.fasta         ← (to add)
│   └── target_protein/
│       ├── construct.fasta       ← (to add — from collaborator)
│       └── full.fasta            ← (to add — from collaborator)
├── known_binders/
│   ├── metadata.csv              ← curated DB: antibodies, drugs, Baker controls
│   └── *.pdb                     ← binder structures (to add as available)
├── failed_drugs/
│   └── *.sdf / *.pdb             ← LMTM etc. (to add)
└── structures/                   ← generated/relaxed PDBs (gitignored)
```

---

## Results Directory Structure

```
results/                          ← gitignored
├── starling/{construct}/         ← ensemble PDBs
├── esmfold/{construct}/          ← folded binder PDBs + ipTM scores
├── rfdiffusion/{construct}/      ← RFDiffusion backbone outputs
├── mpnn/{construct}/             ← ProteinMPNN sequences
├── bindcraft/{construct}/        ← BindCraft outputs
├── forge/{construct}/            ← Forge generated sequences
├── binders/all_candidates.fasta  ← merged filtered candidates
├── foldmason/{run}/              ← FoldMason structural MSA outputs
├── embeddings/
│   ├── embeddings_raw.npz
│   ├── umap_clusters.csv
│   └── cluster_summary.csv
├── foldseek/
│   ├── hits.tsv
│   └── annotated_hits.csv
├── motifs/
│   ├── foldmason_motifs.csv
│   └── integrated_ranking.csv
└── validation/
    └── benchmark_report.csv
```

---

## External Tool Dependencies

| Tool | Version | Install | Used in |
|------|---------|---------|---------|
| Starling | latest | `pip install starling-protein` | Phase 1 |
| AlphaFold2 / ColabFold | 2.3+ | colabfold or local | Phase 1 |
| OpenMM | 8.0+ | `conda install -c conda-forge openmm` | Phase 1 |
| RFDiffusion | latest | GitHub RosettaCommons | Phase 2A |
| BindCraft | latest | GitHub martinpacesa | Phase 2A |
| ProteinMPNN | latest | GitHub dauparas | Phase 2A |
| Forge | e80 checkpoint | HuggingFace yk0/forge-e80 | Phase 2B |
| ESMFold | via ESM lib | `pip install fair-esm` | Phase 2 + 3 |
| ESM2 (650M) | esm2_t33_650M | `pip install fair-esm` | Phase 3 |
| FoldMason | latest | `conda install foldmason` | Phase 3 |
| FoldSeek | latest | `conda install foldseek` | Phase 3 |

---

## Module Interface Summary

| Module | Inputs | Outputs |
|--------|--------|---------|
| `foldmason_refine.py` | Starling ensemble dir | ranked/filtered conformer PDBs + ranking CSV |
| `stitch_constructs.py` | stable PDB, IDP ensemble dir, residue range | stitched PDBs |
| `run_md_relaxation.py` | stitched PDB, IDP residue range | relaxed PDB |
| `motif_scanner.py` | FASTA of binder sequences | UMAP CSV, cluster CSV, raw embeddings NPZ |
| `foldmason_parser.py` | FoldMason MSA FASTA | motifs CSV |
| `foldseek_utils.py` | ESMFold PDB dir, known binder DB | annotated hits CSV |
| `benchmark_binders.py` | ranked candidates CSV, embeddings CSV | benchmark report CSV |
