# CHEM 280 Final Project — Pipeline Overview

A plain-language guide to the four notebooks, why each step exists,
and what gets passed between sessions.

---

## High-Level Flow

```
NB01: TauK18 Conformer Ensemble
        ↓  top 5 conformer PDBs
NB02b: BindCraft Binder Design  (×5 sessions, one per conformer)
NB02c: Forge Binder Design      (×1 session, all lengths)
        ↓  passing CSVs + CIF structures
NB03:  Binder Comparison Analysis
        ↓  final ranked binder list
```

All four notebooks run on **Google Colab with a T4 GPU**.
The only files you need to save between sessions are the **CSVs** and the **conformer PDBs**.

---

## NB01 — TauK18 Conformer Ensemble Pipeline

**Goal**: Generate a realistic ensemble of Tau K18 backbone conformations and
select the 5 most structurally informative ones to use as binder design targets.

### Why this matters
Tau K18 is an intrinsically disordered protein (IDP) — it has no single fixed
structure. Different conformations expose different epitopes. Designing binders
against a single static structure would miss most of the accessible binding surface.
By sampling an ensemble and selecting diverse, well-structured conformers, we
give the downstream binder design tools a realistic and varied set of targets.

### Steps
| Step | Tool | What it does |
|------|------|-------------|
| 1 | **Starling** | Generates 400 IDP conformer structures from the K18 sequence using a generative model trained on disordered proteins |
| 2 | **AFRC** (Analytical Flory Random Coil) | Computes sequence-specific polymer reference values (Rg, end-to-end distance). Filters out conformers with unphysical chain dimensions |
| 3 | **FoldMason** | Performs a structural multiple sequence alignment (MSA) on all conformers using the 3Di alphabet — a 20-state encoding of local backbone geometry. Produces a per-conformer lDDT score (structural self-consistency) |
| 4 | Conformer scoring | Each conformer scored as `0.6 × lDDT + 0.4 × epitope_structure`. Epitope score rewards conformers where PHF6\* (VQIINK), PHF6 (VQIVYK), and jR2R3 residues adopt non-disordered (non-D) 3Di states |
| 5 | Output | Top 5 conformers saved as PDB files + downloaded as zip |

### Epitope regions tracked
| Epitope | K18 residues | Weight | Why |
|---------|-------------|--------|-----|
| PHF6\* (VQIINK) | 2–7 | 1.0 | Primary aggregation nucleation site |
| PHF6 (VQIVYK) | 33–38 | 1.0 | Primary aggregation nucleation site |
| jR2R3 | 52–70 | 0.5 | R2/R3 junction, Sugiyama et al. 2025 |

### Output to save
```
tau_k18_conformer_results.zip
├── rank01_filtered_XXXX.pdb   ← upload to NB02b session 1
├── rank02_filtered_XXXX.pdb   ← upload to NB02b session 2
├── rank03_filtered_XXXX.pdb
├── rank04_filtered_XXXX.pdb
└── rank05_filtered_XXXX.pdb
```

---

## NB02b — BindCraft TauK18 Binder Design

**Goal**: Use structure-based hallucination to generate binders that physically
dock against a specific Tau K18 conformer, then verify and refine them.

**Run this notebook 5 times** — once per conformer PDB from NB01.

### Pipeline logic

```
Upload rank0X conformer PDB
        ↓
BindCraft  →  complex PDB (binder + Tau, docked)
        ↓
ProteinMPNN  →  redesigned binder sequences (8 per backbone)
        ↓
Protenix  →  re-predicts complex from sequence alone
        ↓
Epitope contact filter  →  passing binders
```

### Step-by-step explanation

**BindCraft**
Takes the Tau conformer PDB and hallucinate a binder from scratch using
AF2-multimer internally. You specify hotspot residues (PHF6\*, PHF6) and a
binder length range. BindCraft outputs a complex PDB containing both the Tau
backbone and the hallucinated binder backbone, with an ipTM score estimating
binding confidence. This is the creative step — it generates novel binder
geometries that fit the target surface.

**ProteinMPNN on the complex**
BindCraft's hallucinated binder backbone is structurally reasonable but the
sequence is not optimally designed for it. ProteinMPNN reads the full complex
backbone (both chains) and asks: *what amino acid sequence would pack most
stably onto this backbone at this interface?* Tau (chain A) is fixed — only the
binder (chain B) is redesigned. Because it sees the interface geometry, it
can place complementary residues (hydrophobic patches, salt bridges, H-bonds)
that BindCraft's hallucination may have missed. The output is a set of
optimized binder sequences — **not new coordinates, just sequences**.

**Why extract binder-only sequences if Protenix reassembles them?**
Once ProteinMPNN changes the sequence, the BindCraft coordinates are stale —
the side chains no longer match the backbone. The new sequence needs to be
evaluated on its own merits. Protenix is a blind test: given only the two
sequences (K18 + redesigned binder), does AF3-level structure prediction
predict a stable complex? If yes, the sequence is independently validated.
This avoids the risk of inflated scores from BindCraft's self-consistent
optimization.

**Protenix (AF3-level complex prediction)**
ByteDance's open-source AlphaFold3 implementation. Takes two protein sequences,
predicts the complex structure from scratch with no MSA, outputs a CIF file
and a confidence JSON containing:
- `iptm` — interface predicted TM-score (how confident the chains are in the
  right relative orientation)
- `ptm` — global predicted TM-score
- `ranking_score = 0.8 × iptm + 0.2 × ptm`

No coordinates are carried forward from BindCraft or ProteinMPNN. This is
deliberately a clean-room prediction.

**Epitope contact analysis**
Parses the Protenix CIF to find binder residues within 8 Å (Cα–Cα) of
PHF6\*, PHF6, or jR2R3 residues on K18. A binder that scores well on ipTM
but doesn't actually contact the epitopes isn't useful.

**Filter**: `ranking_score ≥ 0.5 AND epitope_contacts ≥ 1`

### What each score validates
| Score | Source | What it measures |
|-------|--------|-----------------|
| BindCraft ipTM | AF2-multimer (internal) | Does the hallucinated backbone bind at all? |
| ProteinMPNN score | Implicit in sequence | Is this sequence the best fit for this backbone + interface? |
| Protenix ranking_score | AF3-level, blind | Does the final sequence fold and bind from scratch? |
| Epitope contacts | CIF geometry | Does the binder actually touch the target epitopes? |

Each step strips away artifacts from the previous step and validates the thing
that matters downstream.

### Output to save (per conformer)
```
bindcraft_{conformer_name}.zip
├── bindcraft_{conformer_name}_passing.csv   ← critical — upload to NB03
├── cifs/*.cif                               ← optional for NB03 structural analysis
└── bindcraft_qc_{conformer_name}.png
```

---

## NB02c — Forge TauK18 Binder Design

**Goal**: Use sequence-only generative design to produce binders across multiple
lengths, independently of any specific conformer PDB.

**Run this notebook once.** All binder lengths are handled in a single session.

### Why Forge is different from BindCraft
Forge is a **sequence-to-sequence** model — it only needs the target protein
sequence as input. No PDB upload, no hotspot residues. It uses flow matching
on ESM2 embeddings to generate binder sequences that are likely to bind the
target based on learned sequence-level co-evolution patterns. This means it
can sample a different region of sequence space than structure-based methods.

### Pipeline logic

```
K18 sequence (hardcoded)
        ↓
Forge  →  300 binder sequences (100 × 3 lengths: 50, 65, 80 aa)
        ↓
ESMFold pass 1  →  backbone PDB per sequence + pLDDT filter (≥70)
        ↓
ProteinMPNN (monomer mode)  →  4 redesigned sequences per backbone
        ↓
free GPU (del esmfold)
        ↓
Protenix  →  binder + K18 complex, no MSA
        ↓
Epitope contact filter  →  passing binders
```

### Why ESMFold pass 1 is kept
Forge generates sequences only — no structure. ProteinMPNN needs a backbone
to work with. ESMFold converts each Forge sequence into a backbone PDB so
ProteinMPNN has something to optimize against. ESMFold is *not* used for
final validation — that's Protenix's job.

### Why ProteinMPNN runs in monomer mode here
Unlike BindCraft, Forge doesn't produce a complex PDB — it only produces
a binder backbone. There is no Tau chain present in the PDB, so ProteinMPNN
can't use interface geometry. It redesigns the binder sequence purely to be
self-consistent with the backbone. Protenix then evaluates whether the
redesigned sequence actually binds K18.

**ESMFold is unloaded before Protenix** (`del esmfold; torch.cuda.empty_cache()`)
because both models are large and the T4 doesn't have enough VRAM to hold both
simultaneously.

**Filter**: `ranking_score ≥ 0.5 AND epitope_contacts ≥ 1`

### Output to save
```
forge_designs.zip
├── forge_binder_filter_passing.csv   ← critical — upload to NB03
├── structures/*.cif
└── forge_pipeline_qc.png
```

---

## NB03 — Binder Comparison Analysis

**Goal**: Integrate all passing binders from NB02b (×5) and NB02c (×1),
compare them structurally and chemically, and produce a final ranked list.

### Inputs
- `bindcraft_rank01_passing.csv` through `bindcraft_rank05_passing.csv`
- `forge_binder_filter_passing.csv`
- Known control binders (AT8, HJ8.5, PHF1, Baker 2025 miniproteins)

### Pipeline
| Step | Tool | What it does |
|------|------|-------------|
| 1 | ESMFold | Folds any binders missing structures |
| 2 | FoldMason | Structural MSA on all binders — 3Di alphabet, lDDT |
| 3 | SaProt / ESM2 | Joint sequence+structure embeddings (1280-dim) |
| 4 | UMAP + HDBSCAN | Clusters binders by structural/chemical similarity |
| 5 | Contact maps | Checks epitope coverage per binder |
| 6 | Chemical validation | Checks for problematic residues, charge, solubility |
| 7 | Integrated score | Combines all metrics into a final ranking |

### Integrated score formula
```
0.30 × iptm
+ 0.20 × plddt
+ 0.15 × cluster_confidence
+ 0.20 × epitope_contact_score
+ 0.10 × chemical_score
+ 0.05 × known_binder_proximity
```

### Output
```
binder_comparison_results.zip
├── final_ranked_binders.csv
├── umap_clusters.png
├── epitope_coverage.png
└── integrated_scores.csv
```

---

## Data Flow Summary

```
NB01
  └─ rank01–05.pdb ──────────────────────────────┐
                                                  ↓
NB02b (×5)                              bindcraft_rankXX_passing.csv
                                                  │
NB02c (×1)                         forge_binder_filter_passing.csv
                                                  │
                                                  ↓
                                              NB03
                                      final_ranked_binders.csv
```

**Minimum files to keep between sessions:**

| After | Keep |
|-------|------|
| NB01 | `rank01–05_filtered_XXXX.pdb` |
| NB02b (each run) | `bindcraft_{conformer}_passing.csv` |
| NB02c | `forge_binder_filter_passing.csv` |

Everything else (intermediate PDBs, ProteinMPNN run dirs, Protenix run dirs)
lives in `/content/` and is lost when the Colab session ends.

---

## Tool Reference

| Tool | Role | Input type |
|------|------|-----------|
| Starling | IDP conformer generation | Sequence |
| AFRC | Polymer physics filter | Sequence |
| FoldMason | Structural MSA / 3Di alphabet | PDB ensemble |
| BindCraft | Structure-based binder hallucination | PDB + hotspots |
| Forge | Sequence-based binder generation | Sequence |
| ProteinMPNN | Sequence optimization on fixed backbone | PDB (backbone) |
| ESMFold | Fast structure prediction | Sequence |
| Protenix | AF3-level complex prediction (blind test) | Two sequences |
| SaProt | Joint seq+structure embeddings | Sequence + 3Di tokens |
| UMAP / HDBSCAN | Dimensionality reduction + clustering | Embeddings |
