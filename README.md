# CHEM 280 Final Project: IDP Binder Generation Pipeline

> **Generative design of binders for intrinsically disordered proteins and disordered regions of structured proteins**

---

## Scientific Question

Can we systematically discover binders for IDP regions — including conformationally heterogeneous states like monomers, toxic oligomers, and aggregated fibrils — by combining structural ensemble generation, generative binder design, and embedding-based motif discovery?

---

## Target System

### Protein A — [Fragile X Related Protein 1 (human)]
| Construct | Description |
|-----------|-------------|
| A1 | IDP/disordered construct only |
| A2 | Full protein (IDP region + stable domains) |

### Protein B — Tau (MAPT)
| Construct | Description |
|-----------|-------------|
| B1 | Tau monomer |
| B2 | Toxic oligomer |
| B3 | Neurofibrillary tangle (NFT, paired helical filament) |

All Tau constructs focus on the repeat domain (R1–R4, residues ~244–368) which contains the disordered aggregation-prone microtubule-binding region.

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 1: STRUCTURE PREPARATION                                  │
│                                                                  │
│  Sequences (FASTA)                                               │
│       │                                                          │
│       ▼                                                          │
│  Starling (IDP ensemble)  ──▶  MD relaxation (OpenMM)           │
│       │                              │                           │
│       ▼                              ▼                           │
│  Stitch IDP regions ◀──── Stable region structures (AlphaFold2)  │
│       │                                                          │
│       ▼                                                          │
│  Representative structures (one per construct, per state)        │
└───────────────────────────┬─────────────────────────────────────┘
                            │
          ┌─────────────────┴──────────────────┐
          │                                    │
          ▼                                    ▼
┌─────────────────────┐            ┌───────────────────────┐
│  PHASE 2A:          │            │  PHASE 2B:            │
│  Structure-based    │            │  Sequence-based       │
│  Binder Design      │            │  Binder Design        │
│                     │            │                       │
│  RFDiffusion or     │            │  Forge                │
│  BindCraft          │            │  (flow matching on    │
│  (on structures)    │            │  Raygun embeddings)   │
└──────────┬──────────┘            └──────────┬────────────┘
           │                                  │
           └────────────────┬─────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 3: MOTIF ANALYSIS                                         │
│                                                                  │
│  All generated binder structures/sequences                       │
│       │                                                          │
│       ├──▶ FoldMason (structural alignment → enriched motifs)    │
│       │                                                          │
│       ├──▶ ESM2 single-residue embeddings → UMAP → HDBSCAN       │
│       │    (embedding clusters → sequence/structure trends)       │
│       │                                                          │
│       └──▶ FoldSeek (structural similarity → known binder hits)  │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 4: VALIDATION                                             │
│                                                                  │
│  Positive controls (known binders):                              │
│    • Hemoglobin / EGFR (structured protein, ground truth)        │
│    • Baker lab IDP binder benchmarks (2025)                      │
│                                                                  │
│  IDP-specific validation:                                        │
│    • Published Tau binders (antibodies: AT8, PHF1, MC1)          │
│    • Failed clinical trial compounds (LMTM, TRx0237)            │
│    • α-syn, FUS binder literature                                │
│                                                                  │
│  Metrics:                                                        │
│    • Embedding cluster overlap with known binders                │
│    • FoldSeek structural recall vs. validation set               │
│    • ipTM / pLDDT (ESMFold) on top candidates                    │
└─────────────────────────────────────────────────────────────────┘
```

---

## Tool Stack

| Tool | Purpose | How used |
|------|---------|----------|
| **Starling** | IDP ensemble generation | Generate disordered monomer conformers |
| **AlphaFold2** | Stable domain structures | Anchor the IDP stitch |
| **OpenMM** | MD relaxation | Relax stitched IDP+stable structures |
| **RFDiffusion** | Structure-based binder design | Design binders on relaxed IDP structures |
| **BindCraft** | All-atom binder design | Alternative to RFDiffusion |
| **Forge** | Sequence-based design | Flow matching on Raygun embeddings |
| **ESMFold** | Structure prediction | Fold generated binder sequences |
| **ProteinMPNN** | Sequence design | Optimize sequences on top scaffolds |
| **FoldMason** | Structural MSA + motifs | Find enriched structural motifs across binders |
| **FoldSeek** | Structural search | Validate motifs against PDB/known binders |
| **ESM2** | Residue embeddings | Single-residue embedding clusters |
| **UMAP + HDBSCAN** | Dimensionality reduction + clustering | Discover embedding trends |

---

## Repository Structure

```
CHEM_280_Final_Project/
├── README.md
├── workflow/
│   ├── 01_structure_prep.md       ← Starling + stitching + MD protocol
│   ├── 02_binder_generation.md    ← RFDiffusion/BindCraft/Forge protocol
│   ├── 03_motif_analysis.md       ← FoldMason + embedding pipeline
│   └── 04_validation.md          ← controls + benchmark protocol
├── src/
│   ├── structure_prep/
│   │   ├── stitch_constructs.py   ← stitch IDP to stable region PDBs
│   │   └── run_md_relaxation.py   ← OpenMM relaxation wrapper
│   ├── binder_generation/
│   │   ├── run_rfdiffusion.py     ← RFDiffusion config + launcher
│   │   └── run_forge.py           ← Forge inference wrapper
│   ├── motif_analysis/
│   │   ├── motif_scanner.py       ← ESM2 embedding + UMAP + clustering
│   │   ├── foldmason_parser.py    ← parse FoldMason structural MSA output
│   │   └── foldseek_utils.py      ← FoldSeek search + hit parsing
│   └── validation/
│       ├── benchmark_binders.py   ← score candidates vs. validation set
│       └── known_binder_db.py     ← curated known binder/drug database
├── data/
│   ├── constructs/
│   │   ├── tau/                   ← Tau FASTA sequences + PDB inputs
│   │   └── target_protein/        ← Collaborator protein sequences
│   ├── known_binders/             ← Validation set: known binder structures
│   ├── failed_drugs/              ← Clinical trial compound structures
│   └── structures/                ← Generated/relaxed structures (gitignored)
├── notebooks/
│   ├── 01_structure_prep.ipynb
│   ├── 02_binder_generation.ipynb
│   ├── 03_motif_analysis.ipynb
│   └── 04_validation.ipynb
├── results/                       ← Output files (gitignored)
└── requirements.txt
```

---

## Validation Strategy

### Tier 1: Structured protein control
- **Hemoglobin (HBA1/HBB)** or **EGFR kinase domain** — well-characterized binders exist
- Confirms the pipeline produces sensible binders for stable proteins
- Baseline for ipTM and FoldSeek recall metrics

### Tier 2: Published IDP binder controls (Baker lab 2025)
- de novo designed IDP binders from:
  - *Nature* 2025: generalized IDP binder design
  - *Science* 2025: specific IDP/disordered target
- Use as positive class in embedding cluster validation

### Tier 3: Tau-specific validation
- **Antibodies**: AT8 (pS202/pT205), PHF1 (pS396/pS404), MC1 (conformational)
- **Failed drugs**: LMTM/TRx0237 (Tau aggregation inhibitor, Phase 3 fail)
- **Cryo-EM structures**: PDB entries for Tau filament conformations

---

## Key Scientific Insight

> IDPs don't have a single binding-competent conformation. By generating binders across **multiple states** (monomer → oligomer → fibril), we can ask:
> - Do the same structural/sequence motifs recur across states?
> - Are state-specific binders discoverable?
> - Do embedding clusters separate by disease-relevant aggregation state?

This is the core hypothesis that distinguishes this project from single-state binder design.

---

## References

See `workflow/04_validation.md` for full citation list.
Key papers driving the design:
- Baker lab IDP binder design (Nature 2025, Science 2025)
- RFDiffusion (Watson et al., Nature 2023)
- BindCraft (Pacesa et al., 2024)
- Forge (Young et al., unpublished)
- FoldMason (van Kempen et al., 2024)
- FoldSeek (van Kempen et al., Nature Biotechnology 2024)
- ESM2 (Lin et al., Science 2023)
