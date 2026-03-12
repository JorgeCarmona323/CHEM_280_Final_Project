# Phase 1: Structure Preparation

## Goal
Produce one or more representative relaxed structures per construct that capture the conformational heterogeneity of the IDP region while maintaining the correct topology of any stable domains.

---

## 1.1 Input Sequences

### Tau constructs
| Name | Region | Notes |
|------|--------|-------|
| `tau_monomer` | Full-length 2N4R (441 aa) or 0N4R repeat domain (aa 244–368) | Start with repeat domain |
| `tau_oligomer` | Same repeat domain, but MD initialized as pre-aggregated seed | Use REMD or biased MD |
| `tau_nft` | Protofilament core (aa 274–380, C2 polymorph) | Use cryo-EM PDB as input (PDB: 5O3L, 6QJH) |

### Target protein constructs
| Name | Region | Notes |
|------|--------|-------|
| `targetA_construct` | IDP region / disordered loop only | Add once sequences confirmed |
| `targetA_full` | Full protein | Use AlphaFold2 for stable domains |

---

## 1.2 Tools

### Starling — IDP ensemble generation
Starling (Baker lab) generates structural ensembles for disordered proteins using a flow-based model trained on coarse-grained IDP simulations.

```bash
# Install (Python 3.10+, requires GPU)
pip install starling-protein

# Generate ensemble for Tau repeat domain
python -c "
from starling import Ensemble
ens = Ensemble.from_fasta('data/constructs/tau/tau_monomer.fasta')
ens.generate(n_structures=100, output_dir='data/structures/tau_monomer_ensemble/')
"
```

**Output**: 100 PDB structures representing the disordered ensemble.
**For our purpose**: cluster by RMSD and take 3–5 representative conformers per construct.

### AlphaFold2 — stable domain prediction
For the full-protein constructs, predict the structured regions with AF2 and use those as the anchor for stitching.

```bash
# ColabFold (easiest)
# Upload FASTA to https://colabfold.mmseqs.com/
# Or use local ColabFold:
colabfold_batch data/constructs/tau/tau_full.fasta data/structures/tau_af2/ --num-recycle 3
```

---

## 1.3 Stitching IDP regions to stable domains

For constructs where the IDP region is flanked by structured elements (e.g., Tau proline-rich region + microtubule-binding repeat + C-terminal tail):

```bash
# Script: src/structure_prep/stitch_constructs.py
python src/structure_prep/stitch_constructs.py \
    --stable data/structures/tau_af2/tau_full.pdb \
    --idp_ensemble data/structures/tau_monomer_ensemble/ \
    --idp_residues 244-368 \
    --output data/structures/tau_stitched/
```

**Method**:
1. Extract stable N-terminal region (AF2 structure) up to IDP junction
2. Sample from Starling ensemble for the IDP region
3. Align termini using Cα superposition at junction residues ±5
4. Rebuild backbone connectivity with a short energy minimization

---

## 1.4 MD Relaxation (OpenMM)

After stitching, the junction region may have clashes or non-physical geometry. Run a short MD relaxation.

```bash
python src/structure_prep/run_md_relaxation.py \
    --input data/structures/tau_stitched/tau_monomer_stitched.pdb \
    --output data/structures/tau_relaxed/tau_monomer_relaxed.pdb \
    --steps 50000 \          # 100 ps at 2 fs timestep
    --temperature 300 \      # Kelvin
    --forcefield amber14
```

**Protocol**:
- Force field: AMBER14-SB + TIP3P water
- Energy minimize → NVT equilibration (100 ps) → NPT production (100–500 ps)
- Restrain Cα of stable regions, let IDP region relax freely
- Extract lowest-energy frame as representative structure

---

## 1.5 Output Structures (per construct)

| File | Description |
|------|-------------|
| `tau_monomer_relaxed.pdb` | Tau repeat domain, extended monomer |
| `tau_oligomer_relaxed.pdb` | Pre-aggregated seed state |
| `tau_nft_relaxed.pdb` | NFT protofilament conformation |
| `targetA_construct_relaxed.pdb` | Collaborator construct |
| `targetA_full_relaxed.pdb` | Full target protein |

These 5 structures feed into Phase 2 (binder generation).

---

## Notes / Caveats
- Tau NFT structure: use published cryo-EM models (PDB 5O3L or 7P65) directly rather than predicting from scratch — these are experimentally validated
- Toxic oligomers are poorly characterized structurally; consider using multiple Starling conformers biased toward β-sheet-prone configurations
- Starling requires GPU; run on Colab if not available locally
