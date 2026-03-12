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

## 1.3 FoldMason Ensemble Refinement

Before stitching, filter the Starling ensemble to keep only the most structurally consistent conformers. This removes outliers (misfolded or unrealistic structures) and ensures the representative conformers used downstream reflect the true ensemble average.

```bash
python src/structure_prep/foldmason_refine.py \
    --ensemble data/structures/tau_monomer_ensemble/ \
    --output data/structures/tau_monomer_refined/ \
    --top_n 5 \
    --tmp_dir tmp/foldmason_tau_monomer/ \
    --report results/structure_prep/tau_monomer_conformer_ranking.csv
```

**How it works**:
1. Runs `foldmason easy-msa` on all ensemble PDBs
2. Parses the structural MSA to compute a **consensus** character per alignment column
3. Scores each conformer by how closely it matches the column-wise consensus (**consistency score**)
4. Combines with FoldMason per-structure **lDDT** (structural quality) into a single rank
5. Copies the top-N ranked conformers to the refined output directory

**Scoring**:
```
combined_score = 0.6 × consistency + 0.4 × lDDT_normalized
```

Adjust weights with `--consistency_weight` / `--lddt_weight` as needed.

**Why this matters for IDPs**:
Starling ensembles can include low-probability, physically implausible conformers. FoldMason refinement ensures the stitching and binder design steps operate on a clean, representative subset rather than outliers that would produce poor binders.

**Output**: `conformer_ranking.csv` lists every conformer with its rank, consistency score, and lDDT — useful for understanding ensemble diversity.

---

## 1.4 Stitching IDP regions to stable domains

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

## 1.5 FoldMason Pre-MD Quality Gate

Before committing to MD (which is expensive), run FoldMason on the stitched construct(s) to confirm the IDP motifs survived stitching and the overall structure is coherent.

**The question**: *Are the structurally conserved motifs identified in step 1.3 still present and well-formed after stitching?*

If the FoldMason score drops sharply relative to the ensemble → the stitching disrupted the IDP motif → try the next-ranked conformer before spending GPU time on MD.

```bash
# Run FoldMason on all stitched conformers together
python src/structure_prep/foldmason_refine.py \
    --ensemble data/structures/tau_stitched/ \
    --output data/structures/tau_stitched_refined/ \
    --top_n 3 \
    --tmp_dir tmp/foldmason_tau_stitched/ \
    --report results/structure_prep/tau_stitched_conformer_ranking.csv
```

**What to check in `conformer_ranking.csv`**:
- Consistency scores should be close to the ensemble scores from step 1.3
- A large drop (>0.2) in consistency score means the IDP motif is distorted at the junction
- Only pass top-scoring stitched conformers to MD

**Second FoldMason vs first FoldMason** — purpose differs:
| Step | Input | Purpose |
|------|-------|---------|
| 1.3 (first) | Raw Starling ensemble | Pick representative IDP conformers |
| 1.5 (second) | Stitched constructs | Verify motifs survived stitching |

---

## 1.6 MD Relaxation (OpenMM)

After the pre-MD FoldMason quality gate passes, run MD on the top-ranked stitched conformers.

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

## 1.7 Construct Validation via FoldSeek Identity Check

After MD relaxation, validate that the Frankenstein construct is still recognizable as the intended protein. This confirms the stitched IDP conformer is structurally plausible in context — not a phantom structure that drifted away from the real protein.

**The question**: *Does FoldSeek retrieve the correct protein when searching the relaxed construct against PDB?*

If yes → the IDP conformer is realistic and fits coherently with the stable domain.
If no → the chosen IDP conformer is incompatible; try a different Starling ensemble member.

```bash
python src/structure_prep/validate_construct.py \
    --construct data/structures/tau_relaxed/tau_monomer_relaxed.pdb \
    --reference_name "Tau" \
    --foldseek_db pdb \
    --output results/validation/construct_validation/ \
    --min_tm 0.5
```

**How it works**:
1. FoldSeek searches the relaxed stitched structure against PDB
2. Top hits are parsed and filtered by TM-score (≥ 0.5)
3. Checks whether `reference_name` (e.g. "Tau", "MAPT") appears in top-K hits
4. Writes a summary report: confirmed ✓ or failed ✗

**Pass criterion**: TM-score ≥ 0.5 to a known Tau/MAPT structure in top-10 hits.

**If a construct fails**:
- Try a different Starling conformer (next in the `conformer_ranking.csv`)
- Widen the junction window in `stitch_constructs.py` (--junction_window 8)
- Check whether the IDP region FoldMason motif is genuinely compatible with the stable domain termini

**Full loop** (Starling → FoldMason → stitch → MD → FoldSeek → iterate):
```bash
for RANK in 01 02 03 04 05; do
    python src/structure_prep/stitch_constructs.py \
        --stable data/structures/tau_af2/tau_full.pdb \
        --idp_ensemble data/structures/tau_monomer_refined/rank${RANK}_*.pdb \
        --idp_residues 244-368 \
        --output data/structures/tau_stitched/

    python src/structure_prep/run_md_relaxation.py \
        --input data/structures/tau_stitched/stitched_rank${RANK}_*.pdb \
        --output data/structures/tau_relaxed/tau_monomer_rank${RANK}_relaxed.pdb \
        --idp_residues 244-368

    python src/structure_prep/validate_construct.py \
        --construct data/structures/tau_relaxed/tau_monomer_rank${RANK}_relaxed.pdb \
        --reference_name "Tau" \
        --output results/validation/construct_validation/ && break
done
# Stops at first conformer that passes validation
```

---

## 1.8 Output Structures (per construct)

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
