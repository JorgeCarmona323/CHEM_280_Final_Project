# Phase 2: Binder Generation

Two parallel tracks: structure-based (Track A) and sequence-based (Track B).

---

## Track A: Structure-Based Binder Design

### A1 — RFDiffusion

RFDiffusion diffuses over protein backbone coordinates conditioned on a target hotspot region on the IDP structure.

**Input**: Relaxed target structure (PDB) from Phase 1
**Output**: ~100–500 backbone-only binder scaffolds per target

```bash
# RFDiffusion binder design
python /path/to/RFdiffusion/scripts/run_inference.py \
    inference.output_prefix=results/rfdiffusion/tau_monomer/binder \
    inference.input_pdb=data/structures/tau_relaxed/tau_monomer_relaxed.pdb \
    'contigmap.contigs=[B1-100/0 50-100]' \
    'ppi.hotspot_res=[B244,B268,B274,B280]' \
    inference.num_designs=100 \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0
```

**Hotspot residues**: Key residues in the IDP repeat domain known to be important for interactions (R1 hexapeptide: PHF6*, R2–R3 junction).

**Then design sequences with ProteinMPNN**:
```bash
python /path/to/ProteinMPNN/protein_mpnn_run.py \
    --pdb_path results/rfdiffusion/tau_monomer/ \
    --out_folder results/mpnn/tau_monomer/ \
    --num_seq_per_target 8 \
    --sampling_temp 0.1 \
    --batch_size 8
```

**Then fold with ESMFold and filter by ipTM > 0.7**:
```bash
python src/binder_generation/fold_and_filter.py \
    --input results/mpnn/tau_monomer/ \
    --output results/esmfold/tau_monomer/ \
    --iptm_threshold 0.7
```

---

### A2 — BindCraft (Alternative)

BindCraft performs end-to-end hallucination using AF2 multimer + Rosetta refinement. Better for short peptide binders and disordered targets where RFDiffusion hotspot-selection is ambiguous.

```bash
python /path/to/BindCraft/bindcraft.py \
    --target data/structures/tau_relaxed/tau_monomer_relaxed.pdb \
    --target_chains B \
    --binder_len 50-80 \
    --num_designs 50 \
    --output results/bindcraft/tau_monomer/
```

**When to prefer BindCraft over RFDiffusion**:
- IDP target has no clear hotspot (BindCraft optimizes binding site jointly)
- Smaller peptide binders desired (<80 aa)
- Want all-atom quality from the start (no separate MPNN step)

---

## Track B: Sequence-Based Design with Forge

Forge generates novel protein sequences via flow matching on Raygun (ESM2) embeddings. We run Forge on all construct sequences to generate binder-like sequences conditioned on IDP embedding context.

**Note**: Forge is not running locally (requires `yk0/forge-e80` checkpoint + GPU). Use Colab.

### Colab workflow
```python
# Install
!pip install transformers torch
!git clone https://github.com/yk0/forge

# Load model
from forge import ForgeModel
model = ForgeModel.from_pretrained("yk0/forge-e80")

# Generate for each construct sequence
constructs = {
    "tau_monomer": "QGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEG...",
    "tau_oligomer": "...",
    "tau_nft": "...",
    "targetA_construct": "...",
    "targetA_full": "...",
}

for name, seq in constructs.items():
    binders = model.generate_binders(
        target_sequence=seq,
        num_samples=100,
        temperature=1.0,
        length_range=(50, 100)
    )
    # Save to results/forge/{name}/
```

**Output**: 100 generated sequences per construct — these are then:
1. Folded with ESMFold
2. Filtered by pLDDT > 70 (well-folded)
3. Fed into Phase 3 motif analysis

---

## Combined Output for Phase 3

After both tracks, we have:

| Source | Constructs | Count |
|--------|-----------|-------|
| RFDiffusion + ProteinMPNN + ESMFold | 5 constructs × 100 designs × 8 seqs | ~4000 |
| BindCraft | 5 × 50 | 250 |
| Forge | 5 × 100 | 500 |

Filter to high-confidence designs (ipTM/pLDDT thresholds), targeting ~200–500 candidates total entering Phase 3.

---

## Per-Construct Design Notes

### Tau monomer
- Focus hotspots on: PHF6 (VQIVYK, residues 306–311), PHF6* (VQIINK, residues 275–280)
- These hexapeptides are the aggregation-driving cores and likely binding interfaces

### Tau oligomer
- Less defined structure — use BindCraft over RFDiffusion
- Design peptides that can engage partially folded β-sheet conformations

### Tau NFT (protofilament)
- Use cryo-EM structure (PDB 6QJH, AD Tau fold)
- Clear structural target — RFDiffusion with hotspots on filament outer surface
- Avoid the intramolecular core (already packed) — target solvent-exposed grooves

### Target protein constructs
- TBD once sequences are provided
