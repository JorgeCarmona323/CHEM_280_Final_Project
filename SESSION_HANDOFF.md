# CHEM 280 Session Handoff
**Date**: 2026-03-16  |  **Deadline**: Wednesday 2026-03-17

---

## Experimental Design (2 × 3 Factorial)

| | Monomer (STARLING rank01) | Oligomer (AF2-multimer) | Fibril (6HRE) |
|---|---|---|---|
| **BindCraft** | partner runs | partner runs | partner runs |
| **Forge** | CSV in hand | need to run | need to run |

Top 3 by Protenix ranking score per cell → up to 18 candidate binders → NB03 comparison

---

## What Was Done This Session

### NB02c Forge Generation (`colab/NB02c_Forge/02c_Forge_TauK18_Binder_Design.ipynb`)
- Added `TAU_STATE` config variable (`monomer` / `oligomer` / `fibril`)
- Three target sequences defined in config
- Outputs tagged `forge_{tau_state}_top{N}.csv` + `tau_state` column added
- Stripped ~16 junk cells leftover from previous Colab run

### NB02c Downstream (`colab/NB02c_Forge/02c_Forge_downstream_ESMFold_MPNN_Protenix.ipynb`)
- Same `TAU_STATE` + three target sequences added to config
- Output CSV renamed to `forge_{TAU_STATE}_top{TOP_N}.csv`
- `tau_state` column added to output for NB03

### NB03 — NOT YET UPDATED
- Still expects 2 CSVs (old design), needs full update next session

---

## What Needs To Be Done Next Session

### 1. Update NB03 (code)
Config cell should use a dict of 6 CSVs:
```python
CONDITION_CSVS = {
    ('Forge',     'monomer'):  'forge_monomer_top10.csv',
    ('Forge',     'oligomer'): 'forge_oligomer_top10.csv',
    ('Forge',     'fibril'):   'forge_fibril_top10.csv',
    ('BindCraft', 'monomer'):  'bindcraft_monomer_passing.csv',
    ('BindCraft', 'oligomer'): 'bindcraft_oligomer_passing.csv',
    ('BindCraft', 'fibril'):   'bindcraft_fibril_passing.csv',
}
```
Also add: `tau_state` column to load, UMAP tau_state coloring panel,
epitope heatmap faceted by tau_state, 2x3 summary heatmap plot.

### 2. Colab runs needed

**ColabFold K18 homodimer** (run first — needed for BindCraft oligomer)
- https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb
- Sequence: paste K18 twice separated by `:` (K18:K18)
- Settings: `num_models=1`, `num_recycles=3`, `model_type=multimer`
- Save best PDB as `k18_dimer.pdb` to Drive: `CHEM_280/results/nb02b_output/`

**NB02c oligomer run** — set `TAU_STATE = 'oligomer'`, run full notebook
**NB02c fibril run** — set `TAU_STATE = 'fibril'`, run full notebook

**Partner BindCraft** (3 sessions):
- monomer: use STARLING rank01 PDB (already in Drive)
- oligomer: use `k18_dimer.pdb`
- fibril: notebook auto-downloads 6HRE

**NB03** — runs last after all 6 CSVs collected

---

## Target Sequences Quick Reference

```
Monomer (168 aa):
KVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGG
QVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPV
VSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL

Oligomer (79 aa):
KVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKGGG
KVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYK

Fibril (38 aa):
KVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYK
```

---

## Key File Locations

| Item | Path |
|------|------|
| Forge monomer CSV (done) | `results/forge_sequences_clean.csv` (local, gitignored) |
| STARLING conformers | `data/constructs/tau/tau_k18_phase1_results_final.zip` |
| Forge generation NB | `colab/NB02c_Forge/02c_Forge_TauK18_Binder_Design.ipynb` |
| Forge downstream NB | `colab/NB02c_Forge/02c_Forge_downstream_ESMFold_MPNN_Protenix.ipynb` |
| Analysis NB | `colab/NB03_BinderComparison/03_Binder_Comparison_Analysis.ipynb` |
| Drive base | `/content/drive/MyDrive/CHEM_280/results/` |

---

## Git Status
All NB02c changes committed. NB03 update pending next session.
