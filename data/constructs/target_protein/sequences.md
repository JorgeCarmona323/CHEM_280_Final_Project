# Target Protein A: FMRP (Fragile X Mental Retardation Protein)

**Gene**: FMR1 | **UniProt**: Q06787 (human, 621 aa isoform used here)
**Project role**: "Fragile protein × protein that has an IDP region"

FMRP is an RNA-binding protein with a well-folded N-terminal domain (NLS, Tudor, KH1, KH2)
and a long disordered C-terminal region containing the RGG box — a phase-separation
and aggregation driver implicated in Fragile X syndrome.

---

## Disorder Profile Summary (IUPred2A, long-range)

| Region          | Residues | IUPred score      | Interpretation                         |
|-----------------|----------|-------------------|----------------------------------------|
| N-terminal / Tudor | 1–100 | 0.02–0.25       | Highly structured                      |
| KH1 domain      | 100–270  | 0.05–0.40         | Structured RNA-binding domain          |
| KH2 domain      | 270–380  | 0.10–0.50         | Structured (some flexible loops)       |
| Transition zone | 380–425  | 0.50–0.76         | Semi-disordered; avoid for Starling    |
| **IDR / RGG box** | **425–621** | **0.76–0.91** | **→ Starling input region**        |
| ANCHOR peaks    | 460–475, 610–621 | high ANCHOR | Binding-competent sub-regions     |

**Recommended Starling trim: residues 430–621 (~191 aa)**

Rationale:
- Cleanly past the structured KH2 domain
- Avoids the 380–425 transition zone (semi-structured, confuses Starling)
- ~191 aa is under the 200 aa default max_length
- Captures the full RGG box + C-terminal disordered tail
- ANCHOR peaks at 460–475 and 610–621 indicate two binding-competent sub-regions

---

## Constructs

### A1: FMRP IDR Construct (Starling input)

```
>fmrp_idr_430_621
```
*(Residues 430–621 of Q06787; verify exact sequence from UniProt before use)*

Run IUPred parse:
```bash
python src/structure_prep/parse_iupred.py \
    --input data/iupred/fmrp_iupred.txt \
    --threshold 0.5 \
    --padding 15 \
    --max_length 200 \
    --fasta_out data/constructs/target_protein/fmrp_idr_trimmed.fasta \
    --plot
```

### A2: FMRP Full Protein (stable domain anchor)

Full 621 aa sequence — use AlphaFold2/ColabFold to predict stable N-terminal structure.
The AF2 output (structured KH1/KH2 domains, aa 1–425) serves as the anchor for stitching.

```bash
# Run ColabFold on full sequence
colabfold_batch data/constructs/target_protein/fmrp_full.fasta \
    data/structures/fmrp_af2/ --num-recycle 3
```

After AF2: stitch the AF2 stable domain (aa 1–429) + Starling IDR ensemble (aa 430–621).

---

## Phase 1 Pipeline (FMRP)

```bash
# Step 1: Generate IDR ensemble with Starling
starling <fmrp_430_621_sequence> --outname fmrp_idr -c 400 -r

# Step 2: Extract individual frames
python src/structure_prep/extract_starling_frames.py \
    --topology fmrp_idr_STARLING.pdb \
    --trajectory fmrp_idr_STARLING.xtc \
    --output data/structures/fmrp_idr_ensemble/ \
    --n_frames 100

# Step 3: FoldMason ensemble QC
python src/structure_prep/foldmason_refine.py \
    --ensemble data/structures/fmrp_idr_ensemble/ \
    --output data/structures/fmrp_idr_refined/ \
    --top_n 5 \
    --tmp_dir tmp/foldmason_fmrp_idr/ \
    --report results/structure_prep/fmrp_idr_conformer_ranking.csv

# Step 4: Stitch IDR to AF2 stable domain
python src/structure_prep/stitch_constructs.py \
    --stable data/structures/fmrp_af2/fmrp_full.pdb \
    --idp_ensemble data/structures/fmrp_idr_refined/ \
    --idp_residues 430-621 \
    --output data/structures/fmrp_stitched/

# Step 5: FoldMason pre-MD quality gate
python src/structure_prep/foldmason_refine.py \
    --ensemble data/structures/fmrp_stitched/ \
    --output data/structures/fmrp_stitched_refined/ \
    --top_n 3 \
    --tmp_dir tmp/foldmason_fmrp_stitched/ \
    --report results/structure_prep/fmrp_stitched_conformer_ranking.csv

# Step 6: MD relaxation
python src/structure_prep/run_md_relaxation.py \
    --input data/structures/fmrp_stitched/fmrp_stitched.pdb \
    --output data/structures/fmrp_relaxed/fmrp_construct_relaxed.pdb \
    --idp_residues 430-621

# Step 7: FoldSeek identity validation
python src/structure_prep/validate_construct.py \
    --construct data/structures/fmrp_relaxed/fmrp_construct_relaxed.pdb \
    --reference_name "FMRP" \
    --foldseek_db pdb \
    --output results/validation/construct_validation/ \
    --min_tm 0.5
```

---

## Biology Notes

- **RGG box** (~aa 527–552 in 2N4R-style numbering): arginine-glycine-glycine repeats; mediates
  RNA binding, protein-protein interaction, liquid-liquid phase separation
- **Phase separation**: FMRP C-terminal IDR drives condensate formation; relevance to
  synaptic protein translation in neurons
- **Disease link**: Fragile X syndrome = CGG repeat expansion in FMR1 → FMRP loss-of-function
- **Binding partners**: NUFIP1, CYFIP1/2, FXR1/2; the IDR mediates many of these interactions
- **Therapeutic angle**: binders targeting the FMRP IDR could modulate phase separation or
  restore protein-protein interactions lost in FX syndrome
- **Known structures**: AF2 model available; KH1/KH2 co-crystal with RNA (PDB: 2QND); full
  C-terminal IDR is unresolved experimentally

---

## Comparison with Tau Constructs

| Feature              | Tau MTBR (B1)         | FMRP IDR (A1)            |
|----------------------|-----------------------|--------------------------|
| Region               | 583–720 (Big Tau)     | 430–621                  |
| Length               | ~137 aa               | ~191 aa                  |
| Disorder type        | Aggregation-prone IDP | Phase-separating IDR     |
| Pathology            | Neurodegeneration     | Fragile X syndrome       |
| Key sub-motifs       | PHF6* (VQIINK), PHF6 (VQIVYK) | RGG box, ANCHOR regions  |
| Starling ready?      | Yes (trimmed MTBR)    | Yes (trimmed aa 430–621) |
