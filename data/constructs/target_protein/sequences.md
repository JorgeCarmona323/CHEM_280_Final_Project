# Target Protein A: FMRP (Fragile X Mental Retardation Protein)

**Gene**: FMR1 | **UniProt**: Q06787 (human)
**Project role**: "Fragile protein × protein that has an IDP region"

FMRP is an RNA-binding protein. The N-terminal half is well-folded (structured domains),
while the C-terminal region is a long intrinsically disordered region (IDR) containing
the RGG box — a key phase-separation and aggregation driver.

---

## Domain Architecture (InterPro / NCBI CDD)

| Domain | InterPro / CDD accession | Residues (approx) | Notes |
|--------|--------------------------|--------------------|-------|
| SMART/Tudor | IPR000870 / SM00333 | ~1–70 | Chromatin/RNA reader |
| KH domain 1 | IPR004087 / cd00105 | ~200–270 | RNA-binding |
| KH domain 2 | IPR004087 / cd00105 | ~270–340 | RNA-binding |
| **C-terminal IDR** | *(no annotated domain)* | **~420–621** | **→ Starling input** |

Verify exact residue boundaries using:
- InterPro: https://www.ebi.ac.uk/interpro/protein/UniProt/Q06787/
- NCBI CDD: https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi

**IDP construct boundary logic**: the Starling input starts where the last annotated domain
(KH2) ends and the unstructured C-terminal IDR begins. Validate against:
1. PDB structures — KH2 is resolved; C-terminal IDR has no experimental structure
2. AF2 model — pLDDT drops below 70 at the KH2/IDR boundary
3. Literature — RGG box and C-terminal IDR well-established as disordered (Darnell et al.)

---

## Constructs

### A1: FMRP IDR Construct (Starling input)

Disordered C-terminal IDR — boundaries set by InterPro KH2 domain end + literature.

```
>fmrp_idr
```
*(Paste residues from KH2 domain end to C-terminus; verify exact coordinates from InterPro
before use — do not rely on estimated numbers above)*

Starling run:
```bash
starling <fmrp_idr_sequence> --outname fmrp_idr -c 400 -r
```

### A2: FMRP Full Protein (stable domain anchor)

Full sequence for AF2/ColabFold prediction. The structured N-terminal half (Tudor + KH1 + KH2)
is the anchor for stitching after Starling generates the IDR ensemble.

```bash
colabfold_batch data/constructs/target_protein/fmrp_full.fasta \
    data/structures/fmrp_af2/ --num-recycle 3
```

---

## Phase 1 Pipeline (FMRP-specific bash)

```bash
# Step 1: Starling on IDR construct
starling <fmrp_idr_sequence> --outname fmrp_idr -c 400 -r

# Step 2: Extract frames
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

# Step 4: Stitch to AF2 stable domain
# --idp_residues: use actual residue numbers from InterPro KH2 end
python src/structure_prep/stitch_constructs.py \
    --stable data/structures/fmrp_af2/fmrp_full.pdb \
    --idp_ensemble data/structures/fmrp_idr_refined/ \
    --idp_residues <START>-621 \
    --output data/structures/fmrp_stitched/

# Steps 5–7: FoldMason pre-MD QC → MD relax → FoldSeek identity check
# (see workflow/01_structure_prep.md §1.5–1.7 for full loop)
```

---

## Biology Notes

- **RGG box**: arginine-glycine-glycine repeats in the IDR; mediates RNA binding, LLPS,
  and protein-protein interactions with partners including NUFIP1, CYFIP1/2
- **Disease**: Fragile X syndrome — CGG repeat expansion → FMRP loss-of-function;
  disrupted synaptic translation regulation
- **Therapeutic angle**: binders to the FMRP IDR could modulate phase separation
  or rescue protein-protein interactions lost in FX syndrome
- **Known PDB structures**: KH1/KH2 co-crystal with RNA (PDB 2QND); C-terminal IDR
  has no resolved experimental structure
