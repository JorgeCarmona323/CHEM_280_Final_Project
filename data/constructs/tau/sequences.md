# Tau Construct Sequences

## Tau Isoforms
Tau has 6 isoforms from the MAPT gene (0N3R, 1N3R, 2N3R, 0N4R, 1N4R, 2N4R).
The most studied is 2N4R (441 aa, full-length) and 0N4R (383 aa).

## Constructs for this project

### B1: Tau Repeat Domain (K18 fragment, residues 244–372, 0N4R numbering)
- Contains R1–R4 microtubule-binding repeats
- Includes the aggregation-driving hexapeptides PHF6* (VQIINK) and PHF6 (VQIVYK)
- This is the minimal aggregation-competent fragment
- Source: UniProt P10636-8 (2N4R), residues 244–372

```
>tau_monomer_K18
KVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL
```
*(Verify against UniProt P10636 before use)*

### B2: Tau Toxic Oligomer
- Same sequence as B1 (K18), but structure generated under aggregation-promoting conditions
- Use Starling with β-sheet biased conformations, or use published oligomer cryo-EM seed (if available)
- Key references: Lasagna-Reeves et al. 2012 (FASEB J), Bhatt et al. 2020

### B3: Tau NFT (Protofilament core, Alzheimer's disease fold)
- AD Tau protofilament core: residues 274–380 (C2 fold)
- Source PDB: **6QJH** (straight filament), **5O3L** (Fitzpatrick et al. 2017)
- Use experimental cryo-EM structure directly — do NOT predict from sequence

```
>tau_nft_core_AD
SKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYK
```
*(Protofilament core only; residues 274-380 approx)*

## Target Protein A
To be added once sequences are confirmed with collaborator.

## Notes
- All sequences should be verified against UniProt before generating structures
- Use canonical 2N4R numbering for residue positions throughout the project
- Phosphorylation sites: S202, T205 (AT8 epitope), S396, S404 (PHF1 epitope)
