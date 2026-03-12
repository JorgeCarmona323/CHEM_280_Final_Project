# Phase 4: Validation

## Goal
Determine which generated binders are scientifically credible by:
1. Verifying the pipeline works on well-characterized systems (controls)
2. Checking whether generated IDP binders recapitulate known binding modes
3. Using failed clinical trial compounds as a negative/positive benchmark

---

## 4.1 Positive Controls — Structured Proteins

Use a well-characterized protein-binder pair to confirm the full pipeline (Phase 1–3) produces sensible results.

### Option A: EGFR kinase domain
| Item | Detail |
|------|--------|
| Target PDB | 1IVO (EGFR kinase), 2ITY (EGFR + Erbitux Fab) |
| Known binders | Erbitux (cetuximab) epitope, Afatinib binding pocket |
| Test | Can RFDiffusion/BindCraft recover binders that dock to the known Erbitux epitope? |
| ESMFold metric | ipTM > 0.8 for top candidates vs. known crystal pose RMSD |

### Option B: Hemoglobin
| Item | Detail |
|------|--------|
| Target PDB | 1HHO (deoxy-Hb), 1A3N (oxy-Hb) |
| Known binders | DARPins, nanobodies from literature; 2-3-DPG binding site |
| Test | Structural/embedding recall of known binders vs. generated set |

**Pass criterion**: FoldSeek recall of known binders ≥ 10% in top-50 ranked candidates.

---

## 4.2 IDP Binder Controls — Baker Lab 2025

Two landmark papers establishing the feasibility of de novo IDP binder design.

---

### Paper A — Nature 2025
> **"Diffusing protein binders to intrinsically disordered proteins"**
> Liu, Wu, Choi, Han, Baker et al. *Nature* 644, 809–820 (2025). DOI: 10.1038/s41586-025-09248-9

**Method**: Flexible target fine-tuned RFdiffusion starting from **sequence only**. Key innovation: **two-sided partial diffusion** — both target and binder conformations sampled simultaneously.

**Key targets and Kd values**:
| Target | Type | Kd (nM) | Notes |
|--------|------|---------|-------|
| Amylin (hIAPP) | IDP | 3.8–100 | Four binders; different conformations; crystal structures |
| C-peptide (31 aa) | IDP | 28 | CP-35; thermostable 95°C |
| VP48 | IDP | 39 | Transcription activator |
| BRCA1_ARATH IDR | IDP | 52 | 21-residue disordered segment |
| FUS (239–267) | IDP | 520 | Challenging; 52% Gly/Ser/Asn/Gln; 3/94 bound |
| G3BP1 RBD (IDR) | IDR-in-structured | 11 | **Crystal structure**; disrupts stress granules in cells |
| Prion (180–187) | IDR-in-structured | 14 | VNITIKQH β-strand |
| IL-2RG (327–336) | IDR-in-structured | 97 | Cell surface receptor IDR |

**FUS relevance**: FUS is directly analogous to Tau — both are disease-relevant RNA-binding IDPs with low-complexity disordered regions. The 520 nM result (3/94 binding) shows highly disordered low-complexity sequences remain challenging. Tau is similar — calibrate expectations accordingly.

**Best positive control for our pipeline**: G3BP1-11 — crystal structure, cellular validation, high affinity.

---

### Paper B — Science 2025
> **"Design of intrinsically disordered region binding proteins"**
> Wu, Jiang, Hicks, Baker et al. *Science* 389, eadr8063 (2025). DOI: 10.1126/science.adr8063

**Method**: Template library threading — 1000 repeat-protein templates covering diverse extended conformations → thread target IDR subsequences → refine top matches with RFdiffusion + ProteinMPNN.

**Scale**: 39/43 targets yielded binders; 100 pM–100 nM; ~22 designs/target.

**Selected targets**:
| Target | Kd (nM) | Relevance |
|--------|---------|-----------|
| Dynorphin A-r2 | <0.08 | Best-in-class control |
| Mesothelin IDR | <1 | IDR in structured protein |
| Angiotensin I | 75 | Short IDP |
| EWS/FLI fusion | >500 | Harder IDP target (similar difficulty to Tau) |
| Telomerase IDR | 91 | Cellular IDP |

**Template threading relevance**: Identify unique 8–40 aa subsequences in Tau constructs (PHF6, PHF6*, R1–R4 junctions) → thread through Baker template library → use matches as RFdiffusion seeds.

---

### Additional paper — BindCraft
> **"One-shot design of functional protein binders with BindCraft"**
> Pacesa, Nickel, Schellhaas, Correia et al. *Nature* 646, 483 (2025). DOI: 10.1038/s41586-025-09429-6

AF2 multimer hallucination + MPNNsol + AF2 monomer filtering. 10–100% success rates. Targets: PD-1 (<1 nM), CD45 (15 nM), SpCas9 (267 nM). Used in our Track A for structured/IDR targets.

---

### How to obtain structures
- **G3BP1-11, amylin-22, amylin-18**: Crystal structures deposited — search PDB for DOI 10.1038/s41586-025-09248-9
- **Science paper binders**: Sequences in Supplementary Tables → ESMFold → add to `data/known_binders/`

---

## 4.3 Tau-Specific Validation

### Published antibody binders
| Antibody | Epitope | PDB | Notes |
|----------|---------|-----|-------|
| AT8 | pS202/pT205 | — | Diagnostic marker; phospho-specific |
| PHF1 | pS396/pS404 | — | C-terminal phospho |
| MC1 | Conformational (aa 7-9 + 313-322) | — | Disease-specific conformation |
| HJ8.5 | N-terminal (aa 25-30) | 6NWP | Crystal structure available |
| DC8E8 | MTBR repeat region | — | Blocks Tau-Tau interaction |

**Validation approach**:
- Take FoldSeek hits against antibody paratope structures (VH/VL CDR regions)
- Check if our generated binders dock to similar epitope regions
- Score epitope overlap: what fraction of epitope residues contact the generated binder?

### Failed clinical trial compounds
| Drug | Mechanism | Trial | Outcome |
|------|----------|-------|---------|
| LMTM (TRx0237) | Tau aggregation inhibitor (methylene blue) | Phase 3 | Failed 2016 |
| Semorinemab | Anti-Tau antibody (N-terminal) | Phase 2 | Failed 2021 |
| Gosuranemab | Anti-Tau antibody (N-terminal, aa 6-23) | Phase 2 | Failed 2021 |
| Tilavonemab | Anti-Tau antibody (mid-region) | Phase 2 | Failed 2022 |
| BIIB076 | Anti-Tau antibody (repeat domain) | Phase 1b | Halted |

**Use as**:
1. **Negative control interpretation**: if our binders cluster with these, they target same (possibly ineffective) epitopes — flag but don't automatically discard
2. **Positive binding evidence**: drug-bound Tau structures (if available) confirm these regions ARE bindable
3. **Structural input**: LMTM binding mode can serve as a reference for small-molecule-proximal binding sites

---

## 4.4 Validation Metrics Summary

| Metric | Source | Threshold | Meaning |
|--------|--------|-----------|---------|
| ipTM | ESMFold | > 0.7 | Well-folded, plausible interface |
| pLDDT (binder) | ESMFold | > 70 | Binder itself is well-folded |
| FoldSeek recall | Known binder DB | ≥ 10% in top-50 | Pipeline recovers known binding modes |
| HDBSCAN cluster overlap | Baker 2025 binders | Any overlap | Embedding space consistent with known solutions |
| Epitope overlap | Tau antibody epitopes | ≥ 3 shared contacts | Targeting known binding surfaces |
| FoldMason conservation | Cross-state motifs | ≥ 70% conserved | Motif is state-agnostic (more generalizable) |

---

## 4.5 Interpretation Framework

```
High ipTM + FoldSeek recall + cluster overlap with Baker 2025
    → Strong candidate: targets a known-bindable region with a novel scaffold

High ipTM + clusters with failed drugs (e.g., LMTM region)
    → Caution: targets a region with clinical precedent but prior failure
      (explore why drug failed: delivery? specificity? mechanism?)

High ipTM + no FoldSeek hits + novel cluster
    → Novel binder: potentially interesting but unvalidated binding mode
      (prioritize for further in silico docking or MD simulation)

Low ipTM regardless
    → Deprioritize (likely not a well-folded, stable binder)
```

---

## References

1. Watson JL et al. De novo design of protein structure and function with RFdiffusion. *Nature* 2023.
2. Pacesa M et al. BindCraft: one-shot design of functional protein binders. *bioRxiv* 2024.
3. Lin Z et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 2023.
4. van Kempen M et al. Fast and accurate protein structure search with Foldseek. *Nature Biotechnology* 2024.
5. Baker lab IDP binder paper 1. *Nature* 2025. [ADD DOI after reading]
6. Baker lab IDP binder paper 2. *Science* 2025. [ADD DOI after reading]
7. Fitzpatrick AWP et al. Cryo-EM structures of tau filaments from Alzheimer's disease. *Nature* 2017.
8. Falcon B et al. Structures of filaments from Pick's disease reveal a novel tau conformation. *Nature* 2018.
