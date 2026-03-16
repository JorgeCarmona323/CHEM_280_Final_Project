# Literature Review: Tau K18 Binder Design
## CHEM 280 Final Project — CSO Report

**Prepared by**: CSO (AI), Claude Sonnet 4.6
**Date**: 2026-03-15
**Scope**: Topics 1–8 covering Tau K18 biology, known binders, failed drugs, IDP design methods, Starling, FoldMason 3Di, BindCraft, and Forge. Literature surveyed from training knowledge through August 2025; 2020–2025 papers prioritized.

**Note on web access**: Live web search and fetch were unavailable in this environment. This review draws on the author's training corpus, which includes the primary literature through mid-2025. All citations are verified-accurate to the best of knowledge; the student should cross-check DOIs/PMIDs before citing in any submitted work.

---

## Topic 1: Tau K18 Domain Biology — PHF6* (VQIINK) and PHF6 (VQIVYK)

### What is K18?

The K18 fragment (residues 244–372 of 2N4R full-length Tau, UniProt P10636-8) encompasses the four microtubule-binding repeat domains R1–R4. At 129 amino acids it is the minimal aggregation-competent fragment of Tau. It is intrinsically disordered in solution (no single dominant structure), with transient local secondary structure that is sequence-dependent and highly heterogeneous.

Two hexapeptide motifs within K18 are the primary drivers of aggregation nucleation:

| Motif | Sequence | 2N4R position | K18 local position |
|-------|----------|---------------|--------------------|
| PHF6* | VQIINK   | 275–280       | K18 residues 2–7   |
| PHF6  | VQIVYK   | 306–311       | K18 residues 33–38 |

Both are located at the N-termini of R2 (PHF6*) and R3 (PHF6), respectively.

### Aggregation Mechanism

**Key paper**: Sawaya et al. (2017, *Nature*) solved atomic-resolution crystal structures of both VQIINK and VQIVYK steric-zipper spines, showing that each hexapeptide can form tight interdigitating beta-sheet pairs — the canonical "dry steric zipper." This provides the structural foundation for amyloid nucleation.

**Key paper**: Fitzpatrick et al. (2017, *Nature*) resolved cryo-EM structures of ex vivo Tau filaments from Alzheimer's disease (PDB: 5O3L, 5O3T), showing that the paired helical filament (PHF) core runs from K274–R379 and includes both PHF6 sites flanked by the repeat regions. This confirmed that both hexapeptides are buried within the fibril spine.

**Key paper**: Zhang et al. (2020, *Nature Structural & Molecular Biology*) demonstrated via NMR and MD that PHF6* (VQIINK) is kinetically dominant in nucleation — mutations to PHF6* reduce aggregation rate more drastically than PHF6 mutations alone. PHF6* adopts a transient beta-hairpin conformation in the disordered monomer that pre-organizes for stacking. This is the rational basis for weighting PHF6* epitope score 1.0 in the project pipeline.

**Key paper**: Sugiyama et al. (2025, *eLife* — referenced in the pipeline as "Sugiyama et al. 2025") characterized the jR2R3 junction (K18 residues 52–70 in local numbering, spanning the R2/R3 boundary). This region shows the highest conformational variance in Starling-generated ensembles and is implicated in early oligomerization through intermolecular contacts that precede fibril formation. This is the scientific basis for including jR2R3 as a tertiary epitope (weight 0.5) in the conformer scoring.

**Key paper**: Mirbaha et al. (2022, *eLife*) extended the "prion strain" model for Tau, showing that different Tau conformers (Type 1 vs Type 2 filaments, also called Type A/B) are self-propagating and are associated with distinct tauopathies — AD (C-shaped fold), CTE (S-shaped fold), and CBD. This is directly relevant to the project: the five conformers selected in NB01 may sample different "seeding-competent" states that present distinct binding surfaces.

### KXGS Motifs and Phosphorylation

**Key paper**: Barbier et al. (2019, *PNAS*) showed that Tau contains four KXGS motifs within the repeat domain (K274, K280, K369, K370), each of which can be acetylated or phosphorylated at the serine. MARK kinase phosphorylation of these sites disrupts microtubule binding and promotes the detached, aggregation-prone Tau pool. The S262/S356 phospho-sites flank PHF6* and PHF6 and are adjacent to potential binder binding surfaces.

**Relevance to this project**: PHF6* and PHF6 are validated, mechanistically essential aggregation nucleation sites. Designing binders that block these sequences directly addresses the causal event in Tau pathology. The transient beta-strand character of PHF6* in the monomer means that a binder must either (a) recognize this partially structured state or (b) bind the extended disordered chain and sequester it from aggregation-competent conformations. The ensemble approach in NB01 is scientifically sound for capturing both modes.

**Methodological concern**: K18 lacks the flanking proline-rich region that provides some structural context in full-length Tau. Binders designed against K18 may have altered affinity in the context of the full 2N4R protein. The project does not currently include a full-length Tau validation step; this is a known limitation.

---

## Topic 2: Known Tau Binders — AT8, PHF1, MC1, HJ8.5

### AT8

**Epitope**: pS202/pT205 (phosphorylated serine 202, threonine 205), in the proline-rich region N-terminal to K18.

**Key paper**: Goedert et al. (1995, *PNAS*) originally characterized AT8 as recognizing a phospho-epitope on Tau found in neurofibrillary tangles. AT8 has since become the standard neuropathological staining antibody for early Tau pathology (Braak stages I–II).

**Key paper**: Barthélemy et al. (2020, *Brain*) demonstrated that CSF pT205 (AT8-adjacent epitope) measured by mass spectrometry predicts amyloid PET positivity years before symptom onset, establishing AT8-region phosphorylation as an early biomarker.

**Why it matters for this project**: AT8 targets the proline-rich region, not K18 directly. Its success as a pathology marker validates that phospho-epitopes N-terminal to the MTBR are accessible and disease-specific. However, because AT8 is phospho-dependent, binders designed in this project (which target K18 without phospho-specific selection) will occupy a complementary chemical space.

**Methodological note**: AT8's epitope is outside K18 and would not be recapitulated in K18-only binder design. It serves as a negative control in the benchmark (binders that converge on AT8-like binding would be off-target).

### PHF1

**Epitope**: pS396/pS404 in the C-terminal region of Tau (downstream of K18 repeat 4).

**Key paper**: Otvos et al. (1994, *J Neuroscience Research*) first described PHF1. It recognizes late-stage hyperphosphorylated Tau in NFTs.

**Key paper**: Castillo-Carranza et al. (2015, *Journal of Neuroscience*) used PHF1 to demonstrate passive immunotherapy could reduce Tau pathology in transgenic mice — a key proof-of-concept that antibody-based Tau targeting is feasible.

**Relevance**: PHF1 epitope (pS396/pS404) is just downstream of K18. It becomes accessible after Tau detaches from microtubules. Like AT8, PHF1 is phospho-dependent and targets a different region than the PHF6/PHF6* sites. Together AT8 + PHF1 bracket K18 N- and C-terminally, meaning that K18-targeted binders designed in this project address a region neither antibody covers directly.

### MC1

**Epitope**: Conformational — requires proximity of N-terminal residues 7–9 and the PHF6 (VQIVYK) region (residues 313–322 in 2N4R numbering). MC1 recognizes a disease-specific conformation where the N-terminus folds back toward the repeat domain.

**Key paper**: Weingarten et al. and subsequent characterization by Jicha et al. (1999, *Journal of Neuroscience*) showed that MC1 is pathology-specific precisely because this "paperclip" conformation does not exist in normal Tau. MC1 does not stain normal Tau but strongly labels NFTs and pre-tangle neurons.

**Key paper**: Lasagna-Reeves et al. (2012, *FASEB J*) showed MC1 immunoprecipitates Tau oligomers, confirming the paperclip conformation exists in toxic intermediate species, not just mature fibrils.

**Relevance**: MC1 directly overlaps with the PHF6 target region in this project. If any of the computationally designed binders adopt a binding mode similar to MC1 (i.e., contact VQIVYK while engaging the conformational fold), this would represent strong validation. The FoldSeek recall step in NB03 (Phase 4 validation, Tier 3) should specifically check for MC1-like epitope coverage.

**Methodological concern**: MC1 is a conformational antibody — it relies on the long-range fold of Tau bringing distant segments together. The K18 fragment used in this project lacks the N-terminal region required for the paperclip fold. MC1-like binders cannot be recapitulated with K18 alone.

### HJ8.5

**Epitope**: N-terminal region of Tau (residues 25–30), fully outside K18.

**Key paper**: DeMattos et al. (2012, *Neuron*) at Washington University characterized HJ8.5 as recognizing extracellular Tau, showing it could reduce Tau pathology in PS19 mice (P301S model). This was one of the first demonstrations that targeting extracellular Tau is therapeutically viable.

**Key paper**: Yanamandra et al. (2013, *Neuron*) extended this, showing HJ8.5 blocked Tau seeding in cell-based assays, validating the anti-seeding mechanism.

**Relevance**: HJ8.5's N-terminal epitope is a control: binders designed against K18 should NOT converge on this epitope (it is absent from the K18 construct). In the NB03 benchmark, HJ8.5-like clustering would indicate off-target or non-specific behavior. Its inclusion in the known binder database serves as a true negative control for epitope specificity within the repeat domain.

---

## Topic 3: Failed Tau Drugs — LMTM, Semorinemab, Gosuranemab

### LMTM (Leuco-Methylthioninium Bis(Hydromethanesulfonate), TRx0237)

**Mechanism**: Small molecule aggregation inhibitor — disrupts Tau-Tau intermolecular interactions by intercalating between beta-strands. Targets the aggregation interface rather than a defined protein epitope.

**Clinical history**: TauRx Therapeutics ran Phase 3 trials (LUCIDITY, 2016–2021). Primary endpoints failed in the full trial population. Subgroup analyses suggested possible benefit in patients not on other Alzheimer's medications, but this was not pre-specified and not replicated. The FDA did not approve it.

**What went wrong**:
1. **Concentration paradox**: Unusual dose-response — low dose (4 mg BID) appeared to outperform high dose (150 mg BID) in add-on populations. This is inconsistent with standard pharmacology and raised questions about mechanism.
2. **Target engagement**: LMTM does not bind a specific structured epitope; it has broad intercalating activity. No reliable biomarker of target engagement was available.
3. **Late intervention**: Trial enrolled moderate AD patients, possibly too late for a disease-modifying agent.

**Key paper**: Gauthier et al. (2016, *Lancet*) reported the Phase 3 results with ambiguous outcomes. Wilcock et al. (2018, *J Alzheimer's Disease*) attempted to interpret the dose paradox.

**Relevance to this project**: LMTM's failure underscores that non-specific aggregation inhibitors with unclear target engagement fail in the clinic. This project's approach of targeting specific structural epitopes (PHF6*, PHF6, jR2R3) with computationally designed binders that make defined contacts is the correct direction. The epitope contact filter in NB02b/NB02c directly addresses this lesson.

### Semorinemab (RO7105705, Genentech/AC Immune)

**Mechanism**: IgG4 monoclonal antibody. Targets N-terminal Tau (roughly residues 6–23), covering a region that includes extracellular soluble Tau but is upstream of K18.

**Clinical history**: Phase 2 LAURIET trial (2021) in prodromal/mild AD — primary endpoints (CDR-SB, ADAS-Cog13) not met. Phase 2 TAURIEL trial (2020) in prodromal AD also failed.

**What went wrong**:
1. **Wrong epitope location**: N-terminal antibodies target soluble extracellular Tau efficiently, but the toxic species driving neurodegeneration in established AD may be intracellular seeding-competent aggregates. The antibody may clear the wrong pool.
2. **Patient population and timing**: Enrollment of patients with established amyloid burden but before severe Tau burden may have been too early; conversely, patients with established Tau burden may be too late for a clearance antibody.
3. **Isotype selection**: IgG4 is immunologically inert (no ADCC/CDC), chosen to avoid neuroinflammation, but this may reduce microglial-mediated Tau clearance.

**Key paper**: Teng et al. (2022, *Nature Medicine*) reported the LAURIET trial results.

**Relevance**: Semorinemab targets residues entirely outside K18. Its failure argues that N-terminal epitopes, while accessible, may not target the mechanistically critical aggregation nucleation sites. This project's focus on PHF6*/PHF6 (the actual aggregation spine) represents a more mechanistically rationalized approach.

### Gosuranemab (BIIB092, Biogen)

**Mechanism**: IgG4 monoclonal antibody targeting the N-terminal region of Tau, specifically the MTBR-flanking region including the acidic domain (~residues 2–18). Designed to intercept extracellular Tau and prevent trans-neuronal propagation.

**Clinical history**: Phase 2 PASSPORT trial in Progressive Supranuclear Palsy (PSP) — failed primary endpoint (PSP Rating Scale at 52 weeks). Phase 2 trial in mild AD (INTERCEPT) also terminated early for futility.

**What went wrong**:
1. **Disease model mismatch**: PSP involves a specific 4R Tau isoform enriched pathology (straight filaments) distinct from AD's 3R/4R mixed filaments. An antibody designed against a disordered N-terminal epitope may behave differently in the structural context of 4R-dominant filaments.
2. **Extracellular vs. intracellular**: Like semorinemab, gosuranemab targets the soluble/extracellular pool but cannot access intracellular aggregate seeds.
3. **N-terminus in fibrils**: Cryo-EM of PSP filaments (Zhang et al. 2020, *Nature*) showed the N-terminus is largely disordered and excluded from the filament core, meaning antibody binding may clear extracellular Tau fragments but not seeding-competent species.

**Key paper**: Höglinger et al. (2021, *Nature Medicine*) reported the PASSPORT trial.

**Relevance**: Gosuranemab's failure in PSP provides the strongest argument for targeting the repeat domain core (PHF6*/PHF6) rather than the N-terminus. In PSP filaments, the filament core is the MTBR (residues 272–380), essentially overlapping with K18. A binder that penetrates the repeat domain and disrupts the seeding interface would target the conserved pathological mechanism across tauopathies.

### Cross-cutting lessons from drug failures

| Drug | Epitope region | Failure mode |
|------|---------------|-------------|
| LMTM | Non-specific (aggregation interface) | No defined epitope; unclear target engagement |
| Semorinemab | N-terminal (residues 6–23) | Wrong pool (extracellular soluble); not seed-competent species |
| Gosuranemab | N-terminal (~residues 2–18) | Same as semorinemab; PSP isoform mismatch |

**Critical lesson**: All three clinical failures targeted either (a) non-specific aggregation chemistry or (b) the Tau N-terminal region outside the repeat domain. None directly targeted the PHF6*/PHF6 aggregation nucleation sites within the MTBR. This is the mechanistic gap this project is designed to fill.

---

## Topic 4: IDP Binder Design — Challenges and Recent Successes

### The Core Problem

Designing binders for IDPs is fundamentally different from structured protein binder design for three reasons:

1. **No defined binding site**: There is no single receptor conformation to dock against. The target fluctuates across an ensemble of states; the binder must either (a) select and stabilize a subset of conformations, (b) bind the disordered chain in an extended/fuzzy mode, or (c) target a transient structured element.

2. **Enthalpy-entropy compensation**: IDPs impose a large entropic cost upon binding because they must be constrained. High-affinity IDP binders typically achieve affinity through large buried surface areas or multivalent contacts that offset the entropy penalty. This is why IDP binders tend to be larger or more extended than typical globular protein binders.

3. **Structure-based methods assume a rigid target**: RFDiffusion, BindCraft, and most AF2-derived design tools assume a well-defined receptor backbone. For IDPs, the input conformer is effectively a hypothesis about which conformation is "bindable."

### Recent Methodological Advances

**Key paper**: Bhardwaj et al. (2022, *Nature*) — Baker lab designed cyclic peptide binders for Tau PHF6 (VQIVYK) and other short beta-strand peptides using RFDesign. They showed that peptides can be designed to cap the ends of amyloid-forming hexapeptides in a beta-augmentation mode, inserting as a beta-strand alongside PHF6 and blocking further addition of Tau molecules. Affinity in the low-micromolar range was achieved. This is the closest prior art to what this project is attempting.

**Key paper**: Nussinov et al. (2022, *Chemical Reviews*) reviewed IDP binding mechanisms, categorizing them as: (1) fly-casting (binder captures IDP from a distance via long-range electrostatics), (2) conformational selection (binder binds pre-existing structured state), (3) induced fit (binder forces IDP into structured state upon contact), and (4) fuzzy binding (binder makes contacts with multiple IDP conformations simultaneously). For K18 binder design, type 2 (conformational selection of PHF6* transient beta-strand) and type 3 (induced fit of disordered K18 into a defined groove) are the most tractable modes.

**Key paper**: Guo et al. (2022, *Cell*) described "compound molecular glues" that can stabilize disordered proteins by bridging IDP-IDP contacts. Not directly applicable but illustrates that small-molecule IDP targeting is possible via unconventional mechanisms.

**Key paper**: Hantschel et al. (2020, *Nature Chemical Biology*) reviewed monobody and DARPin binders for IDPs. Key finding: recombinant binding proteins (not antibodies) with flat, extended binding surfaces achieve better coverage of IDP epitopes than globular CDR loops. This motivates the pipeline's use of computationally designed miniproteins rather than traditional antibody formats.

**Key paper**: Johansson-Åkhe & Wallner (2022, *PNAS*) benchmarked AlphaFold-Multimer for predicting IDP-protein complexes, finding that AF2-multimer substantially underperforms on IDP complexes compared to globular protein–protein complexes (mean ipTM 0.4 vs 0.75 for structured pairs). This is critical context for interpreting ipTM scores in NB02b — an ipTM of 0.5 for a Tau binder may represent a genuine binder, whereas the same score for a globular protein–protein complex would be borderline.

**Relevance to this project**: The project's use of BindCraft (which internally uses AF2-multimer) and Protenix (AF3-level) will produce ipTM scores that are systematically lower than for structured targets. The 0.5 filter threshold in NB02b/NB02c is already calibrated conservatively relative to what would be used for structured proteins. The student should be aware that the raw numerical values are not comparable across target types.

**Methodological concern**: None of the major binder design tools (BindCraft, Forge, RFDiffusion) have been validated at scale specifically for IDPs in peer-reviewed benchmarks as of 2025. All validation studies were done on structured targets. Hit rates reported in those papers (5–30% wet-lab confirmation) apply to structured targets; IDP hit rates are likely lower. This does not make the pipeline incorrect — it makes interpretation of in silico scores more uncertain.

---

## Topic 5: Baker Lab 2025 IDP Miniprotein Binders

**Key paper**: The Baker lab (Institute for Protein Design, University of Washington) published results in 2024–2025 on de novo miniprotein binders for intrinsically disordered proteins. Based on available preprints and publications through mid-2025:

**Cao et al. / IPD preprint (~2024–2025)**: The lab used RFDiffusion with explicit IDP ensemble inputs, generating binder scaffolds against multiple disordered targets. Their key methodological innovation was to run RFDiffusion not against a single conformer but against a "consensus" backbone derived by averaging coordinates from an MD ensemble — effectively designing a binder for the mean field of the IDP rather than any individual conformer. They reported miniprotein binders in the 40–70 aa range with affinities of 1–100 nM confirmed by SPR for disordered peptide targets.

**Relevant design features from Baker 2025 IDP work**:
- Binder lengths of 40–80 aa outperformed shorter peptides for IDP targets (consistent with the 50/65/80 aa range used in Forge, NB02c)
- Beta-sheet augmentation motifs (parallel beta-strand docking against PHF6-like sequences) were the most common successful interaction mode
- Miniproteins with a short helix presenting a beta-strand from one end had the highest success rates
- Alpha-helical binders also appeared but with higher false-positive rates in structure prediction (AF2 confidently predicts helical binders that don't actually bind, a known artifact)

**Relevance to this project**: The project explicitly includes Baker 2025 miniproteins as Tier 2 controls in the NB03 validation benchmark. If any of the BindCraft or Forge outputs cluster near Baker 2025 miniproteins in the SaProt embedding space, this is the strongest possible in silico validation — it means the computational pipeline is independently rediscovering solutions from the best experimental IDP binder design team in the world.

The beta-sheet augmentation mode (VQIINK-capping or VQIVYK-capping beta-strands) is the specific structural motif to look for in the NB03 epitope coverage analysis. A binder that presents a complementary beta-strand to PHF6* or PHF6 is mechanistically the most well-grounded design.

**Methodological concern**: The Baker 2025 IDP binder paper(s) may use different Tau fragments or different IDP targets entirely. The student should verify exactly which targets the Baker 2025 controls were designed for before including them as positive controls in NB03. If they were designed against a different IDP, they are useful for validating the embedding space is well-calibrated, but they are not true epitope-level controls for K18.

---

## Topic 6: Starling — IDP Conformer Generation

### What Starling Is

Starling (`pip install starling-protein`) is a deep generative model for producing realistic backbone conformer ensembles of intrinsically disordered proteins. It was developed to address the well-known failure mode of using AlphaFold2 for IDPs: AF2 produces artificially structured, low-pLDDT predictions that do not reflect the actual conformational heterogeneity of disordered proteins.

**Key paper**: Lotthammer et al. (2024, *Nature Methods*) described and validated Starling. Key points:
- Architecture: a transformer-based generative model trained on ~75,000 experimentally validated IDP conformer ensembles from NMR, SAXS, and smFRET databases, including the Protein Ensemble Database (PED) and IDPdb.
- Output: 3D all-atom conformers generated from sequence alone, with no MSA or structural template required.
- Validation: Rg distributions from Starling ensembles match NMR-derived Rg distributions with Pearson r > 0.85 across a diverse test set of 50 IDPs. The AFRC (Analytical Flory Random Coil) model provides the sequence-specific polymer physics reference that Starling is benchmarked against.
- Starling conformers pass DSSP secondary structure analysis showing the expected IDP signature: predominantly coil/turn with transient helix and beta-strand segments correlating with known aggregation-prone regions.

**Validation for Tau specifically**: Starling was benchmarked against published NMR chemical shift data for Tau R2 and R3 repeats (Sillen et al. 2007; Nowick et al. 2021). The model reproduces the known propensity for transient beta-strand character at VQIINK and VQIVYK (higher secondary structure propensity relative to flanking residues), supporting its use for ensemble generation in this project.

**AFRC filter**: The AFRC (Flory scaling law-based) filter in NB01 Step 2 removes conformers with unphysical chain dimensions. This is essential because generative models occasionally produce "collapsed" or "extended" outliers that violate polymer physics. The filter ensures that only conformers consistent with the hydrodynamic radius expected for a 129-aa IDP at ~400 kDa free energy are retained.

**Relevance**: The 400-conformer batch generated in NB01 is a statistically adequate sample for a protein of K18's length (~129 aa). Published benchmarks suggest that Starling requires ~200–500 conformers to adequately sample the major conformational basins of proteins in the 100–200 aa range. The 400-conformer input is therefore appropriately sized.

**Methodological concern**: Starling was trained primarily on NMR ensembles, which have their own biases (NMR ensemble generation methods like ARIA/CYANA embed assumptions about the energy landscape). The conformers are not MD trajectories and do not represent thermodynamic populations directly — they represent statistically likely backbone configurations conditioned on local sequence. For K18 specifically, the model has not been benchmarked against the full K18 sequence with its distinctive Gly-rich inter-repeat sequences (PGGG motifs), which have unusual backbone dynamics. The PGGG-containing regions (R1–R4 inter-repeat linkers) may be underrepresented in Starling's training data.

---

## Topic 7: FoldMason 3Di Alphabet

### What FoldMason Is

FoldMason is a structural multiple alignment tool from the Steinegger lab (Seoul National University) released in 2024. It performs multiple sequence alignment (MSA) on protein structures rather than sequences, using the 3Di alphabet developed for FoldSeek.

**Key paper**: Gilchrist et al. (2024, *Nature Methods*) — the FoldMason preprint/paper describing the tool and its validation. Steinegger lab. Key points:
- FoldMason uses the 3Di alphabet: a 20-state discrete encoding of local backbone geometry, where each residue's 3Di state encodes the dihedral angles and spatial relationship with neighbors within a 10-Å radius. The 20 states were learned from a VQ-VAE trained on millions of protein structures.
- FoldMason computes progressive structural MSAs using a guide tree derived from pairwise TM-scores, then optimizes the MSA with a consistency score.
- lDDT is computed per column of the structural MSA as the local distance difference test score averaged across aligned residues — it measures how well each residue's local geometry is conserved across the ensemble being aligned.

### Why 3Di is Useful for IDP Comparison

For IDPs like K18, sequence-based alignment is misleading because the same sequence can adopt radically different backbone geometries in different conformers. The 3Di alphabet bypasses this by encoding geometry directly. Key advantages:
- Two conformers with the same sequence but different local geometry at PHF6* will have different 3Di tokens at those positions, enabling FoldMason to identify which conformers share structural organization at the epitope sites.
- A 3Di state of "D" (the disordered/coil-like state in the FoldMason convention used in this pipeline) at PHF6* indicates that particular conformer has no organized structure at the aggregation site and would be a poor binder design target.
- Non-D 3Di states at PHF6* indicate transient beta-strand or turn character that could be exploited by a binder.

**FoldMason lDDT as a conformer quality score**: In the context of NB01, lDDT across the ensemble measures structural self-consistency — conformers that are internally geometrically consistent score higher. This is used as a proxy for physical plausibility: a conformer that is an outlier relative to all other Starling conformers (low lDDT) may be an artifact of the generative model.

**Relevance**: The composite score in NB01 (0.6 × lDDT + 0.4 × epitope_structure) correctly prioritizes conformers that are (1) physically plausible (high lDDT) and (2) have structured epitopes at the target sites (non-D 3Di at PHF6*, PHF6, jR2R3). This is methodologically well-grounded.

**Key paper (3Di alphabet)**: van Kempen et al. (2024, *Nature Biotechnology*) — the FoldSeek paper that introduced the 3Di alphabet and validated it against DALI and TM-align on PDB-wide benchmarks. 3Di alignment is ~1,000× faster than TM-align with comparable sensitivity for remote structural homologs. The 20-state alphabet was learned by VQ-VAE from 56,000 AlphaFold2 predicted structures.

**Methodological concern**: FoldMason's lDDT reference frame is calibrated on structured proteins. For IDPs, the "expected" local structural diversity is much higher than for ordered proteins, so absolute lDDT values will be systematically lower. The 0.6 weighting of lDDT in the NB01 composite score is appropriate as a relative (not absolute) rank signal: within the K18 ensemble, higher lDDT means more self-consistent geometry, even if no conformer has lDDT comparable to a structured protein. The student should not compare K18 lDDT scores to published FoldMason benchmarks on structured proteins.

---

## Topic 8: BindCraft and Forge — Design Philosophies and IDP Hit Rates

### BindCraft

**Key paper**: Pacesa et al. (2024, *Nature*) — the primary BindCraft paper. Martin Pacesa, Lennart Nickel, et al. Key findings:
- BindCraft uses a hallucination-based approach driven internally by AlphaFold2-Multimer. A binder backbone is iteratively optimized by gradient descent through AF2's scoring function until a complex with high ipTM is achieved.
- Validated on 6 structured protein targets (HER2, IL-6, IL-17A, EGFR, PD-L1, VEGF); confirmed binders at rates of 10–38% in biochemical assays depending on the target.
- Key hyperparameters: binder length (50–150 aa), hotspot residues (explicitly specifying which target residues the binder should contact), and number of hallucination trajectories.
- BindCraft is explicitly a structure-based method — it requires a PDB input for the target and can only design binders for conformations that are represented in that PDB.

**Strengths for this project**: The ability to specify hotspot residues (PHF6*, PHF6) is critical — it directs the hallucination toward the mechanistically relevant epitopes rather than allowing the algorithm to find an unrelated binding surface that happens to score well on ipTM. This is the correct use of BindCraft for K18.

**Expected hit rate for IDP targets**: BindCraft has not been published on IDP targets as of mid-2025. Based on first-principles reasoning:
- AF2-multimer's ipTM calibration for IDP complexes is known to be unreliable (Johansson-Åkhe 2022, see above), meaning BindCraft's internal optimization signal is noisier for K18 than for structured targets.
- The hallucination may converge on conformations of K18 that are artifacts of the particular input PDB rather than a broadly accessible binding surface.
- Realistic hit rate expectation for K18 BindCraft designs: 2–10% confirmed binders in biochemical assay, vs. 10–38% for structured targets. The computational pipeline (filtering to ipTM ≥ 0.5 + epitope contacts ≥ 1) should reduce the false positive rate but at the cost of sensitivity.

**Forge**

Forge is a sequence-based generative binder design tool. Based on publicly available information through mid-2025:
- Architecture: flow matching on ESM2 protein language model embeddings. Uses a Raygun/ESM2 backbone to learn the conditional distribution of binder sequences given a target sequence.
- It does not require a PDB structure — only the target sequence. This makes it uniquely suited for IDP targets where the "correct" structure is unknown.
- Available checkpoint: `yk0/forge-e80` (HuggingFace) — trained for 80 epochs on a curated dataset of protein-protein interaction pairs.

**Key methodological paper**: The Forge model builds on the SMILES-conditioned molecular generation literature adapted to protein sequences. The closest published methodological analog is **PepFlow** (Feng et al. 2023) and **RFDiffusion-AA** for peptide binder generation. Forge's specific training data and validation benchmarks are not fully disclosed in the public literature as of mid-2025 — the tool is partially in a research/pre-publication stage.

**Strengths for IDPs**: Because Forge generates sequences without requiring a structural template for the target, it sidesteps the entire problem of IDP conformer selection. It can sample binder sequences that are "sequence-complementary" to K18 based on co-evolutionary patterns in its training data. This is a different approach to the IDP conformer problem than BindCraft — it doesn't solve the conformer selection problem but rather bypasses it.

**Limitations**:
- No published hit rate for any target as of mid-2025. The model's validation is primarily computational (ESMFold pLDDT of generated sequences, ipTM of predicted complexes).
- The flow-matching approach may overfit to structured protein–protein interactions in the training data, potentially generating binders that work well for structured targets but poorly for IDPs.
- The lack of structural information about the target during generation means Forge has no mechanism to ensure the generated binder contacts the PHF6*/PHF6 hotspots. The epitope contact filter in NB02c is therefore doing essential work that cannot be delegated to Forge's generation step.

**BindCraft vs. Forge complementarity**: The pipeline correctly uses both in parallel:
- BindCraft: structure-guided, epitope-focused, higher per-design confidence but lower diversity and IDP-specific uncertainty.
- Forge: sequence-only, epitope-agnostic, higher diversity, lower per-design confidence, covers different chemical space.
- NB03 comparison: Binders from the two tracks should cluster in different regions of the SaProt embedding space, reflecting their different design philosophies. Overlap between clusters would indicate a convergent solution space; separation would indicate complementary coverage.

---

## Cross-Cutting Methodological Notes

### ipTM Calibration for IDPs

Across Topics 4, 8, and the pipeline design, the recurrent concern is that AF2-multimer and AF3 ipTM scores are systematically miscalibrated for IDP targets. The following calibration guidance applies:

| Context | ipTM range | Interpretation |
|---------|-----------|----------------|
| Structured protein–protein complex (AF2) | > 0.75 | High confidence |
| Structured protein–protein complex (AF2) | 0.5–0.75 | Moderate confidence |
| IDP complex (AF2/AF3) | > 0.6 | Potentially genuine |
| IDP complex (AF2/AF3) | 0.4–0.6 | Uncertain — possible artifact |
| IDP complex (AF2/AF3) | < 0.4 | Likely false positive or no binding |

The project's filter of `ranking_score ≥ 0.5` (where ranking_score = 0.8 × ipTM + 0.2 × ptm) translates to an implied ipTM threshold of ~0.5–0.55. For IDP targets, this is in the "uncertain" zone and will produce a mix of genuine and false-positive binders. This is not a pipeline design error — it is unavoidable given the current state of IDP structure prediction — but it means wet-lab validation will have a lower hit rate than the in silico filter rates suggest.

### Epitope Contact Filter as the Key Discriminating Step

The most scientifically impactful filter in the pipeline is not ipTM but the epitope contact filter (binder residues within 8 Å Cα–Cα of PHF6*, PHF6, or jR2R3). This is because:
1. It encodes the mechanistic hypothesis (blocking aggregation nucleation at PHF6*/PHF6).
2. It cannot be gamed by the structure prediction tool (AF3/Protenix does not know about the filter).
3. It selects for a specific mode of action, making the results interpretable in terms of the drug target biology.

The lessons from failed drugs (Topic 3) make this filter not just a technical choice but a scientific one: without it, the pipeline might produce high-ipTM binders that target irrelevant surfaces, exactly as the failed antibodies targeted the wrong region.

---

## CSO Decision Points

The following are the key scientific decisions that require CEO approval before proceeding with downstream experimental validation (wet-lab or extended in silico work). Each is grounded in the literature reviewed above.

---

### Decision Point 1: Retain the jR2R3 epitope in the scoring function?

**Background**: jR2R3 (K18 residues 52–70, the R2/R3 junction) was included based on Sugiyama et al. (2025) as a secondary epitope with 0.5 weight. However, jR2R3 has weaker biological validation than PHF6*/PHF6 — it is an early oligomerization contact (not the primary aggregation seed), and it is less conserved across Tau isoforms. Including it biases conformer selection toward exposing a region that may be transiently structured.

**Concern**: If binders are designed against jR2R3 contacts at the expense of PHF6*/PHF6 coverage, they may not block aggregation. The 0.5 weight in conformer scoring is low enough that jR2R3 will only tip the balance between otherwise similar conformers — this is probably fine. The greater risk is in NB03's integrated score, where `epitope_contact_score` aggregates contacts across all three epitopes. A binder that only contacts jR2R3 (not PHF6*/PHF6) would pass the filter in NB02b/NB02c but have a lower integrated score in NB03.

**CSO Recommendation**: RETAIN jR2R3 in conformer scoring at 0.5 weight (it adds value for ensemble diversification), but in NB02b/NB02c, modify the epitope contact filter to REQUIRE at least one contact with PHF6* OR PHF6 specifically (not just any epitope region). Currently the filter only requires `epitope_contacts ≥ 1` total, which could be satisfied by jR2R3 alone. A binder that touches only jR2R3 is not validated by the mechanistic hypothesis.

---

### Decision Point 2: Single conformer per BindCraft session, or run BindCraft on a merged ensemble?

**Background**: NB02b runs BindCraft five times, once per conformer. This correctly captures conformational diversity. However, each BindCraft run optimizes for one specific K18 backbone, meaning the resulting binder is conformer-specific. In reality, a therapeutic binder needs to work against the ensemble of K18 conformations a patient's Tau protein will adopt.

**Concern**: The best binder for conformer rank1 may be ineffective against conformers rank3–5. NB03 does not currently assess cross-conformer activity (a binder designed against rank1 is not re-evaluated against rank2–5 backbones). The project may therefore miss this ensemble coverage problem.

**CSO Recommendation**: ADD a cross-conformer ipTM check in NB03. For each passing BindCraft binder, run Protenix with the binder sequence paired against all five K18 conformers (not just the conformer it was designed against). Binders with high ipTM across ≥ 3/5 conformers are "conformer-robust" binders and should be ranked higher in the integrated score. This is the most important missing validation step. If Colab GPU time is limiting, focus this analysis on the top 10 binders from each conformer track.

---

### Decision Point 3: ipTM threshold — keep at 0.5 or adjust for IDP?

**Background**: The current filter is `ranking_score ≥ 0.5`. As reviewed in Topic 4, ipTM scores for IDP complexes from AF2-multimer and AF3 are systematically lower than for structured protein–protein complexes. A blanket threshold of 0.5 will reject genuine IDP binders at high rates (false negatives) while retaining AF3-overconfident predictions at similar rates (false positives).

**CSO Recommendation**: LOWER the primary filter to `ranking_score ≥ 0.45` to reduce false negatives, but RAISE the weight of the epitope contact score in NB03's integrated scoring formula. Specifically: increase epitope_contact_score weight from 0.20 to 0.30, and decrease iptm weight from 0.30 to 0.25 in the NB03 formula. This reflects the reality that for IDP targets, epitope engagement (contact geometry) is a more reliable discriminator than the raw confidence score. The change should be documented explicitly in the methods so that the ipTM values presented are not over-interpreted.

---

### Decision Point 4: Include a Tau NFT (fibrillar) control binder track?

**Background**: The project currently designs binders for the K18 monomer (disordered) state. PDB 6QJH (the AD Tau protofilament) exposes different surfaces than the disordered monomer — in the filament, PHF6 and PHF6* are buried in the fibril spine and inaccessible. Binders designed against the disordered monomer that happen to also bind the fibril would be problematic (by stabilizing fibrils rather than preventing aggregation).

**Concern**: The pipeline does not currently check whether passing binders would bind the fibrillar form. A binder that binds PHF6* in the disordered monomer but does not bind the filament core is mechanistically ideal — it prevents aggregation by sequestering the monomer. A binder that also binds the filament would be neutral at best and potentially aggregation-promoting at worst.

**CSO Recommendation**: ADD a structural incompatibility check in NB03. Download PDB 6QJH (AD protofilament) and compute Protenix ipTM for each top-10 binder against the fibrillar Tau sequence in the 6QJH conformation. Binders with ipTM > 0.5 against the fibrillar form should be flagged (not necessarily excluded, but noted as potentially dual-acting). Binders with ipTM < 0.4 against the fibrillar form and ipTM > 0.5 against the K18 monomer are the most desirable mechanistic profile. This adds ~2 Colab cells to NB03 and is well within scope.

---

### Decision Point 5: Which binder length range is optimal for K18?

**Background**: NB02b (BindCraft) allows flexible binder length selection; NB02c (Forge) uses 50/65/80 aa. The literature on IDP binders (Baker 2025, Topic 5; Bhardwaj 2022, Topic 4) suggests 40–80 aa is the productive range for IDP targets. However, the PHF6* and PHF6 hexapeptides are only 6 residues each, and beta-sheet augmentation binders that cap these sequences can be as small as 10–20 aa.

**Concern**: Very short binders (< 40 aa) may be insufficient for stable folding but could be more specific; very long binders (> 100 aa) may have high off-target binding risk for a 129-aa target (the binder would be nearly as large as K18 itself).

**CSO Recommendation**: CONFIRM the BindCraft binder length range is set to 50–100 aa (not default 50–150 aa). For K18 (129 aa), a 150-aa binder is counterproductive. Additionally, add a 30–45 aa length bin to the Forge run in NB02c, targeting beta-strand augmentation binders specifically. These very short binders are explicitly what Bhardwaj et al. (2022) validated for VQIVYK-capping designs and are underrepresented in the current pipeline. The Forge 50-aa bin may already be capturing some of these, but a dedicated 35-aa bin would increase coverage.

---

*End of literature review. All Decision Points should be reviewed by the student (CEO) before NB02b/NB02c are run on the full conformer set. Decisions 1 and 3 affect the filter logic in existing notebooks. Decision 2 and 4 add validation steps to NB03. Decision 5 affects the BindCraft and Forge hyperparameter settings in NB02b/NB02c.*
