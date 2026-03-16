# Strategic Synthesis & Final Recommendations
**CHEM 280 Final Project — Tau K18 Binder Design Pipeline**
**Synthesized**: Sunday March 15, 2026, 11 PM
**Deadline**: Wednesday March 18, 2026 EOD
**Time remaining**: ~59 hours total; ~47 hours to data-complete

---

## Section 1: Cross-Report Conflict and Reinforcement Analysis

### Reinforcing Signals (All Three Reports Agree)

**The epitope contact filter is the most important scientific discriminator.**
The literature review, consultant risk assessment, and data engineer review all converge on this. The literature review identifies it explicitly: "the most scientifically impactful filter in the pipeline is not ipTM but the epitope contact filter." The data engineer flags two bugs (BUG-20, BUG-23) that could silently zero out this filter entirely. The consultant lists zero passing binders as Risk 3. These are the same problem seen from three angles — if the epitope contact computation is wrong, the entire scientific thesis of the pipeline collapses. This is the highest-priority correctness concern in the system.

**NB02b is the chokepoint and the highest-probability failure point.**
The consultant identifies this as the single bottleneck. The data engineer independently identifies four separate bugs in NB02b (BUG-17 through BUG-20) that all produce silent failures — BindCraft crash ignored, wrong binder chain extracted from ProteinMPNN, ranking_score silently zero, epitope contacts silently zero. Every one of these bugs routes to the same outcome: an empty CSV sent to NB03 with no warning. The consultant says one crash costs 3–5 hours; the data engineer says these bugs mean you may not know a run failed until NB03 reads an empty file. This is a compounding problem.

**ipTM thresholds are miscalibrated for IDPs and the literature justifies relaxing them.**
The consultant recommends lowering the threshold if pass rates are low (to 0.4). The literature review independently derives the same number from first principles (Johansson-Åkhe 2022 benchmark on AF2 IDP complexes) and frames it as methodologically justified, not a shortcut. The data engineer's BUG-19 (NB02b always returns ranking_score=0.0 because the key doesn't exist) means the 0.5 threshold currently rejects every binder from NB02b regardless. Fix the bug first; adjust the threshold second only if genuine pass rates remain low.

**lDDT from NB01 is likely None and this is acceptable.**
The consultant flags lDDT parsing as Risk 2 (moderate probability). The data engineer flags the regex as BUG-11/12. The context statement confirms NB01 has run and "lDDT is None due to mode 1, but the run was accepted — top 5 conformers selected by 3Di epitope conservation." All three reports reach the same conclusion: the fallback to 3Di-only ranking is scientifically defensible and the pipeline should proceed without retrying NB01.

### Conflicts and Resolutions

**Conflict: Scope of NB03 enhancements.**
The literature review proposes two new NB03 analyses as CSO Decision Points: (1) cross-conformer ipTM validation (Decision Point 2) and (2) fibrillar Tau incompatibility check using PDB 6QJH (Decision Point 4). The consultant's schedule shows NB03 running Tuesday morning with a 30–60 minute runtime. Adding two Protenix loops to NB03 would extend runtime substantially — potentially 3–6+ hours depending on binder count — and would require data from all five NB02b runs plus internet access to fetch PDB 6QJH on Colab.

**Resolution**: Both enhancements are scientifically compelling but out of scope for Wednesday. See Section 5 (CSO Decision Points). The correct move is to note both as explicitly planned future work in the report, not to implement them now.

**Conflict: How many NB02b sessions to run (5 vs. 3).**
The consultant's primary schedule runs all 5 conformers but explicitly identifies "cut to 3 conformers" as the first and least-impactful scope reduction. The literature review says 40–80 aa binders covering diverse conformers is correct. There is no scientific conflict here — 3 conformers is adequate; 5 is better if time allows.

**Resolution**: Decide Monday morning based on the score spread across the top 5 conformers. If ranks 4 and 5 score < 0.7 × rank 1's combined score, drop them. This is a time management decision, not a scientific one.

---

## Section 2: Final Call — Specific, Actionable, Prioritized

### Immediate Actions (Before Sleeping Tonight)

1. **Verify both Colab sessions (NB01, NB02c) are still running.** Open both tabs. Confirm cells are executing, not disconnected. If either has disconnected, restart it now.

2. **Confirm Google Drive is mounted in both sessions** and that output files are appearing at the expected Drive paths. This is the only protection against losing all overnight compute. If Drive is not mounted in either session, this is a fire — fix it before sleeping.

3. **Do not start NB02b tonight.** No monitoring capability overnight = high probability of losing a session without knowing it.

### Monday Morning (8 AM — First 30 Minutes)

4. **Download and verify NB01 outputs.** Confirm exactly 5 PDB files in `binder_ready/`. Check that conformers have plausible scores. The lDDT=None situation is already accepted; proceed.

5. **Download and verify NB02c outputs.** Check that `forge_binder_filter_passing.csv` has rows. If it has zero rows, lower `MIN_RANKING_SCORE` to 0.45 in NB02c config and inspect the epitope contact function — BUG-23 may have silently zeroed all contacts.

6. **Before starting NB02b run 1: apply the three code fixes in Section 3.** BUG-19 is already fixed per the context. Verify BUG-17 (BindCraft returncode) and BUG-18 (ProteinMPNN chain order) are also addressed. These two bugs can silently produce empty CSVs. Ten minutes of verification now saves a 4-hour re-run.

7. **Start NB02b run 1 (rank01 PDB).** Keep the browser tab visible. Do not start the next run until this one finishes and you have downloaded the CSV.

### Monday Decision Gate

8. **After NB02b run 1 completes:** Look at how many binders passed. If zero passed: check BUG-20 (epitope residue numbering) and BUG-18 (chain order). Do not start run 2 until you understand why run 1 produced zero results. A known-bad configuration run 5 times is worse than diagnosing it once.

9. **Decide 3 vs. 5 conformers.** Look at the conformer ranking CSV from NB01. If ranks 4 and 5 score significantly below ranks 1–3, stop at 3 conformers. This decision should be made Monday by noon, not Tuesday night.

### Tuesday (Target: All CSVs in Hand by 2 AM, NB03 Running by 9 AM)

10. **Start NB03 as soon as all CSVs are collected.** Do not wait for perfect results. If some BindCraft runs produced zero passing binders, proceed with what you have. A real null result for some conformers is a result.

11. **Begin writing Methods Monday night.** You know the pipeline well enough to write it before results are complete. Methods + Introduction can be done before Tuesday morning.

---

## Section 3: High-Priority Bug Fixes (Data Engineer's Remaining HIGH Bugs)

The following remaining HIGH bugs are worth fixing before the deadline, ordered by impact. BUG-19 is already fixed; BUG-11, BUG-28, BUG-30 are already fixed. The bugs below are the next tier.

### Fix 1: BUG-17 — BindCraft Subprocess Failure Is Silent (NB02b cell `cell-11`)

**Why it matters**: If BindCraft crashes (PyRosetta license, AF2 weights, OOM), execution silently continues to a glob that returns empty, then raises a misleading `FileNotFoundError` pointing at the CSV rather than BindCraft. You diagnose the wrong thing at 1 AM.

**Current code:**
```python
subprocess.run(['python', '/content/BindCraft/bindcraft.py', '--settings', settings_path])
```

**Fix:**
```python
result = subprocess.run(
    ['python', '/content/BindCraft/bindcraft.py', '--settings', settings_path],
    capture_output=True, text=True
)
if result.returncode != 0:
    print("=== BindCraft STDERR ===")
    print(result.stderr[-3000:])  # last 3000 chars to avoid Colab cell overflow
    raise RuntimeError(f"BindCraft exited with code {result.returncode}. See stderr above.")
```

---

### Fix 2: BUG-18 — ProteinMPNN Chain Order Assumed (NB02b cell `cell-13`)

**Why it matters**: If ProteinMPNN outputs chain B (binder) before chain A (target), `parts[1]` extracts the K18 sequence (129 aa) as the "binder" and submits it to Protenix. Protenix predicts a K18:K18 homodimer. The ipTM will be meaningless and all binders will fail the filter — or worse, some will pass with the wrong sequence. This is a science-destroying silent error.

**Current code:**
```python
parts = lines[i+1].strip().split('/')
binder_seq = parts[1] if len(parts) > 1 else parts[0]
```

**Fix:**
```python
parts = lines[i+1].strip().split('/')
# Validate: binder must be in the expected length range, not the full K18 length
binder_candidates = [p for p in parts if BINDER_LEN_MIN <= len(p) <= BINDER_LEN_MAX]
if len(binder_candidates) == 0:
    print(f"  WARNING: No part of MPNN output is in binder length range "
          f"[{BINDER_LEN_MIN}, {BINDER_LEN_MAX}]. Parts: {[len(p) for p in parts]}. Skipping.")
    continue
binder_seq = binder_candidates[0]
```

This relies on `BINDER_LEN_MIN` and `BINDER_LEN_MAX` being defined in the config cell (they are, from the design settings). If somehow both sequences are in range, it takes the first — which is still better than blindly taking index 1.

---

### Fix 3: BUG-32 — NB03 Relative Paths Don't Resolve to Drive (NB03 cell `cell-7`)

**Why it matters**: NB03 never `%cd`s to Drive. Its relative paths (`results/esmfold/forge`, `tau_k18_binder_ready/rank01_filtered_0001.pdb`) resolve to `/content/` which contains nothing. NB03's folding loop will find no existing PDBs and will attempt to re-fold everything from scratch — writing to `/content/` (ephemeral, not Drive-backed) and silently replacing the NB01/NB02b outputs. This is not a minor issue; it makes NB03 run much slower and destroys result traceability.

**Current code (cell `cell-7`):**
```python
FORGE_PDB_DIR     = 'results/esmfold/forge'
TARGET_PDB        = 'tau_k18_binder_ready/rank01_filtered_0001.pdb'
BINDCRAFT_PDB_DIR = 'results/esmfold/bindcraft'
```

**Fix:** Replace with absolute Drive paths, consistent with NB01/NB02b/NB02c path conventions. Add this immediately after the Drive mount cell in NB03:

```python
import os

# -- Set CWD to Drive output root so relative paths resolve correctly --
_NB03_DRIVE_ROOT = '/content/drive/MyDrive/CHEM_280/results'
assert os.path.isdir(_NB03_DRIVE_ROOT), f"Drive not mounted or path wrong: {_NB03_DRIVE_ROOT}"
os.chdir(_NB03_DRIVE_ROOT)
print(f"CWD set to: {os.getcwd()}")

# Now relative paths will resolve against Drive
FORGE_PDB_DIR     = 'nb02c_output/esmfold'
TARGET_PDB        = 'nb01_output/binder_ready/rank01_filtered_0001.pdb'
BINDCRAFT_PDB_DIR = 'nb02b_output/esmfold'
```

**Note**: Verify the actual Drive subdirectory names match what NB01/NB02b/NB02c actually wrote. The path stems above should match the `NB01_OUT`, `NB02B_OUT`, `NB02C_OUT` variables in those notebooks. If they differ, update accordingly.

---

### Remaining HIGH Bugs: Worth Noting but Not Fixing Before Wednesday

**BUG-29** (FoldMason JSON schema assumed): NB01 has already run and produced results. The JSON parsing will either work or it won't — there is no time to make it version-robust before NB03 runs. Add a `try/except KeyError` guard if NB03 fails on the JSON step, printing `fm_data.keys()` to see what FoldMason actually produced.

**BUG-31** (known_proximity NaN fill): Scientifically valid concern, but the impact is that some binders score neutral instead of low on one metric. Will not produce incorrect rankings at the top, only at the bottom. Not worth fixing under deadline pressure.

**BUG-13** (PDB lookup glob ambiguity): Only matters if `N_CONFORMERS >= 1000`. With 400 conformers, the `filtered_0001` through `filtered_0400` range has no ambiguous substring matches. Safe to leave.

---

## Section 4: Follow-Up Questions for the Consultant

The consultant's timeline estimate assumes specific runtime parameters that need confirmation before Monday morning planning is finalized.

**Question 1: What is the actual BindCraft runtime for a K18 conformer on a Colab T4, given 50 hotspot residues and a 50–100 aa binder length range?**

The consultant's schedule assumes "3–5h per conformer (BindCraft dominates)." This estimate drives the entire Monday–Tuesday schedule. The range is wide — 3h vs. 5h per session is the difference between finishing NB02b Tuesday at 2 AM vs. Tuesday at 10 AM or later. Has this been timed on the actual NB02b settings (specifically: how many hallucination trajectories are configured, and what is `n_designs`/`n_mpnn_seqs` set to)? The consultant needs to look at the NB02b config cell and give a tighter estimate before the student commits to the serial-5-session plan.

**Question 2: Is the Google Drive continuous write integration in NB02b actually tested and working, or is it assumed?**

The consultant says "verify Drive is mounted so intermediate CSV outputs are written to Drive continuously." Drive mounting is confirmed in the setup cell. But Drive-continuous writing means every BindCraft output CSV is written to Drive in real-time, not buffered to `/content/` and copied at the end. If NB02b only writes to Drive at the final cell, a crash at hour 4 of a 5-hour run loses everything. Has the Drive write pattern in NB02b been verified to be continuous? This changes the risk calculus for overnight runs significantly.

**Question 3: Given that NB02b has four separate silent-failure bugs (BUG-17, -18, -19 [fixed], -20), what is the actual probability that a first run of NB02b produces zero passing binders, and does that change the recommendation to start Monday 9 AM vs. doing a 30-minute diagnostic test run first?**

The consultant recommends starting NB02b run 1 (rank01) Monday at 9 AM without qualification. But if the likelihood of a silent-fail run is moderate (due to BUG-18 or BUG-20 being unpatched), the rational move may be to run a 10-design test batch first (set `n_designs=10` in the config, run to completion, verify the CSV has rows and the binder sequences are the right length) before committing a full 3–5 hour session. The consultant should give an explicit recommendation on test-run vs. full-run Monday morning based on whether the fixes in Section 3 are applied.

---

## Section 5: CSO Decision Points — Final Verdicts

### Decision Point 1: Retain jR2R3 epitope in scoring function?

**Verdict: (a) IMPLEMENT NOW — retain as-is, with one targeted fix.**

The jR2R3 weight of 0.5 in conformer scoring is correct and already baked into NB01's completed run. Do not change the conformer scores retroactively. However, the CSO's specific concern — that a binder contacting only jR2R3 (not PHF6*/PHF6) passes the filter — is valid and easy to address. In NB02b and NB02c, modify the epitope contact pass condition from `epitope_contacts >= 1` to `phf6_contacts + phf6star_contacts >= 1`. The contact dictionaries are already computed per-epitope in both notebooks; this is a one-line logical change in the filter condition. This is the mechanistically correct filter and takes 5 minutes to implement before Monday's NB02b run 1.

### Decision Point 2: Cross-conformer ipTM validation in NB03?

**Verdict: (c) OUT OF SCOPE for Wednesday.**

Running Protenix for each top binder against all 5 K18 conformers multiplies Protenix compute by 5×. If NB03 receives ~50 passing binders total (10 per conformer), cross-conformer validation means 250 Protenix jobs. At ~3–5 minutes per job, that is 12–20 GPU hours — more compute than the entire rest of the pipeline. Even if prioritized to the top 10 binders only, this adds 30–50 minutes of Protenix time plus substantial code additions to NB03. The science is excellent but the implementation is not a Tuesday-morning task.

**Recommended disposition**: Add one paragraph to the Discussion section noting this as the primary limitation and explicit next step.

### Decision Point 3: Adjust ipTM threshold for IDP calibration?

**Verdict: (b) IMPLEMENT IF TIME ALLOWS — but fix BUG-19 first and see what you get.**

The current NB02b ranking_score=0.0 bug (BUG-19, already fixed) means all NB02b binders have been rejected at the 0.5 threshold. With the fix applied, actual pass rates from Monday's NB02b runs will tell you whether 0.5 is too aggressive. Do not pre-emptively lower the threshold before seeing real data. If Monday's run 1 produces zero passing binders after BUG-19 is fixed, then lower to 0.45 and add a methods note citing Johansson-Åkhe 2022.

The weight rebalancing the literature review recommends (epitope weight 0.20 → 0.30 in NB03 integrated score) is low-risk and scientifically grounded. This can be done in NB03's config cell in 2 minutes. **Implement this change when setting up NB03 on Tuesday morning.**

### Decision Point 4: Add fibrillar Tau incompatibility check (PDB 6QJH)?

**Verdict: (c) OUT OF SCOPE for Wednesday.**

This is a compelling validation that would strengthen the paper significantly. It requires: fetching PDB 6QJH from the RCSB on Colab, extracting the filament chain, running Protenix for each top binder against the filament, and adding analysis cells to NB03. The compute is manageable (~10–20 Protenix jobs for the top binders) but the coding work adds 2–3 notebook cells and new analysis logic to NB03. More importantly, NB03 will run Tuesday morning with potentially no margin for debugging new cells.

**Recommended disposition**: Frame this as the most important future experimental validation direction in the Discussion. The mechanistic argument (monomer-selective vs. fibril-binding) is already clearly supported by the literature review and the Gosuranemab/Semorinemab failure analysis.

### Decision Point 5: Adjust binder length range for K18?

**Verdict: (c) OUT OF SCOPE for Wednesday — but verify BindCraft settings before run 1.**

Adding a 30–45 aa Forge length bin in NB02c is not actionable; NB02c is already running. BindCraft's length range is set in the NB02b config cell. The CSO recommendation to cap at 50–100 aa (not 50–150 aa default) is worth checking Monday morning before run 1, but only if the NB02b config cell shows a length upper bound > 100 aa. A quick look at the `binder_length` parameter and a one-number change if needed. Do not let this become a 30-minute investigation — check the number, fix it if wrong, start run 1.

---

## Section 6: The Single Most Important Thing in the Next 2 Hours

**Verify Google Drive is writing continuously in both active Colab sessions.**

Here is why this is the one thing above all others:

Everything downstream — NB02b Monday morning, NB03 Tuesday morning, the report, the submission — depends on NB01 and NB02c outputs surviving overnight. NB01 has already produced valid results (top 5 conformers selected). NB02c is running on an A100. If either session disconnects tonight and Drive is not mounted, those results are gone and the only way to recover is a full re-run tomorrow, consuming 4–6 hours that belong to NB02b.

Checking Drive takes 2 minutes per session. In each session, run:
```python
import os
print(os.listdir('/content/drive/MyDrive/CHEM_280/results/'))
```

If this lists files that have been written, Drive is mounted and writing. If it errors or shows an empty directory, mount Drive immediately via the Drive cells at the top of the notebook.

The bugs, the science decisions, the timeline planning — all of that matters but all of it can be addressed Monday morning. Losing tonight's compute is the one thing that cannot be recovered from in time.

After verifying Drive: set an alarm for 7:30 AM Monday and go to sleep.

---

## Appendix: Known-Fixed Bugs Summary (Context)

The following 6 critical bugs are already fixed and do not require action:

| Bug | Location | Fix Applied |
|-----|----------|-------------|
| BUG-24 | NB02c cell-17 | `np.linalg.norm()` replacing `__abs__()` |
| BUG-28 | NB03 cell-5 | `import subprocess` added |
| BUG-19 | NB02b cell-15 | `ranking_score` computed as `0.8*iptm + 0.2*ptm` |
| BUG-11 | NB01 cell-02d9bbe9 | `--report-mode 2` (JSON) now used |
| BUG-30 | NB03 cell-20 | HDBSCAN guard added for zero-cluster case |
| BUG-02/10 | NB01 cell-f0b4cf53 | Empty `keep_idx` fallback added |

NB01 has run and produced valid top-5 conformers. NB02c is currently running on Colab A100. Do not rerun either.

---

*Synthesis complete. Three reports cross-referenced; no fundamental conflicts in scientific direction. Execution risk is concentrated in NB02b reliability and cross-notebook data flow. Pipeline science is sound.*
