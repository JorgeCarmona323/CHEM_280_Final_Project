# Deadline Risk Assessment — CHEM 280 Final Project
**Prepared**: Sunday March 15, 2026, 10:00 PM
**Deadline**: Wednesday March 18, 2026 (end of day)
**Hard cutoff for data generation**: Wednesday March 18 morning (~9 AM)
**Time available**: ~59 hours total; ~47 hours to data-complete

---

## Executive Summary

The pipeline is achievable before the deadline, but it requires zero unrecovered failures
in NB01 and NB02b. The schedule has almost no slack. The single largest threat is the
5 × serial NB02b sessions — each runs 3–5 hours and cannot overlap on a single Colab
account. Any one session that hangs, crashes, or produces zero passing binders requires a
re-run that will consume time that does not exist. Mitigation for that specific scenario
must be decided tonight, not Tuesday.

---

## 1. Realistic Hour-by-Hour Schedule

All times are **calendar time**, not compute time. Parallelism is limited by having one
active Colab A100 (NB02c) and one T4 (NB01) right now.

### Sunday March 15

| Time | NB01 (T4) | NB02c (A100) | Action Required |
|------|-----------|--------------|-----------------|
| 10 PM | RUNNING | RUNNING | Monitor. Check GPU dashboards for disconnection. Do not close browser tabs. |

### Monday March 16

| Time | NB01 (T4) | NB02c (A100) | Action Required |
|------|-----------|--------------|-----------------|
| ~12–2 AM | DONE | RUNNING | NB01 likely finishes overnight (400 conformers, Starling + FoldMason on T4 ≈ 2–4h). Download `tau_k18_phase1_results.zip` immediately. Verify 5 PDBs are in `binder_ready/`. |
| ~3–4 AM | — | DONE | NB02c likely finishes (Forge pipeline estimated 2–2.5h). Download `forge_designs.zip`. |
| 8 AM | — | — | **Wake up. Verify both downloads exist and are non-empty.** If either failed: restart immediately. |
| 9 AM | NB02b run 1 | — | Start NB02b on rank01 PDB. Needs new T4 session (or A100 if available). |
| ~12–2 PM | NB02b run 1 DONE | — | Download `bindcraft_rank01_passing.csv`. Start NB02b run 2 (rank02). |
| ~3–5 PM | NB02b run 2 DONE | — | Download CSV. Start NB02b run 3 (rank03). |
| ~6–8 PM | NB02b run 3 DONE | — | Download CSV. Start NB02b run 4 (rank04). |
| ~9–11 PM | NB02b run 4 DONE | — | Download CSV. Start NB02b run 5 (rank05). |

### Tuesday March 17

| Time | NB02b | Action Required |
|------|-------|-----------------|
| ~12–2 AM | Run 5 DONE | Download final CSV. All 6 CSVs now in hand. |
| 8 AM | — | **Start NB03.** All inputs available. |
| ~9–11 AM | — | NB03 finishes (30–60 min estimated). Download `tau_binder_comparison_results.zip`. |
| Rest of day | — | **Write report and slides.** Data is complete. |

### Wednesday March 18

| Time | Action |
|------|--------|
| All day | Report writing, figure polish, presentation prep. No compute needed. |
| End of day | Submit. |

### Summary

| Stage | Estimated Duration | Completion Target |
|-------|--------------------|-------------------|
| NB01 | 2–4h | Mon ~2 AM |
| NB02c | 2–2.5h | Mon ~4 AM |
| NB02b ×5 (serial) | 3–5h each = 15–25h total | Mon 9AM → Tue 2 AM |
| NB03 | 0.5–1h | Tue ~11 AM |
| Report writing | ~1.5 days | Wed EOD |

The schedule is tight but executable if NB02b runs proceed without major failures.

---

## 2. Critical Path

**The single bottleneck is NB02b — five serial BindCraft sessions, each 3–5 hours,
each dependent on the previous session completing and yielding a downloaded CSV.**

Why this is the chokepoint:

- NB01 and NB02c are already running in parallel and will likely finish tonight.
- NB03 is fast (30–60 min) and depends only on having the CSVs.
- The only thing preventing NB03 from running Tuesday morning is NB02b finishing
  5 complete sessions in roughly 20 hours of calendar time starting Monday morning.
- BindCraft runtime on T4 is the most variable step in the pipeline — it runs
  AF2-multimer internally and design quality varies. The notebook itself states
  "~3–5h per conformer (BindCraft dominates)". On a shared Colab T4, actual
  runtime can be higher.
- Five sessions × worst-case 5h = 25 hours. Starting Monday 9 AM, that is
  Tuesday 10 AM. That leaves barely enough time for NB03 before the Tuesday
  "done with data" target — and no recovery buffer if even one session fails.

**If the BindCraft runs average 5h each rather than 3h each, the data will not
be complete until Tuesday afternoon, leaving only Wednesday for writing, which
is survivable but stressful.**

---

## 3. Top 3 Risks

### Risk 1 — NB02b runtime overrun or session crash
**Probability: HIGH**
Colab sessions disconnect on inactivity, after 12h, or when GPU quota is exhausted.
BindCraft runs AF2-multimer and is memory-intensive. If a session crashes mid-run,
all progress for that conformer is lost and the session must restart from the beginning.
Five serial sessions with no overlap means one crash costs 3–5 hours of real time.

### Risk 2 — NB01 lDDT parsing failure on Colab
**Probability: MODERATE**
The session handoff document explicitly flags that lDDT parsing from FoldMason HTML
was not confirmed in local testing: *"The regex pattern may not match FoldMason's
current HTML format."* If lDDT is unavailable, the scoring falls back to 3Di epitope
consistency only. This is non-blocking for the pipeline, but it means the top-5
conformer selection is less scientifically differentiated. More importantly, if this
causes an unhandled exception rather than a graceful fallback, the notebook could
fail and require manual debugging at 2 AM.

### Risk 3 — Zero or very few passing binders from NB02b/NB02c
**Probability: LOW TO MODERATE**
The filter requires `ranking_score ≥ 0.5 AND epitope_contacts ≥ 1`. For IDPs like
Tau K18, Protenix may consistently predict weak binding (ranking_score < 0.5) because
the target has no stable fold. If multiple conformers produce zero passing binders,
the NB03 integrated analysis has little to work with, weakening the scientific
narrative even if everything runs on time.

---

## 4. Mitigation Strategies

### Mitigation for Risk 1 — BindCraft session crashes

**Before starting Monday**: Enable Google Drive integration in NB02b. The notebook
already has Drive mount code (`drive-setup` cell). Verify Drive is mounted so that
intermediate CSV outputs are written to Drive continuously, not just at the end.
This means a mid-run crash does not lose the BindCraft designs already generated.

**Use Colab Pro or Pro+ if you have it.** Longer session timeouts and background
execution are the only real defenses against disconnection during overnight runs.

**Start sessions during waking hours only.** Do not kick off a new NB02b session
right before going to sleep unless you have Drive integration confirmed working.
A crashed session you don't notice for 8 hours is a 8-hour setback.

**Keep a browser tab open and disable computer sleep** during each session. Colab
disconnects idle tabs.

**Fallback option**: If two NB02b sessions fail to produce passing binders and
you are running short on time, stop at 3 conformers (rank01–03). Three BindCraft
CSVs plus one Forge CSV is sufficient for NB03 to run and for a credible comparison.
The scientific story changes slightly (sampling 3 conformers rather than 5) but
the methodology is intact.

### Mitigation for Risk 2 — lDDT parsing failure

**Before running NB01 on Colab**: Check whether lDDT shows up after the first
10–20 conformers by inspecting `tau_k18_foldmason_msa.html` in the Drive output
directory. If the `"lddt"` key is absent from the HTML, update the regex per the
instructions in SESSION_HANDOFF.md cell-20. This is a 5-minute fix but must be
done proactively, not at 3 AM when NB01 finishes.

Since NB01 is already running right now, you cannot change it mid-run. After it
completes, inspect the output immediately. If lDDT is missing, verify that the
`binder_ready/` folder has 5 PDBs anyway — the fallback ranking is scientifically
defensible for a course project.

### Mitigation for Risk 3 — Low pass rate from binder design

**If NB02c produces zero passing binders**: Relax the threshold to
`ranking_score ≥ 0.4` in the config cell. The threshold is arbitrary; 0.4 is
reasonable for an IDP target where predictions are inherently noisier.

**If multiple NB02b runs produce zero passing binders**: Same threshold relaxation.
Also consider that for a course project, presenting the pipeline with honest
reporting of low pass rates is scientifically stronger than inflating thresholds
to generate fake positives. If the filter is the problem, say so in the report.

**Prepare a fallback narrative**: Even if binder design yields few or no high-
confidence candidates, the NB01 conformer analysis (AFRC + FoldMason + epitope 3Di
scoring) is publishable on its own as a pipeline contribution. The comparison in
NB03 can still be framed as a null result or negative control if needed.

---

## 5. Scope Reduction Options

Ordered from least to most impact on the science:

### Cut 1 — Run only 3 conformers in NB02b (not 5)
**Time saved**: 2 × 3–5h = 6–10h
**Science impact**: Minimal. Three diverse conformers adequately sample the IDP
ensemble. The paper's narrative changes from "5 conformers" to "top 3 conformers,"
which is still defensible. NB03 still gets a meaningful comparison.
**Trigger**: If NB01 finishes and the top-5 conformers look highly similar in the
UMAP/ranking plot, this is even easier to justify.

### Cut 2 — Skip ProteinMPNN redesign in NB02b, submit BindCraft sequences directly to Protenix
**Time saved**: ~30–60 min per conformer (ProteinMPNN step)
**Science impact**: Moderate. ProteinMPNN redesign is scientifically important —
it removes the artifact of BindCraft's self-consistent optimization. However, for
a course project, you can acknowledge this limitation and present results as
"BindCraft designs validated by Protenix without MPNN redesign."
**Trigger**: Use this only if ProteinMPNN is failing to install or runs are timing out.

### Cut 3 — Skip NB02b entirely; use NB02c (Forge) + known binders only for NB03
**Time saved**: ~15–25h (all 5 BindCraft sessions)
**Science impact**: HIGH. Losing the structure-based track removes the key
method comparison. NB03 becomes a Forge-only analysis benchmarked against
published binders. This is still a real result, but the comparative angle
(sequence-based vs structure-based binder design for IDPs) is the most scientifically
interesting aspect of the project.
**Trigger**: Use only if it is Tuesday afternoon and NB02b is less than 50% done.

### Cut 4 — Simplify NB03 (skip SaProt, use ESM2; skip HDBSCAN, use k-means)
**Time saved**: ~20–30 min
**Science impact**: Low. ESM2 is a legitimate fallback embedding. k-means is
less principled than HDBSCAN but produces interpretable clusters. The notebook
already has fallback code for ESM2. Not worth doing unless NB03 fails to install
SaProt.

---

## 6. Action Items — Priority Ordered

### Tonight (Sunday, before sleeping)

1. **Verify both Colab sessions are still running.** Open each tab. Check that
   cells are executing, not disconnected.

2. **Confirm Google Drive is mounted in both sessions** and that output files are
   appearing in `MyDrive/CHEM_280/results/`. If Drive is not mounted, the session
   will write to `/content/` only and all output will be lost on disconnect.

3. **Do not start NB02b tonight.** It needs NB01 to finish first, and you cannot
   monitor a new session overnight without risk of losing it.

4. **Set an alarm for 7–8 AM Monday.** NB01 and NB02c will likely have finished.

### Monday morning (8 AM)

5. **Download and verify outputs immediately.** Check that `tau_k18_phase1_results.zip`
   contains exactly 5 PDB files in `binder_ready/`. Check that
   `forge_designs.zip` contains `forge_binder_filter_passing.csv` with at least
   some rows.

6. **Inspect the lDDT column** in `tau_k18_conformer_ranking.csv`. If all values
   are `None`, check the HTML and update the regex if needed — but do not block
   NB02b on this. Proceed with the 5 PDBs regardless.

7. **Start NB02b run 1 (rank01) immediately.** Do not do anything else first.
   Mount Drive. Run all cells. Keep the tab visible.

8. **Decide now whether to run 3 conformers or 5.** Look at the ranking scores
   of the top-5 conformers. If ranks 4 and 5 have significantly lower combined
   scores than rank 1–3, dropping them is justified and saves up to 10 hours.

### Monday through Tuesday

9. **Run NB02b sessions back-to-back**, downloading each CSV before starting the next.
   Do not queue multiple sessions — you can only reliably monitor one at a time.

10. **Save every CSV to Drive immediately** after it downloads. Name them clearly:
    `bindcraft_rank01_passing.csv`, etc. You will need all of them simultaneously
    for NB03.

### Tuesday morning

11. **Run NB03 as soon as all CSVs are collected.** Do not wait. NB03 takes 30–60
    minutes; the sooner it runs, the more time you have for writing.

12. **Generate all figures from NB03** (UMAP, epitope heatmap, ranking plot).
    These are the core figures for the report. Download the zip immediately.

### Tuesday afternoon and Wednesday

13. **Write the report first, then polish.** Structure: Introduction (IDP binder
    challenge) → Methods (pipeline) → Results (conformer ensemble, binder design,
    comparison) → Discussion (what worked, limitations, future directions).

14. **Do not wait for perfect results before writing.** Write the methods section
    Monday night using what you already know. Fill in results as data arrives.

---

## Bottom Line

You will make the deadline if:
- NB01 and NB02c finish tonight without crashes (likely — they are already running)
- NB02b runs proceed at 3–4h average per conformer with no re-runs needed
- You start NB02b runs Monday 9 AM and work through them continuously

You will miss the data deadline if:
- Two or more NB02b sessions crash without Drive backup
- BindCraft runtime averages 5h+ per conformer (pushing completion to Tue afternoon)
- You sleep through the NB01/NB02c completion and lose Tuesday morning

The most important thing you can do right now is verify Drive is mounted in both
running sessions. Everything else follows from the data surviving tonight.
