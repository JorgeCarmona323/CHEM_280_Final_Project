# Data Engineer Code Review — CHEM 280 Final Project
**Reviewer:** Senior Data Engineer (15 years scientific computing)
**Date:** 2026-03-15
**Scope:** Full pipeline — scripts, Colab notebooks, and environment files

---

## Summary Table

| # | File | Location | Severity | Category | Short Description |
|---|------|----------|----------|----------|-------------------|
| 1 | `run_nb01_local.py` | line 217–224 | **CRITICAL** | Wrong science | AFRC distance map last row/col always zero — biases deviation map |
| 2 | `run_nb01_local.py` | line 246 | **CRITICAL** | Crash | `np.array([ca_arrays[i] for i in keep_idx])` crashes if keep_idx is empty |
| 3 | `run_nb01_local.py` | line 430 | **HIGH** | Wrong output | `rank{rank_idx:02d}` uses DataFrame index (1-based), not rank order |
| 4 | `run_nb01_local.py` | line 176 | **MEDIUM** | Crash | `foldmason version` returncode not checked — silent failure |
| 5 | `run_nb01_local.py` | line 189 | **MEDIUM** | Security / correctness | `shell=True` with unsanitized `SEQUENCE` string |
| 6 | `run_foldmason_k18.sh` | line 90–108 | **HIGH** | Wrong output | JSON rename logic fragile — `msa_report.json` assumed, but FoldMason output name is version-dependent |
| 7 | `run_foldmason_k18.sh` | line 32 | **HIGH** | Crash | Path contains a literal space and parentheses; will break without quoting in callers |
| 8 | `run_foldmason_k18.sh` | line 52 | **MEDIUM** | Edge case | `find | wc -l` silently counts 0 if `FRAMES_DIR` has subdirectory `frame_*.pdb` files |
| 9 | `NB01` cell `f0b4cf53` | AFRC loop | **CRITICAL** | Wrong science | Same as issue #1 — AFRC distance map last residue is always zero (shared defect) |
| 10 | `NB01` cell `f0b4cf53` | keep_idx empty | **CRITICAL** | Crash | If Rg std=0 or all frames fail filter, `positions = np.array([...])` → IndexError |
| 11 | `NB01` cell `02d9bbe9` | FoldMason call | **HIGH** | Wrong science | `--report-mode 1` (HTML only) used in notebook; `run_foldmason_k18.sh` uses `--report-mode 2` (JSON) — these are inconsistent |
| 12 | `NB01` cell `694567f4` | `parse_lddt_from_html` | **HIGH** | Wrong science | Regex `[^}]*?` in `parse_lddt_from_html` is non-greedy but non-anchored; will silently return empty dict if FoldMason HTML format changes |
| 13 | `NB01` cell `c7f308cf` | PDB lookup | **HIGH** | Wrong output | FoldMason names sequences by full path stem, not `filtered_XXXX`; glob fallback `*{stem}*.pdb` may match wrong file |
| 14 | `NB01` cell `7c954c77` | `%cd` magic | **HIGH** | Colab-specific | `%cd {_nb01_out}` changes the notebook's CWD; all subsequent `os.path.abspath()` calls resolve relative to Drive path — safe only if `%cd` runs first every time |
| 15 | `NB01` cell `10b1bb53` | install | **MEDIUM** | Colab reliability | `pip install -q` suppresses errors; broken installs silently pass, failures only surface 3 cells later |
| 16 | `NB02b` cell `cell-7` | `EPITOPE_RANGES` | **CRITICAL** | Wrong science | NB02b redefines `EPITOPE_RANGES` with 4-tuple `(start, end, label, include_bool)` — incompatible with NB01/NB03 format `(start, end, label, weight_float)` |
| 17 | `NB02b` cell `cell-11` | BindCraft subprocess | **HIGH** | Crash | `subprocess.run(...)` without `check=True` or returncode check — BindCraft failure is silently ignored, then `glob` returns empty list → `FileNotFoundError` with misleading message |
| 18 | `NB02b` cell `cell-13` | MPNN chain extraction | **HIGH** | Wrong science | `parts[1]` assumes binder is always second in ProteinMPNN FASTA output; chain order can differ depending on input JSONL; no validation |
| 19 | `NB02b` cell `cell-15` | Protenix `run_protenix` | **HIGH** | Wrong science | `conf.get('ranking_score', 0.0)` — Protenix may not emit `ranking_score` directly; NB02c manually computes `0.8*iptm + 0.2*ptm` but NB02b trusts the JSON key, producing inconsistent scores between notebooks |
| 20 | `NB02b` cell `cell-17` | `epitope_contacts_from_cif` | **HIGH** | Wrong science | Uses `epitope_res_map` which is defined in `cell-7` config — if cell-7 is re-run after cell-17 is defined, dict is stale; also `epitope_res_map` uses 1-indexed keys but CIF residue numbers may differ |
| 21 | `NB02b` cell `cell-19` | output CSV | **MEDIUM** | Data flow | `pass_df` may be empty if all binders fail filter; `pass_df[out_cols].to_csv(...)` succeeds but NB03 then receives an empty CSV with no warning |
| 22 | `NB02b` cell `drive-setup` | path dependency | **MEDIUM** | Colab-specific | `NB01_BINDER_READY` assumes Drive path from NB01; if NB01 was run in a different session, this path may not exist and the upload fallback silently uses the first PDB in an empty glob |
| 23 | `NB02c` cell `cell-7` | `EPITOPE_RANGES` 0-indexed | **CRITICAL** | Wrong science | NB02c defines `EPITOPE_RANGES` as 0-indexed `(1, 7, ...)` / `(32, 38, ...)` — but `epitope_contacts_from_cif` in the same notebook slices `target_residues[start:end]` (Python list slice). For PHF6* this gives residues at list index 1–6, which is K18 residues 2–7 correctly only if chain A residue numbering starts at 1 AND list order matches residue order. Not verified. |
| 24 | `NB02c` cell `cell-17` | distance computation | **CRITICAL** | Wrong science | `(e_atom.coord - b_atom.coord).__abs__()` returns a `numpy.ndarray`, not a scalar. The `< cutoff` comparison on an array is ambiguous and will raise `ValueError: The truth value of an array is ambiguous` in NumPy ≥1.25 |
| 25 | `NB02c` cell `cell-12` | ESMFold loop indexing | **HIGH** | Wrong science | `for i, row in forge_df.iterrows()` — `i` is the DataFrame index, not a sequential counter; `if (i + 1) % 50 == 0` may never trigger if index is non-contiguous |
| 26 | `NB02c` cell `cell-14` | GPU release | **MEDIUM** | Resource | `del esmfold; torch.cuda.empty_cache()` — `del` removes the Python reference but if `esmfold` is still referenced by any closure or variable in Jupyter's namespace, VRAM is not freed |
| 27 | `NB02c` cell `cell-10` | Forge API | **MEDIUM** | Colab reliability | `wrapper.generate_binder(target_sequence=..., binder_length=..., n_samples=...)` — parameter names are inferred from wrapper usage; actual `InferenceWrapper` API may differ from `yk0/forge-e80`, causing `TypeError` |
| 28 | `NB03` cell `cell-5` | missing `subprocess` import | **CRITICAL** | Crash | `subprocess` is used in `cell-5` (`subprocess.run(['./foldmason/bin/foldmason', 'version', ...])`) but is never imported in NB03 — `NameError` at runtime |
| 29 | `NB03` cell `cell-16` | FoldMason JSON schema | **HIGH** | Wrong science | `fm_data['entries']`, `fm_data['scores']`, `fm_data['statistics']['msaLDDT']` — assumes a specific JSON schema. FoldMason `--report-mode 2` JSON structure is not officially documented; field names may differ between versions |
| 30 | `NB03` cell `cell-20` | HDBSCAN `all_points_membership_vectors` | **HIGH** | Crash | `hdbscan.all_points_membership_vectors(clusterer)` requires `prediction_data=True` to have been set **and** the model to have been fit on > 1 cluster. If all points are noise (cluster=-1), this raises `ValueError` |
| 31 | `NB03` cell `cell-28` | `known_proximity` merge | **HIGH** | Wrong science | `all_entries` is merged with `emb_df` on `name`, but `emb_df` was built from `valid_rows` which may have dropped entries (no AA sequence). Merged rows get `NaN` for `known_proximity`, which is then `fillna(0.5)` — silently neutral-scores binders that failed embedding |
| 32 | `NB03` cell `cell-7` | path assumptions | **HIGH** | Colab-specific | `FORGE_PDB_DIR = 'results/esmfold/forge'` and `TARGET_PDB = 'tau_k18_binder_ready/rank01_filtered_0001.pdb'` are relative paths; NB03 never sets `%cd`, so these resolve against Colab's default `/content/` — will not find Drive-based outputs from NB01/NB02 |
| 33 | `NB03` cell `cell-26` | chemical score formula | **MEDIUM** | Wrong science | `chemical_score` uses a simple flag-count / total formula, but `instability_index > 40` is a loose flag from a 1991 heuristic. More critically, the score is applied equally to control binders (`iptm=1.0`) which artificially inflates their integrated score |
| 34 | `NB03` cell `cell-19` | SaProt truncation | **MEDIUM** | Wrong science | `max_length=512` truncation will silently cut SaProt interleaved tokens. For a 130 aa binder, the interleaved format `'K v Q a ...'` produces 259 tokens (2 per residue + special tokens), which fits in 512. But for sequences >254 aa the embedding is of a truncated sequence with no warning |
| 35 | `environment.yml` | pytorch channel | **MEDIUM** | Environment | `pytorch>=2.0` is listed under `conda-forge` but PyTorch's official conda package is on the `pytorch` channel; `conda-forge` ships a CPU-only build, causing silent CPU fallback for all GPU ops in local testing |
| 36 | `requirements.txt` | `fair-esm>=2.0.0` | **MEDIUM** | Environment | `fair-esm` 2.x requires `torch<2.0` for some operations; pairing with `torch>=2.0` can produce silent numerical differences or import errors depending on resolver order |

---

## Detailed Findings

---

### `scripts/run_nb01_local.py`

#### BUG-01 — CRITICAL | Wrong Science
**Lines 217–224 — AFRC distance map misses last residue**

```python
for i in range(1, n):          # i goes 1 … n-1
    for j in range(i + 1, n): # j goes up to n-1
        afrc_dist[i-1, j-1] = mean_d
        afrc_dist[j-1, i-1] = mean_d
```

When `n = 130` (K18 length), `i` runs from 1 to 129 and `j` runs from `i+1` to 129. The largest pair is `(128, 129)` → stored at `afrc_dist[127, 128]`. The pair `(129, 130)` i.e. `i=129, j=130` never executes because `j < n = 130`. This means `afrc_dist[128, :]` and `afrc_dist[:, 128]` (0-indexed last row/column, corresponding to the final residue K18 position 130) remain all zeros.

The `deviation` formula divides by `ad` and uses `np.where(ad > 0, ...)` to suppress zero-denominator positions. This means the last residue is mapped to deviation=0 regardless of what the observed distance map shows there. The per-residue `marginal_dev` plot will show an artificial zero at the C-terminus, hiding any real structural signal. This directly affects the "trim guide" interpretation.

**The same bug is present identically in NB01 cell `f0b4cf53` (the production notebook).**

Fix: change outer loop to `range(1, n+1)` and inner to `range(i+1, n+1)` then store at `afrc_dist[i-1, j-1]` — but verify AFRC accepts residue index `n` (check library docs for whether `get_interresidue_distance_distribution` is 1-indexed with inclusive or exclusive upper bound).

---

#### BUG-02 — CRITICAL | Crash
**Line 246 — Division by zero / IndexError when all frames fail the Rg filter**

```python
positions = np.array([ca_arrays[i] for i in keep_idx])
obs_dist  = np.zeros((n_residues, n_residues))
for k in range(len(positions)):
    ...
obs_dist /= len(positions)    # ZeroDivisionError if keep_idx is empty
```

If the Rg filter passes 0 frames (e.g. `rg_std = 0` for a degenerate ensemble, or all frames are genuinely outliers), `keep_idx` is empty. `np.array([])` produces a shape-`(0,)` array, the loop body never runs, and `obs_dist /= 0` raises `ZeroDivisionError`. The same pattern exists in NB01 cell `f0b4cf53`.

Additionally, if `frame_files` is empty (e.g. FoldMason temp dir collision), `n_residues` stays `None`, causing `np.zeros((None, None))` → `TypeError`.

---

#### BUG-03 — HIGH | Wrong Output
**Line 430 — Output filenames use DataFrame index, not sequential rank**

```python
for rank_idx, row in df.iloc[:TOP_N].iterrows():
    ...
    dst = os.path.join(binder_dir, f'rank{rank_idx:02d}_{os.path.basename(src)}')
```

`iterrows()` yields the DataFrame's **index** as `rank_idx`. After `df.index += 1` (line 374), `iloc[:TOP_N]` iterates indices 1, 2, 3, 4, 5 — so filenames are `rank01_…` through `rank05_…`. This happens to be correct for a clean run, but if the DataFrame index is ever non-contiguous (e.g. after filtering, concatenation, or a re-run), the rank numbers in filenames will be wrong. NB03 `cell-7` hardcodes `rank01_filtered_0001.pdb` as `TARGET_PDB`, making this a latent dependency.

A safer pattern is `for rank_idx, (_, row) in enumerate(df.iloc[:TOP_N].iterrows(), start=1)`.

---

#### BUG-04 — MEDIUM | Crash
**Line 176 — FoldMason version check not validated**

```python
r = subprocess.run([FM_BIN, 'version'], capture_output=True, text=True)
print(f'FoldMason  : {r.stdout.strip()[:60]}')
```

`r.returncode` is never checked. If `FM_BIN` exists on disk but is not executable (e.g. wrong architecture, missing AVX2), this silently prints an empty string and the script continues until FoldMason actually fails at step 6, far from the root cause.

---

#### BUG-05 — MEDIUM | Security / Correctness
**Line 189 — `shell=True` with unsanitized sequence string**

```python
cmd = f'starling {SEQUENCE} --outname {PROTEIN_NAME} -c {args.n_conformers} -r'
result = subprocess.run(cmd, shell=True, cwd=str(out), ...)
```

`SEQUENCE` is hardcoded, so there is no immediate injection risk. However, `shell=True` with a formatted string means any whitespace or shell metacharacter in the sequence (e.g. if a user pastes a sequence with a newline) would cause silent misparse or command injection. The local script should use a list-form call: `subprocess.run(['starling', SEQUENCE, '--outname', PROTEIN_NAME, ...])`.

---

### `scripts/run_foldmason_k18.sh`

#### BUG-06 — HIGH | Wrong Output
**Lines 90–108 — JSON rename logic is fragile and version-dependent**

```bash
JSON_SRC="${OUTPUT_DIR}/msa_report.json"
JSON_DST="${OUTPUT_DIR}/foldmason.json"
if [[ -f "$JSON_SRC" ]] && [[ ! -f "$JSON_DST" ]]; then
    cp "$JSON_SRC" "$JSON_DST"
...
else
    JSON_FOUND=$(find "$OUTPUT_DIR" -name "*.json" | head -1)
```

The FoldMason JSON output filename changes between versions (`report.json`, `msa_report.json`, `msa.json` have all been observed in different releases). The fallback `find ... | head -1` is non-deterministic if multiple JSON files exist (e.g. from a prior partial run). If the wrong JSON is copied, NB02/NB03 will silently ingest stale or wrong data.

The script should enumerate the expected name from `--report-mode 2` documentation or check the FoldMason version and branch accordingly. At minimum, the fallback should error if more than one JSON is found.

---

#### BUG-07 — HIGH | Crash
**Line 32 — Path with space and parentheses**

```bash
FRAMES_DIR="${PROJECT_ROOT}/data/constructs/tau/tau_k18_starling_results (1)/tau_k18_frames"
```

This path contains a literal space and parentheses. While it is double-quoted here, any caller that uses `$FRAMES_DIR` without quoting will word-split. Additionally, FoldMason receives the directory as a positional argument on line 76 — it is quoted there, so FoldMason itself is safe. However, the `find "$FRAMES_DIR" -name "frame_*.pdb"` on line 52 is correctly quoted. The risk is primarily for downstream scripts or documentation examples that copy this path unquoted.

More importantly: the directory name contains `(1)` which suggests it was created by an OS download deduplication. The pipeline depends on this non-reproducible artifact path. It should be renamed and the canonical path documented.

---

#### BUG-08 — MEDIUM | Edge Case
**Line 52 — Frame count via `find | wc -l`**

```bash
N_FRAMES=$(find "$FRAMES_DIR" -name "frame_*.pdb" | wc -l)
```

`find` includes files in subdirectories by default (`-maxdepth` not set). If FoldMason or a prior run created subdirectories containing `frame_*.pdb` files, `N_FRAMES` will be inflated and the preflight check will pass falsely. Add `-maxdepth 1` to restrict to top-level files.

---

### `colab/NB01_ConformerEnsemble/01_TauK18_Conformer_Ensemble_Pipeline.ipynb`

#### BUG-09 — CRITICAL | Wrong Science (cell `f0b4cf53`)
Same as BUG-01. The AFRC distance map loop is identical. Last residue row/col = 0, creating a false zero-deviation signal at the C-terminus. See BUG-01 for full analysis.

#### BUG-10 — CRITICAL | Crash (cell `f0b4cf53`)
Same as BUG-02. If `keep_idx` is empty, `obs_dist /= len(positions)` is `obs_dist /= 0` → `ZeroDivisionError`.

#### BUG-11 — HIGH | Wrong Science (cell `02d9bbe9`)
**`--report-mode 1` in notebook vs `--report-mode 2` in shell script**

```python
cmd = [FM_BIN, 'easy-msa', *pdb_files, fm_msa_out, fm_msa_tmp,
       '--report-mode', '1', '--threads', '2']
```

`--report-mode 1` generates only the HTML report (used for lDDT parsing via regex). `--report-mode 2` generates the JSON file used by `run_foldmason_k18.sh` and expected by NB03. The notebook uses mode 1 and extracts lDDT via `parse_lddt_from_html`, while the shell script uses mode 2 and produces a JSON. If the HTML regex fails (BUG-12), all lDDT scores are silently set to `1.0` (fallback in cell `fo45m38jvqh`), making lDDT meaningless in ranking. Using `--report-mode 2` consistently would give a structured, parseable JSON and eliminate the fragile regex.

#### BUG-12 — HIGH | Wrong Science (cell `694567f4`)
**`parse_lddt_from_html` regex is fragile**

```python
pattern = r'"name"\s*:\s*"([^"]+)"[^}]*?"lddt"\s*:\s*([\d.]+)'
```

This regex relies on `"name"` and `"lddt"` appearing in the same JSON object within the HTML, with `"name"` always preceding `"lddt"` and no intervening `}`. If FoldMason emits the fields in a different order, or if there are nested objects, or if the HTML minification changes, this silently returns an empty dict. The fallback path sets all `lddt_norm` to `1.0`, meaning all conformers score equally on 60% of the combined ranking score — the lDDT component becomes meaningless and the ranking is driven entirely by epitope 3Di score.

#### BUG-13 — HIGH | Wrong Output (cell `c7f308cf`)
**PDB lookup by stem name may not match filtered filenames**

```python
stem = os.path.splitext(row['name'])[0]
candidates = (glob.glob(os.path.join(filtered_dir, f'{stem}.pdb')) +
              glob.glob(os.path.join(filtered_dir, f'*{stem}*.pdb')))
```

FoldMason names sequences in the MSA by the input PDB filename stem (e.g. `filtered_0047`). The filtered frames are named `filtered_0001.pdb`, `filtered_0002.pdb`, etc. If the stem lookup works at all, it is because FoldMason preserves the filename. But the glob `*{stem}*.pdb` is a substring match and could match multiple files if, for example, `filtered_0001.pdb` and `filtered_00010.pdb` both exist (with N_CONFORMERS ≥ 1000). The first candidate is taken silently.

#### BUG-14 — HIGH | Colab-specific (cell `7c954c77`)
**`%cd` magic creates a hard CWD dependency for all subsequent cells**

```python
%cd {_nb01_out}
```

All later cells use `os.path.abspath(...)` which resolves relative to the current working directory. If a user re-runs sections out of order (e.g. re-runs just the FoldMason section), they must re-run this cell first. If `NB01_OUT` points to a Drive path and Drive is not mounted (e.g. the Drive setup cell failed silently), `%cd` will create a local directory at the relative path. All files will be written to `/content/content/drive/...` literally, and the Drive copy cell at the end will fail.

A guard should be added: `assert os.path.isdir(_nb01_out), f'Output dir not ready: {_nb01_out}'`.

#### BUG-15 — MEDIUM | Colab Reliability (cell `10b1bb53`)
**Silent pip install failures**

```python
!pip install -q idptools-starling MDAnalysis afrc
```

`-q` suppresses all output including errors. If `idptools-starling` fails to compile (it has a Cython extension), the cell exits 0 but `import` fails three cells later with a cryptic `ModuleNotFoundError`. The `-q` flag should be replaced with `--quiet` only on the dependency noise level, or the install cell should check `import` immediately after.

---

### `colab/NB02b_BindCraft/02b_BindCraft_TauK18_Binder_Design.ipynb`

#### BUG-16 — CRITICAL | Wrong Science (cell `cell-7`)
**`EPITOPE_RANGES` format incompatible with NB01 and NB03**

In NB01 and NB03:
```python
EPITOPE_RANGES = [(2, 7, 'PHF6* VQIINK', 1.0), ...]  # (start, end, label, weight_float)
```

In NB02b:
```python
EPITOPE_RANGES = [(2, 7, 'PHF6* VQIINK', True), ...]  # (start, end, label, include_bool)
```

The 4th element is used differently. In NB02b `cell-7`, `True/False` controls whether the residue is included in `hotspot_res`. But `epitope_res_map` is built with `for r in range(start, end+1)` regardless of the `include` flag (all epitopes are added to the map). If NB02b code is ever consolidated with NB01/NB03 utility functions, passing a bool where a float weight is expected will silently compute `score = True * fraction` = `1.0 * fraction` for included epitopes — accidentally correct — but `False * fraction` = `0.0` would silently zero out `jR2R3` scoring in any shared code.

This inconsistency should be resolved into a single canonical format across all notebooks.

#### BUG-17 — HIGH | Crash (cell `cell-11`)
**BindCraft subprocess returncode not checked**

```python
subprocess.run(['python', '/content/BindCraft/bindcraft.py', '--settings', settings_path])
```

No `check=True`, no returncode inspection. If BindCraft crashes (common — PyRosetta licensing, AF2 weights missing, OOM on T4), execution continues to `glob.glob(f'{bc_out}/**/*.csv', ...)` which returns an empty list → `raise FileNotFoundError('No BindCraft CSV found')`. The error message points to the CSV, not to BindCraft failing, obscuring the true cause. Add `check=True` or check `result.returncode != 0` with a diagnostic print of `result.stderr`.

#### BUG-18 — HIGH | Wrong Science (cell `cell-13`)
**ProteinMPNN chain order assumed to be target=A, binder=B**

```python
parts = lines[i+1].strip().split('/')
binder_seq = parts[1] if len(parts) > 1 else parts[0]
```

ProteinMPNN outputs the full complex sequence separated by `/`, ordered by the chain order in the input JSONL. The JSONL in this cell sets `fixed_chains: [TARGET_CHAIN]` and `designed_chains: [BINDER_CHAIN]` with `TARGET_CHAIN = 'A'` and `BINDER_CHAIN = 'B'`. ProteinMPNN's output order depends on how it internally sorts chains — not guaranteed to match JSONL order. If chain B is output first, `parts[1]` extracts the target sequence (K18) instead of the binder, and that K18 sequence is submitted to Protenix as the "binder," silently producing completely wrong complex predictions.

Validation: after extraction, assert `len(binder_seq)` is within the configured binder length range `[BINDER_LEN_MIN, BINDER_LEN_MAX]`.

#### BUG-19 — HIGH | Wrong Science (cell `cell-15`)
**Protenix `ranking_score` key assumed to exist in confidence JSON**

In NB02b `run_protenix`:
```python
return (conf.get('iptm', 0.0),
        conf.get('ranking_score', 0.0),
        cif_files[0] if cif_files else None)
```

In NB02c `run_protenix`:
```python
iptm    = conf.get('iptm', 0.0)
ptm     = conf.get('ptm',  0.0)
ranking = 0.8 * iptm + 0.2 * ptm
```

NB02b trusts a `ranking_score` key in Protenix's JSON; NB02c computes it manually from `iptm` and `ptm`. If Protenix does not emit `ranking_score` in its confidence JSON (which is the case for the `protenix_mini` model), NB02b silently returns `0.0` for all `ranking_score` values, and every binder fails the `MIN_RANKING_SCORE = 0.5` filter. The notebook produces an empty `pass_df` and an empty CSV — passed to NB03 without any error.

#### BUG-20 — HIGH | Wrong Science (cell `cell-17`)
**`epitope_res_map` uses 1-indexed keys but CIF residue sequence numbers may differ**

```python
lbl = epitope_res_map.get(t_res)
```

`t_res` is `res.get_id()[1]` from BioPython — the residue **sequence number** in the CIF file, which for Protenix outputs is typically 1-indexed from 1. `epitope_res_map` is built from `range(start, end+1)` using the `EPITOPE_RANGES` start/end values (1-indexed K18 positions). This is correct **if and only if** Protenix numbers chain A residues starting from 1 with no gaps. If Protenix uses 0-indexed residues, or if there are insertion codes, all `epitope_res_map.get(t_res)` lookups return `None` and `epi_hits` stays empty — all binders report 0 epitope contacts and all are filtered out.

#### BUG-21 — MEDIUM | Data Flow (cell `cell-19`)
**Empty `pass_df` produces empty CSV without error, silently breaking NB03**

```python
pass_df = mpnn_df[mpnn_df['verdict']=='PASS'].copy()
...
pass_df[out_cols].to_csv(passing_csv, index=False)
```

If 0 binders pass the filter (e.g. due to BUG-19 above), a valid-looking CSV with headers but 0 rows is written. NB03 `cell-9` reads this CSV, computes `len(df)` = 0, prints the method's describe stats (which will be NaN-filled), and continues silently. The downstream embedding and ranking code will run on zero binders — all computations on empty DataFrames — and produce empty output files that look like successful runs.

#### BUG-22 — MEDIUM | Colab-specific (cell `drive-setup`)
**Drive path fallback logic uses first PDB from a potentially stale glob**

```python
_drive_pdbs = sorted(_glob.glob(os.path.join(NB01_BINDER_READY, '*.pdb')))
if _drive_pdbs:
    TARGET_PDB_PATH = _drive_pdbs[0]
```

If NB01 was run in a previous session with different conformers and the Drive directory was not cleared, `_drive_pdbs[0]` may be from an older run. The sorted glob picks alphabetically-first PDB, which may not be the best-ranked conformer. This will produce valid-looking BindCraft output on the wrong structure.

---

### `colab/NB02c_Forge/02c_Forge_TauK18_Binder_Design.ipynb`

#### BUG-23 — CRITICAL | Wrong Science (cell `cell-17`)
**`epitope_contacts_from_cif` indexes residues by list position, not sequence number**

```python
target_residues = list(chain_A.get_residues())
...
for res in target_residues[start:end]:
    epi_atoms.extend(res.get_atoms())
```

`EPITOPE_RANGES` in NB02c uses 0-indexed values: `(1, 7, 'PHF6star_VQIINK', 1.0)`. The slice `target_residues[1:7]` gives list positions 1 through 6. This happens to match K18 residues 2–7 **only if** chain A has no HETATM, no missing residues, no insertion codes, and Protenix writes residues starting at residue number 1 in sequence order. In any other case, the slice selects the wrong residues. The correct approach is to filter by residue sequence number, not list index.

#### BUG-24 — CRITICAL | Crash (cell `cell-17`)
**`(e_atom.coord - b_atom.coord).__abs__()` is not a scalar distance**

```python
if (e_atom.coord - b_atom.coord).__abs__() < cutoff:
```

`e_atom.coord` and `b_atom.coord` are `numpy.ndarray` of shape `(3,)`. Subtraction gives a `(3,)` array. `__abs__()` on a NumPy array is `np.abs(array)` — which returns a `(3,)` array of absolute values of each component, **not the Euclidean norm**. Comparing a `(3,)` array to a scalar float with `<` triggers `ValueError: The truth value of an array with more than one element is ambiguous` in NumPy ≥ 1.25 (which is the version pinned in `requirements.txt`). This cell will crash at runtime for every CIF processed.

The correct computation is `np.linalg.norm(e_atom.coord - b_atom.coord) < cutoff`.

#### BUG-25 — HIGH | Wrong Science (cell `cell-12`)
**`iterrows()` progress counter uses non-sequential index**

```python
for i, row in forge_df.iterrows():
    ...
    if (i + 1) % 50 == 0:
        print(f'  {i+1}/{len(forge_df)}')
```

After `pd.DataFrame(forge_records)`, the DataFrame index is 0-based and sequential (0, 1, 2, …) — so `i` starts at 0. `(0+1) % 50 == 0` is False, `(49+1) % 50 == 0` is True. This is actually correct for a fresh DataFrame. However, if `forge_df` is ever filtered or reset between cells (e.g. if the user re-runs only certain cells), the index becomes non-contiguous and the progress print either never triggers or triggers at wrong intervals. Use `enumerate` instead of relying on DataFrame index for progress counting.

#### BUG-26 — MEDIUM | Resource (cell `cell-14`)
**GPU memory not reliably freed by `del esmfold`**

```python
del esmfold
torch.cuda.empty_cache()
```

If any other variable in the Jupyter session holds a reference to the ESMFold model or its parameters (e.g. through a closure in `fold_and_plddt`), `del esmfold` only removes one name binding and the model stays in VRAM. On a T4 (16 GB), ESMFold occupies ~12 GB. Protenix mini requires ~8 GB. Running both simultaneously causes OOM. A safe pattern is to also `del fold_and_plddt` (the closure captures `esmfold`) or define the function inside a scope that is explicitly deleted.

#### BUG-27 — MEDIUM | Colab Reliability (cell `cell-10`)
**Forge `InferenceWrapper.generate_binder` API not verified**

```python
seqs = wrapper.generate_binder(
    target_sequence=TARGET_SEQUENCE,
    binder_length=binder_len,
    n_samples=N_SAMPLES,
)
```

The `yk0/forge-e80` HuggingFace repo's `InferenceWrapper` API is not documented in the notebook. If the actual method signature is `generate(...)` or uses positional arguments, this will raise `TypeError`. There is no try/except around this call, so a wrong API call will abort the entire Forge generation step with no useful guidance for recovery.

---

### `colab/NB03_BinderComparison/03_Binder_Comparison_Analysis.ipynb`

#### BUG-28 — CRITICAL | Crash (cell `cell-5`)
**`subprocess` is never imported in NB03**

```python
# cell-5
r = subprocess.run(['./foldmason/bin/foldmason', 'version'], capture_output=True, text=True)
```

The imports in `cell-5` are:
```python
import os, json, glob, shutil, zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import MDAnalysis as mda
```

`subprocess` is not imported anywhere in NB03. This cell will raise `NameError: name 'subprocess' is not defined`. It is also used later in `cell-15` (FoldMason run) and `cell-11` (ESMFold install check). All of these will fail.

#### BUG-29 — HIGH | Wrong Science (cell `cell-16`)
**FoldMason JSON schema assumed without validation**

```python
entries   = fm_data['entries']
fm_scores = np.array(fm_data['scores'])
msa_lddt  = fm_data['statistics']['msaLDDT']
...
di3_map[stem] = entry['ss']
lddt_map[stem] = score
aa_map[stem]   = entry['aa']
```

The FoldMason `--report-mode 2` JSON schema is not documented in the FoldMason paper or README. Field names `entries`, `scores`, `statistics.msaLDDT`, `entry['ss']`, `entry['aa']` are inferred from one observed output. If the FoldMason version on Colab differs from the one used to develop this notebook (very likely — `!wget` fetches the latest), any of these keys may be missing or renamed, causing `KeyError` with no fallback. At minimum, add `if 'entries' not in fm_data: raise RuntimeError(...)` with a message pointing to the actual JSON keys found.

#### BUG-30 — HIGH | Crash (cell `cell-20`)
**`hdbscan.all_points_membership_vectors` fails when all points are noise**

```python
cluster_probs = hdbscan.all_points_membership_vectors(clusterer).max(axis=1)
```

`all_points_membership_vectors` requires `prediction_data=True` (set on line above ✓) but also raises if there are zero clusters (all points assigned to cluster `-1`). With a small binder set (e.g. 10–20 binders from a short run), HDBSCAN with `min_cluster_size=10` may legitimately produce no clusters. This will raise `ValueError: No clusters found; prediction is not possible`. The integrated score then has no `cluster_prob` for any binder, causing `fillna(0.0)` in `cell-28` to silently zero out 15% of the integrated score for all binders.

Add a guard: check `len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0) > 0` before calling `all_points_membership_vectors`.

#### BUG-31 — HIGH | Wrong Science (cell `cell-28`)
**`known_proximity` NaN fill for unembedded binders is silently neutral**

```python
all_entries = all_entries.merge(emb_df[merge_cols], on='name', how='left')
...
W['known_proximity'] * all_entries['known_proximity'].fillna(0.5)
```

Binders that failed the embedding step (no AA sequence, wrong length) are not in `emb_df`. After the left merge, their `known_proximity` is `NaN`, filled with `0.5`. This means failed embeddings score 50th percentile on known proximity — indistinguishable from binders that are genuinely moderately similar to controls. Failed embeddings should be explicitly flagged or given `known_proximity = 0.0` (or excluded).

#### BUG-32 — HIGH | Colab-specific (cell `cell-7`)
**Relative paths do not resolve to Drive outputs**

```python
FORGE_PDB_DIR     = 'results/esmfold/forge'
TARGET_PDB        = 'tau_k18_binder_ready/rank01_filtered_0001.pdb'
BINDCRAFT_PDB_DIR = 'results/esmfold/bindcraft'
```

NB03 never issues `%cd` or `os.chdir()`. The Colab default CWD is `/content/`. These relative paths resolve to `/content/results/esmfold/forge` etc. But NB01 writes to `NB01_OUT = /content/drive/MyDrive/CHEM_280/results/nb01_output/`, and NB02b/02c write to corresponding Drive paths. Nothing in NB03 copies these Drive outputs to `/content/`. The Drive setup cell does copy CSVs from Drive to `/content/`, but does **not** copy PDB directories. The result is that `cell-12` (folding loop) will never find existing PDBs, will re-fold everything from scratch, and will write the new PDBs to `/content/results/esmfold/...` which is an ephemeral path not backed by Drive.

#### BUG-33 — MEDIUM | Wrong Science (cell `cell-26`)
**Chemical scoring applied equally to known controls inflates their integrated score**

```python
known_meta['iptm'] = 1.0   # known binders = reference quality
```

Known controls get `iptm = 1.0` (hardcoded) plus real chemical scores computed from their sequences. Their integrated score will naturally dominate the top of the ranking, which is the intended use for proximity reference. However, because `known_proximity` for a known control to itself is 1.0 (cosine similarity with the centroid it helps define), and the centroid is biased by the controls themselves, the known proximity score is not an independent signal — it partially double-counts the controls' own embedding. If controls cluster together and a designed binder lands near them, it will score high on proximity. But if the controls are poor binders biochemically (as some tau antibodies are for K18 specifically), this creates misleading signal.

#### BUG-34 — MEDIUM | Wrong Science (cell `cell-19`)
**SaProt `max_length=512` silently truncates long sequences**

```python
inputs = tokenizer(batch, return_tensors='pt', padding=True,
                   truncation=True, max_length=512).to('cuda')
```

SaProt interleaved format produces approximately 2 tokens per residue plus 2 special tokens. For 130 aa (K18 length): `130*2 + 2 = 262 tokens` — safe. For 80 aa binders: `80*2 + 2 = 162 tokens` — safe. However, the K18 target sequence itself (130 aa) is also embedded as part of the pipeline. The critical concern is if any binder sequence exceeds ~255 aa (unlikely given BINDER_LEN_MAX=80), but the truncation is silent and produces an embedding for a different-length sequence than intended. A warning should be emitted when `len(tokenized_input) > 500`.

---

### `environment.yml`

#### BUG-35 — MEDIUM | Environment
**PyTorch listed under `conda-forge` channel**

```yaml
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - pytorch>=2.0
```

The `pytorch` package under `conda-forge` is a CPU-only build. The official CUDA-enabled PyTorch is distributed from the `pytorch` channel (`https://conda.anaconda.org/pytorch`). In local testing, `conda env create -f environment.yml` will install CPU PyTorch, causing all GPU operations to silently fall back to CPU (no error, just 10–100x slower). Add `- pytorch` to the channels list before `conda-forge`.

---

### `requirements.txt`

#### BUG-36 — MEDIUM | Environment
**`fair-esm>=2.0.0` and `torch>=2.0` may conflict**

```
torch>=2.0
fair-esm>=2.0.0
```

`fair-esm` 2.0.0 was released before `torch` 2.0 and has internal imports that use deprecated or removed PyTorch APIs (specifically `torch.load` behavior changes in 2.0, and `torch.nn.utils.rnn` changes in 2.1). While recent versions of `fair-esm` (2.0.1+) have addressed some of these, the unconstrained `>=2.0.0` on fair-esm paired with `>=2.0` on torch means pip's resolver may pick any combination. Pin to `fair-esm==2.0.0` with `torch>=2.0,<2.2` or test with the specific versions used in development.

---

## Cross-Notebook Data Flow Issues

### Epitope range format inconsistency (BUG-16 cross-ref)

The `EPITOPE_RANGES` variable is redefined independently in every notebook with different conventions:

| Notebook | Format | 4th element |
|----------|--------|-------------|
| NB01 | `(1-indexed start, 1-indexed end, label, float weight)` | weight |
| NB02b | `(1-indexed start, 1-indexed end, label, bool include)` | include flag |
| NB02c | `(0-indexed start, 0-indexed end, label, float weight)` | weight |
| NB03 | `(1-indexed start, 1-indexed end, label, float weight)` | weight |

This is a maintenance trap. Any function copy-pasted between notebooks silently operates on the wrong indexing convention or interprets the 4th element incorrectly.

**Recommendation:** define a canonical `EPITOPE_RANGES` in a shared config cell or a `config.py` file, and validate the format at the top of each notebook with an assertion.

### `--report-mode` inconsistency (BUG-11 cross-ref)

`run_foldmason_k18.sh` uses `--report-mode 2` (JSON output), while NB01 uses `--report-mode 1` (HTML output). NB03 expects a JSON file. If NB01 is used as the production path (as intended), the JSON required by NB03 is never generated.

### Conformer filename coupling (BUG-13 / BUG-03 cross-ref)

NB03 `cell-7` hardcodes:
```python
TARGET_PDB = 'tau_k18_binder_ready/rank01_filtered_0001.pdb'
```

This assumes: (1) rank prefix is `rank01`, (2) the conformer is `filtered_0001`, (3) the file is named exactly this. BUG-03 shows the rank prefix depends on the DataFrame index, not a guaranteed counter. If the top conformer is `filtered_0047`, the NB03 target PDB will silently be the wrong structure.

---

## Priority Fix Order

1. **BUG-24** (NB02c `__abs__()` crash) — will crash 100% of runs, must fix first
2. **BUG-28** (NB03 missing `subprocess` import) — will crash 100% of NB03 runs
3. **BUG-01/09** (AFRC last residue zero) — silently wrong science in every run
4. **BUG-02/10** (empty keep_idx crash) — probable on degenerate ensembles
5. **BUG-16** (EPITOPE_RANGES format inconsistency) — cross-notebook correctness
6. **BUG-19** (Protenix ranking_score key missing) — silent empty output from NB02b
7. **BUG-11** (--report-mode inconsistency NB01 vs shell script) — NB03 gets no JSON
8. **BUG-30** (HDBSCAN crash on no clusters) — likely with small binder sets
9. **BUG-32** (NB03 relative paths don't resolve to Drive) — PDB folding loop broken
10. All remaining HIGH issues before production run

---

*End of review. Total issues found: 36 (6 CRITICAL, 17 HIGH, 11 MEDIUM, 2 LOW).*
