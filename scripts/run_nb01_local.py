#!/usr/bin/env python3
"""
Local test runner for NB01 — TauK18 Conformer Ensemble Pipeline.

Mirrors the notebook logic end-to-end without any Colab-specific calls
(no %cd magic, no files.download, no google.colab imports).

Usage:
    python scripts/run_nb01_local.py [--n_conformers 50] [--out_dir results/nb01_test]

FoldMason binary is resolved in order:
  1. --foldmason argument
  2. FOLDMASON_BIN env var
  3. 'foldmason' on PATH
  4. /tmp/foldmason_local/foldmason/bin/foldmason (local dev download)
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import zipfile
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')           # headless — no display needed
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from afrc import AnalyticalFRC

# ── Config ────────────────────────────────────────────────────────────────────
PROTEIN_NAME = 'tau_k18'

SEQUENCE = (
    'KVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGG'
    'QVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPV'
    'VSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL'
)

EPITOPE_RANGES = [
    (2,  7,  'PHF6* VQIINK', 1.0),
    (33, 38, 'PHF6 VQIVYK',  1.0),
    (52, 70, 'jR2R3',         0.5),
]

TOP_N          = 5
LDDT_WEIGHT    = 0.6
EPITOPE_WEIGHT = 0.4


# ── Helpers ───────────────────────────────────────────────────────────────────

def resolve_foldmason(hint: str | None) -> str:
    candidates = []
    if hint:
        candidates.append(hint)
    if 'FOLDMASON_BIN' in os.environ:
        candidates.append(os.environ['FOLDMASON_BIN'])
    candidates.append('foldmason')
    candidates.append('/tmp/foldmason_local/foldmason/bin/foldmason')
    candidates.append('/content/foldmason/bin/foldmason')
    for c in candidates:
        if shutil.which(c) or os.path.isfile(c):
            return c
    raise FileNotFoundError(
        'FoldMason not found. Set FOLDMASON_BIN, pass --foldmason, or '
        'run: wget https://mmseqs.com/foldmason/foldmason-linux-avx2.tar.gz'
    )


def extract_frames(topology, trajectory, output_dir, n_frames):
    os.makedirs(output_dir, exist_ok=True)
    u = mda.Universe(topology, trajectory)
    total   = len(u.trajectory)
    indices = np.linspace(0, total - 1, n_frames, dtype=int)
    saved   = []
    for i, ts_idx in enumerate(indices):
        u.trajectory[ts_idx]
        out = os.path.join(output_dir, f'frame_{i+1:04d}.pdb')
        u.atoms.write(out)
        saved.append(out)
    print(f'  Extracted {len(saved)} frames → {output_dir}/')
    return saved


def parse_fasta_msa(path):
    records = {}
    header, lines = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    records[header] = ''.join(lines)
                header = line[1:].split()[0]
                lines = []
            else:
                lines.append(line)
    if header is not None:
        records[header] = ''.join(lines)
    return records


def compute_msa_consistency(records):
    names = list(records.keys())
    seqs  = list(records.values())
    n_cols = len(seqs[0])
    consensus = []
    for col in range(n_cols):
        chars   = [s[col] for s in seqs if col < len(s)]
        non_gap = [c for c in chars if c != '-']
        consensus.append(Counter(non_gap).most_common(1)[0][0] if non_gap else '-')
    scores = {}
    for name, seq in zip(names, seqs):
        matches = total = 0
        for char, cons in zip(seq, consensus):
            if char != '-' and cons != '-':
                total += 1
                if char == cons:
                    matches += 1
        scores[name] = matches / total if total > 0 else 0.0
    return scores, consensus


def parse_lddt_from_html(html_path):
    scores = {}
    if not os.path.exists(html_path):
        return scores
    with open(html_path) as f:
        content = f.read()
    pattern = r'"name"\s*:\s*"([^"]+)"[^}]*?"lddt"\s*:\s*([\d.]+)'
    for m in re.finditer(pattern, content):
        scores[m.group(1)] = float(m.group(2))
    return scores


def epitope_3di_score(di_seq, di_consensus):
    epi_fracs, epi_weights = [], []
    for start, end, label, weight in EPITOPE_RANGES:
        matches = total = 0
        for p in range(start - 1, end):
            if p < len(di_seq) and p < len(di_consensus):
                c, con = di_seq[p], di_consensus[p]
                if c != '-' and con != '-':
                    total += 1
                    if c == con:
                        matches += 1
        epi_fracs.append(matches / total if total > 0 else 0.0)
        epi_weights.append(weight)
    total_w = sum(epi_weights)
    score   = sum(f * w for f, w in zip(epi_fracs, epi_weights)) / total_w
    return score, epi_fracs


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n_conformers', type=int, default=50)
    ap.add_argument('--out_dir', default='results/nb01_test')
    ap.add_argument('--foldmason', default=None)
    args = ap.parse_args()

    out = Path(args.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    print(f'\n=== NB01 local test run  ({args.n_conformers} conformers) ===')
    print(f'Output dir : {out}\n')

    FM_BIN = resolve_foldmason(args.foldmason)
    r = subprocess.run([FM_BIN, 'version'], capture_output=True, text=True)
    print(f'FoldMason  : {r.stdout.strip()[:60]}')
    print(f'Sequence   : {len(SEQUENCE)} aa')
    print()

    # ── 1. Starling ───────────────────────────────────────────────────────────
    print('--- Step 1: Starling ---')
    pdb_file = str(out / f'{PROTEIN_NAME}_STARLING.pdb')
    xtc_file = str(out / f'{PROTEIN_NAME}_STARLING.xtc')

    if os.path.exists(pdb_file) and os.path.exists(xtc_file):
        print('  Starling output already exists, skipping.')
    else:
        cmd = f'starling {SEQUENCE} --outname {PROTEIN_NAME} -c {args.n_conformers} -r'
        print(f'  Running: {cmd[:80]}...')
        result = subprocess.run(cmd, shell=True, cwd=str(out), capture_output=True, text=True)
        if result.returncode != 0:
            print(result.stderr[-600:])
            sys.exit(1)

    assert os.path.exists(pdb_file), f'Missing: {pdb_file}'
    assert os.path.exists(xtc_file), f'Missing: {xtc_file}'
    print(f'  PDB: {os.path.getsize(pdb_file)//1024} KB   '
          f'XTC: {os.path.getsize(xtc_file)//1024} KB')

    # ── 2. Extract frames ─────────────────────────────────────────────────────
    print('\n--- Step 2: Extract frames ---')
    frames_dir = str(out / f'{PROTEIN_NAME}_frames')
    frames = extract_frames(pdb_file, xtc_file, frames_dir, n_frames=args.n_conformers)

    # ── 3. AFRC filter ────────────────────────────────────────────────────────
    print('\n--- Step 3: AFRC filter ---')
    afrc    = AnalyticalFRC(SEQUENCE)
    rg_afrc = afrc.get_mean_radius_of_gyration()
    re_afrc = afrc.get_mean_end_to_end_distance()
    print(f'  AFRC reference:  Rg={rg_afrc:.1f} Å   Re={re_afrc:.1f} Å')

    print('  Building AFRC distance map...')
    n = len(SEQUENCE)
    # AFRC uses 1-indexed residues capped at N-1 (bond indexing: N residues → N-1 bonds)
    afrc_dist = np.zeros((n, n))
    for i in range(1, n):
        for j in range(i + 1, n):
            bins, p_r = afrc.get_interresidue_distance_distribution(i, j)
            bins = np.array(bins);  p_r = np.array(p_r)
            p_r /= p_r.sum()
            mean_d = float(np.dot(bins, p_r))
            afrc_dist[i-1, j-1] = mean_d
            afrc_dist[j-1, i-1] = mean_d
    print('  Done.')

    frame_files = sorted(glob.glob(os.path.join(frames_dir, 'frame_*.pdb')))
    rg_all, ca_arrays = [], []
    n_residues = None
    for fn in frame_files:
        u_f = mda.Universe(fn)
        ca_f = u_f.select_atoms('name CA')
        if n_residues is None:
            n_residues = len(ca_f)
        rg_all.append(ca_f.radius_of_gyration())
        ca_arrays.append(ca_f.positions.copy())
    rg_all = np.array(rg_all)

    rg_mean, rg_std = rg_all.mean(), rg_all.std()
    rg_lo, rg_hi    = rg_mean - 2*rg_std, rg_mean + 2*rg_std
    keep_idx        = np.where((rg_all >= rg_lo) & (rg_all <= rg_hi))[0]
    print(f'  Rg: mean={rg_mean:.1f}  std={rg_std:.1f}  AFRC={rg_afrc:.1f} Å')
    print(f'  Passing ±2σ: {len(keep_idx)}/{len(frame_files)}')

    # ── 4. AFRC deviation plots ───────────────────────────────────────────────
    positions = np.array([ca_arrays[i] for i in keep_idx])
if len(positions) == 0:
    print("WARNING: 0 frames passed Rg filter — check ensemble quality. Using all frames.")
    positions = np.array(ca_arrays)
    keep_idx = np.arange(len(ca_arrays))
    obs_dist  = np.zeros((n_residues, n_residues))
    for k in range(len(positions)):
        diff       = positions[k][:, None, :] - positions[k][None, :, :]
        obs_dist  += np.linalg.norm(diff, axis=-1)
    obs_dist /= len(positions)
    ad = afrc_dist[:n_residues, :n_residues]
    with np.errstate(invalid='ignore', divide='ignore'):
        deviation = np.where(ad > 0, (obs_dist - ad) / ad, 0.0)
    np.fill_diagonal(deviation, 0)
    marginal_dev = np.abs(deviation).mean(axis=1)

    epi_colors_map = {'PHF6* VQIINK': 'royalblue', 'PHF6 VQIVYK': 'crimson', 'jR2R3': 'forestgreen'}
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    ax = axes[0, 0]
    ax.hist(rg_all, bins=30, color='steelblue', edgecolor='white', alpha=0.8, label='All')
    ax.hist(rg_all[keep_idx], bins=30, color='seagreen', edgecolor='white', alpha=0.6, label='Passing')
    ax.axvline(rg_mean, color='steelblue', ls='--', lw=1.5, label=f'Obs={rg_mean:.1f} Å')
    ax.axvline(rg_afrc, color='tomato', ls='-', lw=2, label=f'AFRC={rg_afrc:.1f} Å')
    ax.set_xlabel('Rg (Å)'); ax.set_title('Rg vs AFRC reference'); ax.legend(fontsize=7)
    ax = axes[0, 1]
    ax.plot(range(1, n_residues+1), marginal_dev, color='darkorange', lw=1.5)
    ax.axhline(marginal_dev.mean(), color='gray', ls='--', lw=1)
    for start, end, label, *_ in EPITOPE_RANGES:
        ax.axvspan(start, end, alpha=0.2, color=epi_colors_map.get(label, 'gray'), label=label)
    ax.set_xlabel('Residue'); ax.set_ylabel('Mean |deviation|')
    ax.set_title('Per-residue AFRC deviation'); ax.legend(fontsize=7)
    ax = axes[1, 0]
    im = ax.imshow(obs_dist, cmap='viridis_r', aspect='auto')
    plt.colorbar(im, ax=ax, label='Mean Cα–Cα (Å)'); ax.set_title('Observed Cα–Cα distance map')
    ax = axes[1, 1]
    vmax = max(0.3, float(np.abs(deviation).max()) * 0.8)
    im2  = ax.imshow(deviation, cmap='RdBu_r', aspect='auto', vmin=-vmax, vmax=vmax)
    plt.colorbar(im2, ax=ax, label='(obs−AFRC)/AFRC'); ax.set_title('Fractional AFRC deviation')
    plt.suptitle(f'{PROTEIN_NAME} — Ensemble vs AFRC', fontsize=13, y=1.01)
    plt.tight_layout()
    ensemble_png = str(out / f'{PROTEIN_NAME}_ensemble_analysis.png')
    plt.savefig(ensemble_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved {os.path.basename(ensemble_png)}')

    # ── 5. Write filtered frames ──────────────────────────────────────────────
    print('\n--- Step 5: Write filtered frames ---')
    filtered_dir = str(out / f'{PROTEIN_NAME}_frames_filtered')
    os.makedirs(filtered_dir, exist_ok=True)
    for rank, frame_idx in enumerate(keep_idx):
        src = frame_files[frame_idx]
        dst = os.path.join(filtered_dir, f'filtered_{rank+1:04d}.pdb')
        shutil.copy2(src, dst)
    n_filt = len(keep_idx)
    print(f'  {n_filt} filtered frames → {filtered_dir}/')

    # ── 6. FoldMason ─────────────────────────────────────────────────────────
    print('\n--- Step 6: FoldMason easy-msa ---')
    fm_msa_out = str(out / f'{PROTEIN_NAME}_foldmason_msa')
    fm_msa_tmp = str(out / f'{PROTEIN_NAME}_foldmason_tmp')
    os.makedirs(fm_msa_tmp, exist_ok=True)

    pdb_files_filt = sorted(glob.glob(os.path.join(filtered_dir, '*.pdb')))
    print(f'  {len(pdb_files_filt)} PDBs → FoldMason')
    cmd = [FM_BIN, 'easy-msa', *pdb_files_filt, fm_msa_out, fm_msa_tmp,
           '--report-mode', '1', '--threads', '4']
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f'  FoldMason STDERR:\n{r.stderr[-600:]}')
        sys.exit(1)

    msa_aa_path  = f'{fm_msa_out}_aa.fa'
    msa_3di_path = f'{fm_msa_out}_3di.fa'
    html_path    = f'{fm_msa_out}.html'
    for label, path in [('_aa.fa', msa_aa_path), ('_3di.fa', msa_3di_path), ('.html', html_path)]:
        print(f'  {label}: {"OK" if os.path.exists(path) else "MISSING"}')

    # ── 7. Parse FoldMason outputs ────────────────────────────────────────────
    print('\n--- Step 7: Parse FoldMason outputs ---')
    aa_records  = parse_fasta_msa(msa_aa_path)
    di_records  = parse_fasta_msa(msa_3di_path)
    consistency_scores, aa_consensus = compute_msa_consistency(aa_records)
    lddt_scores = parse_lddt_from_html(html_path)

    di_seqs_list = list(di_records.values())
    n_cols_3di   = len(di_seqs_list[0]) if di_seqs_list else 0
    di_consensus = []
    for col in range(n_cols_3di):
        chars   = [s[col] for s in di_seqs_list if col < len(s)]
        non_gap = [c for c in chars if c != '-']
        di_consensus.append(Counter(non_gap).most_common(1)[0][0] if non_gap else '-')

    print(f'  AA MSA  : {len(aa_records)} conformers × {len(aa_consensus)} cols')
    print(f'  3Di MSA : {len(di_records)} conformers × {n_cols_3di} cols')
    if lddt_scores:
        vals = list(lddt_scores.values())
        print(f'  lDDT    : {len(lddt_scores)} structures  '
              f'mean={np.mean(vals):.3f}  min={np.min(vals):.3f}  max={np.max(vals):.3f}')
    else:
        print('  lDDT    : not found in HTML — ranking by consistency + 3Di only')

    # ── 8. Score & rank ───────────────────────────────────────────────────────
    print('\n--- Step 8: Score and rank conformers ---')
    names_all = list(consistency_scores.keys())
    if lddt_scores:
        lddt_vals = np.array([lddt_scores.get(n, 0.0) for n in names_all])
        lo, hi    = lddt_vals.min(), lddt_vals.max()
        lddt_norm_map = {
            n: float((lddt_vals[i]-lo)/(hi-lo)) if hi > lo else 1.0
            for i, n in enumerate(names_all)
        }
    else:
        lddt_norm_map = {n: 1.0 for n in names_all}

    rows = []
    for name in names_all:
        cons     = consistency_scores[name]
        lddt_n   = lddt_norm_map[name]
        raw_lddt = lddt_scores.get(name, float('nan'))
        di_seq   = di_records.get(name, '')
        epitope_struct, epi_fracs = epitope_3di_score(di_seq, di_consensus)
        combined = LDDT_WEIGHT * lddt_n + EPITOPE_WEIGHT * epitope_struct
        row = {
            'name': name, 'lddt_raw': round(raw_lddt, 4) if not np.isnan(raw_lddt) else None,
            'lddt_norm': round(lddt_n, 4), 'consistency': round(cons, 4),
            'epitope_struct': round(epitope_struct, 3), 'combined': round(combined, 4),
        }
        for (start, end, label, weight), frac in zip(EPITOPE_RANGES, epi_fracs):
            row[f'frac_{label.split()[0]}'] = round(frac, 3)
        rows.append(row)

    df = pd.DataFrame(rows).sort_values('combined', ascending=False).reset_index(drop=True)
    df.index += 1
    show_cols = (['name', 'lddt_raw', 'consistency', 'epitope_struct', 'combined'] +
                 [f'frac_{label.split()[0]}' for _, _, label, _ in EPITOPE_RANGES])
    print(df[show_cols].head(10).to_string())

    # ── 9. FoldMason plots ────────────────────────────────────────────────────
    top_names = set(df.iloc[:TOP_N]['name'])
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    ax = axes[0]
    colors = ['tomato' if n in top_names else 'steelblue' for n in df['name']]
    ax.bar(range(len(df)), df['combined'].values, color=colors, width=1.0, edgecolor='none')
    ax.set_xlabel('Conformer'); ax.set_ylabel('Combined score')
    ax.set_title(f'Combined score (red = top {TOP_N})')
    ax = axes[1]
    sc = ax.scatter(df['lddt_norm'], df['epitope_struct'],
                    c=df['combined'], cmap='RdYlGn', s=40, edgecolors='none', alpha=0.85)
    plt.colorbar(sc, ax=ax, label='Combined score')
    top_df = df.iloc[:TOP_N]
    ax.scatter(top_df['lddt_norm'], top_df['epitope_struct'],
               s=100, facecolors='none', edgecolors='black', lw=1.5, label=f'Top {TOP_N}')
    ax.set_xlabel('lDDT (norm)'); ax.set_ylabel('Epitope 3Di conservation')
    ax.set_title('lDDT vs epitope structure'); ax.legend(fontsize=8)
    ax = axes[2]
    cons_arr = []
    for pos in range(n_cols_3di):
        col_chars = [s[pos] for s in di_seqs_list if pos < len(s)]
        non_gap   = [c for c in col_chars if c != '-']
        cons_arr.append(Counter(non_gap).most_common(1)[0][1] / len(non_gap) if non_gap else 0.0)
    ax.imshow(np.array(cons_arr).reshape(1, -1), cmap='RdYlGn', aspect='auto', vmin=0.5, vmax=1.0)
    for start, end, label, *_ in EPITOPE_RANGES:
        ax.axvspan(start-1.5, end-0.5, alpha=0.35,
                   color=epi_colors_map.get(label, 'gray'), label=label)
    ax.set_xlabel('Residue (MSA col)'); ax.set_yticks([])
    ax.set_title('Per-position 3Di conservation'); ax.legend(fontsize=7, loc='lower right')
    plt.suptitle(f'{PROTEIN_NAME} — FoldMason', fontsize=13, y=1.02)
    plt.tight_layout()
    fm_png = str(out / f'{PROTEIN_NAME}_foldmason_ranking.png')
    plt.savefig(fm_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'\n  Saved {os.path.basename(fm_png)}')

    # ── 10. Copy top conformers ───────────────────────────────────────────────
    print('\n--- Step 10: Copy top conformers ---')
    binder_dir = str(out / f'{PROTEIN_NAME}_binder_ready')
    os.makedirs(binder_dir, exist_ok=True)
    csv_path = str(out / f'{PROTEIN_NAME}_conformer_ranking.csv')
    df.to_csv(csv_path)
    copied = []
    for rank_idx, row in df.iloc[:TOP_N].iterrows():
        stem = os.path.splitext(row['name'])[0]
        cands = (glob.glob(os.path.join(filtered_dir, f'{stem}.pdb')) +
                 glob.glob(os.path.join(filtered_dir, f'*{stem}*.pdb')))
        if not cands:
            print(f'  WARNING: PDB not found for {row["name"]}')
            continue
        src = cands[0]
        dst = os.path.join(binder_dir, f'rank{rank_idx:02d}_{os.path.basename(src)}')
        shutil.copy2(src, dst)
        copied.append(dst)
        print(f'  Rank {rank_idx}: {row["name"]}  '
              f'combined={row["combined"]}  lDDT={row["lddt_raw"]}  '
              f'epitope={row["epitope_struct"]}')

    # ── 11. Zip ───────────────────────────────────────────────────────────────
    zip_path = str(out / f'{PROTEIN_NAME}_phase1_results.zip')
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for fn in sorted(os.listdir(binder_dir)):
            zf.write(os.path.join(binder_dir, fn), os.path.join('binder_ready', fn))
        for path in [msa_aa_path, msa_3di_path]:
            if os.path.exists(path):
                zf.write(path, os.path.basename(path))
        if os.path.exists(html_path):
            zf.write(html_path, os.path.basename(html_path))
        zf.write(csv_path, os.path.basename(csv_path))
        for png in [ensemble_png, fm_png]:
            if os.path.exists(png):
                zf.write(png, os.path.basename(png))

    print(f'\n=== DONE ===')
    print(f'Zip  : {zip_path}')
    print(f'CSV  : {csv_path}')
    print(f'Ready: {len(copied)} conformers in {binder_dir}/')


if __name__ == '__main__':
    main()
