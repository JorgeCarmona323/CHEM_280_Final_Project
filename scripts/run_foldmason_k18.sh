#!/usr/bin/env bash
# =============================================================================
# run_foldmason_k18.sh
#
# Runs FoldMason easy-msa on the Tau K18 Starling ensemble (100 frames)
# and produces the JSON report used by Notebook 02 for:
#   - lDDT-based conformer ranking
#   - Per-position 3Di conservation → refined Starling trim boundaries
#   - Epitope structural content scoring (PHF6* / PHF6)
#
# Usage:
#   bash scripts/run_foldmason_k18.sh
#
# Outputs:
#   results/foldmason/tau_k18_raw/
#     msa_aa.fa          <- amino acid structural MSA
#     msa_3di.fa         <- 3Di structural MSA
#     msa.html           <- interactive viewer
#     foldmason.json     <- lDDT + 3Di per conformer  <- Notebook 02 input
#     msa.nw             <- Newick guide tree
#
# Runtime: ~2-5 min on CPU for 100 PDB frames
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FOLDMASON="${PROJECT_ROOT}/foldmason/bin/foldmason"

FRAMES_DIR="${PROJECT_ROOT}/data/constructs/tau/tau_k18_starling_results (1)/tau_k18_frames"
OUTPUT_DIR="${PROJECT_ROOT}/results/foldmason/tau_k18_raw"
TMP_DIR="${PROJECT_ROOT}/tmp/foldmason_tau_k18_raw"
MSA_PREFIX="${OUTPUT_DIR}/msa"

THREADS=4   # match your CPU count (4 cores available)

# ── Preflight checks ──────────────────────────────────────────────────────────

if [[ ! -x "$FOLDMASON" ]]; then
    echo "ERROR: FoldMason binary not found at: $FOLDMASON"
    echo "Run from project root: wget https://mmseqs.com/foldmason/foldmason-linux-avx2.tar.gz && tar xzf foldmason-linux-avx2.tar.gz"
    exit 1
fi

if [[ ! -d "$FRAMES_DIR" ]]; then
    echo "ERROR: Frames directory not found: $FRAMES_DIR"
    exit 1
fi

N_FRAMES=$(find "$FRAMES_DIR" -name "frame_*.pdb" | wc -l)
if [[ "$N_FRAMES" -eq 0 ]]; then
    echo "ERROR: No frame_*.pdb files found in $FRAMES_DIR"
    exit 1
fi

echo "============================================================"
echo "FoldMason — Tau K18 ensemble alignment"
echo "============================================================"
echo "Frames dir : $FRAMES_DIR"
echo "PDB count  : $N_FRAMES frames"
echo "Output dir : $OUTPUT_DIR"
echo "Threads    : $THREADS"
echo "------------------------------------------------------------"

# ── Setup ─────────────────────────────────────────────────────────────────────

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

# ── Run FoldMason easy-msa ────────────────────────────────────────────────────
# --report-mode 2 → generates foldmason.json (required for Notebook 02)

echo "Running FoldMason easy-msa ..."
"$FOLDMASON" easy-msa \
    "$FRAMES_DIR" \
    "$MSA_PREFIX" \
    "$TMP_DIR" \
    --report-mode 2 \
    --threads "$THREADS"

echo ""
echo "============================================================"
echo "Done. Output files:"
echo "============================================================"
ls -lh "$OUTPUT_DIR"/

# ── Copy JSON to notebook-accessible location ─────────────────────────────────

JSON_SRC="${OUTPUT_DIR}/msa_report.json"
JSON_DST="${OUTPUT_DIR}/foldmason.json"

# FoldMason names the JSON 'msa_report.json' — copy to standard name for nb02
if [[ -f "$JSON_SRC" ]] && [[ ! -f "$JSON_DST" ]]; then
    cp "$JSON_SRC" "$JSON_DST"
    echo "Copied: msa_report.json → foldmason.json"
elif [[ -f "$JSON_DST" ]]; then
    echo "JSON: $JSON_DST"
else
    # Try any .json in output dir
    JSON_FOUND=$(find "$OUTPUT_DIR" -name "*.json" | head -1)
    if [[ -n "$JSON_FOUND" ]]; then
        cp "$JSON_FOUND" "$JSON_DST"
        echo "Copied: $(basename $JSON_FOUND) → foldmason.json"
    else
        echo "WARNING: No JSON found in output. Check --report-mode 2 is supported."
    fi
fi

echo ""
echo "Next steps:"
echo "  1. Upload ${JSON_DST} to Colab"
echo "  2. Run Notebook 02 (02_foldmason_conformer_ranking.ipynb)"
echo "     → Plot 3 (per-position 3Di conservation) shows trim boundaries"
echo "     → Top 3-5 conformers selected for binder design"
echo "  3. Use refined trim coordinates to re-run Starling on Colab"
echo "============================================================"
