"""
Parse IUPred2A output to recommend IDP trim boundaries for Starling input.

Finds contiguous disordered regions, flags semi-structured pockets
(local IUPred minima), and suggests Starling-ready trim coordinates.

Usage:
    python parse_iupred.py --input iupred2a_output.txt --threshold 0.5
    python parse_iupred.py --input iupred2a_output.txt --plot  # requires matplotlib
"""

import argparse
import csv
from dataclasses import dataclass


@dataclass
class IUPredResult:
    pos: int
    aa: str
    iupred: float
    anchor: float


def parse_iupred_file(path: str) -> list[IUPredResult]:
    results = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                results.append(IUPredResult(
                    pos=int(parts[0]),
                    aa=parts[1],
                    iupred=float(parts[2]),
                    anchor=float(parts[3]),
                ))
            except ValueError:
                continue
    return results


def find_disordered_regions(
    results: list[IUPredResult],
    threshold: float = 0.5,
    min_length: int = 10,
) -> list[tuple[int, int, float]]:
    """Return (start, end, mean_score) for contiguous disordered stretches."""
    regions = []
    in_region = False
    start = 0
    for r in results:
        if r.iupred >= threshold and not in_region:
            in_region = True
            start = r.pos
        elif r.iupred < threshold and in_region:
            in_region = False
            end = r.pos - 1
            if end - start + 1 >= min_length:
                span = results[start - 1:end]
                mean = sum(x.iupred for x in span) / len(span)
                regions.append((start, end, mean))
    if in_region:
        end = results[-1].pos
        span = results[start - 1:]
        mean = sum(x.iupred for x in span) / len(span)
        if len(span) >= min_length:
            regions.append((start, end, mean))
    return regions


def find_structured_pockets(
    results: list[IUPredResult],
    window: int = 10,
    top_n: int = 10,
) -> list[tuple[int, int, float]]:
    """
    Find local minima in IUPred score — regions that are relatively more
    structured within an otherwise disordered sequence.
    These are prime Starling trim targets (binding-competent conformations
    are more likely here).
    """
    scores = [r.iupred for r in results]
    n = len(scores)
    pockets = []

    i = window
    while i < n - window:
        local_min = min(scores[i - window:i + window + 1])
        if scores[i] == local_min and scores[i] < 0.7:
            # Find extent of this pocket
            start = i
            end = i
            while start > 0 and scores[start - 1] < 0.75:
                start -= 1
            while end < n - 1 and scores[end + 1] < 0.75:
                end += 1
            mean = sum(scores[start:end + 1]) / (end - start + 1)
            pockets.append((results[start].pos, results[end].pos, mean))
            i = end + window  # skip past this pocket
        else:
            i += 1

    # Deduplicate overlapping pockets, keep lowest mean
    merged = []
    for pocket in sorted(pockets, key=lambda x: x[2]):
        if not merged or pocket[0] > merged[-1][1] + 5:
            merged.append(pocket)
        elif pocket[2] < merged[-1][2]:
            merged[-1] = pocket

    return sorted(merged, key=lambda x: x[2])[:top_n]


def recommend_starling_trim(
    results: list[IUPredResult],
    pockets: list[tuple[int, int, float]],
    padding: int = 20,
    max_length: int = 200,
) -> tuple[int, int]:
    """
    Suggest a Starling input trim encompassing the most structured pockets
    (the biologically interesting IDP region) with padding on each side.
    """
    if not pockets:
        # Fallback: use most variable (lowest mean) region
        scores = [r.iupred for r in results]
        window = 150
        best_start = 0
        best_mean = 1.0
        for i in range(len(scores) - window):
            m = sum(scores[i:i + window]) / window
            if m < best_mean:
                best_mean = m
                best_start = i
        return results[best_start].pos, results[min(best_start + window, len(results) - 1)].pos

    first_pocket_start = pockets[0][0]
    last_pocket_end = pockets[-1][1]

    trim_start = max(1, first_pocket_start - padding)
    trim_end = min(results[-1].pos, last_pocket_end + padding)

    # Cap at max_length, centered on pockets
    if trim_end - trim_start + 1 > max_length:
        center = (first_pocket_start + last_pocket_end) // 2
        trim_start = max(1, center - max_length // 2)
        trim_end = min(results[-1].pos, trim_start + max_length - 1)

    return trim_start, trim_end


def plot_disorder(results: list[IUPredResult], trim_start: int, trim_end: int):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot")
        return

    positions = [r.pos for r in results]
    scores = [r.iupred for r in results]
    anchor = [r.anchor for r in results]

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(positions, scores, color="#2196F3", linewidth=1, label="IUPred2 (long)")
    ax.plot(positions, anchor, color="#FF9800", linewidth=0.8, alpha=0.6, label="ANCHOR2")
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, label="Disorder threshold (0.5)")
    ax.axvspan(trim_start, trim_end, alpha=0.15, color="green", label=f"Starling trim ({trim_start}–{trim_end})")
    ax.set_xlabel("Residue position")
    ax.set_ylabel("Disorder score")
    ax.set_title("IUPred2A — Disorder Profile")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    plt.savefig("iupred_disorder_profile.png", dpi=150)
    plt.show()
    print("Plot saved to iupred_disorder_profile.png")


def extract_sequence(results: list[IUPredResult], start: int, end: int) -> str:
    seq = ""
    for r in results:
        if start <= r.pos <= end:
            seq += r.aa
    return seq


def main():
    parser = argparse.ArgumentParser(description="Parse IUPred2A output → Starling trim recommendation")
    parser.add_argument("--input", required=True, help="IUPred2A output text file")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Disorder threshold (default 0.5)")
    parser.add_argument("--padding", type=int, default=20,
                        help="Residues to pad around structured pockets (default 20)")
    parser.add_argument("--max_length", type=int, default=200,
                        help="Max Starling input length (default 200 aa)")
    parser.add_argument("--plot", action="store_true", help="Generate disorder profile plot")
    parser.add_argument("--fasta_out", default=None,
                        help="Optional: write trimmed sequence as FASTA")
    args = parser.parse_args()

    results = parse_iupred_file(args.input)
    total = len(results)
    sequence = "".join(r.aa for r in results)
    n_disordered = sum(1 for r in results if r.iupred >= args.threshold)

    print(f"\nProtein length: {total} aa")
    print(f"Disordered residues (IUPred ≥ {args.threshold}): {n_disordered}/{total} ({100*n_disordered/total:.1f}%)")
    print(f"Mean IUPred score: {sum(r.iupred for r in results)/total:.3f}")

    disordered = find_disordered_regions(results, threshold=args.threshold)
    print(f"\nDisordered regions (≥ {args.threshold}, ≥ 10 aa):")
    for start, end, mean in disordered:
        print(f"  {start:4d}–{end:4d}  ({end-start+1:4d} aa)  mean={mean:.3f}")

    pockets = find_structured_pockets(results)
    print(f"\nMost structured pockets (local IUPred minima < 0.7) — likely binding-relevant:")
    for start, end, mean in pockets:
        span_seq = extract_sequence(results, start, end)
        print(f"  {start:4d}–{end:4d}  ({end-start+1:3d} aa)  mean={mean:.3f}  seq={span_seq[:20]}{'...' if len(span_seq)>20 else ''}")

    trim_start, trim_end = recommend_starling_trim(
        results, pockets,
        padding=args.padding,
        max_length=args.max_length,
    )
    trimmed_seq = extract_sequence(results, trim_start, trim_end)

    print(f"\n{'='*60}")
    print(f"STARLING TRIM RECOMMENDATION: residues {trim_start}–{trim_end} ({len(trimmed_seq)} aa)")
    print(f"  (covers structured pockets ± {args.padding} aa padding)")
    print(f"  Sequence: {trimmed_seq[:40]}{'...' if len(trimmed_seq)>40 else ''}")
    print(f"{'='*60}")

    if args.fasta_out:
        with open(args.fasta_out, "w") as f:
            f.write(f">trimmed_{trim_start}_{trim_end}\n{trimmed_seq}\n")
        print(f"FASTA written to {args.fasta_out}")

    if args.plot:
        plot_disorder(results, trim_start, trim_end)


if __name__ == "__main__":
    main()
