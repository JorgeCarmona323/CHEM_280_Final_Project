# IUPred2A Data Files

Store raw IUPred2A output files here. These are the plain-text downloads from https://iupred2a.elte.hu/

## Files expected

| File | Protein | Notes |
|------|---------|-------|
| `tau_iupred.txt` | Tau (Big Tau 758-aa isoform) | Run with IUPred2 long |
| `fmrp_iupred.txt` | FMRP (Fragile X, 621 aa) | Run with IUPred2 long |

## File format

IUPred2A output is tab-separated with a comment header:

```
# IUPred2A output
# ...
POS AA  IUPRED2 ANCHOR2
1   M   0.1234  0.5678
2   G   0.2345  0.6789
...
```

Lines beginning with `#` are skipped by `parse_iupred.py`.

## How to generate

1. Go to https://iupred2a.elte.hu/
2. Paste protein sequence → select **IUPred2 long** + **ANCHOR2**
3. Download the text output → save here as `{protein}_iupred.txt`

## Parse and get Starling trim

```bash
# FMRP
python src/structure_prep/parse_iupred.py \
    --input data/iupred/fmrp_iupred.txt \
    --threshold 0.5 --padding 15 --max_length 200 \
    --fasta_out data/constructs/target_protein/fmrp_idr_trimmed.fasta \
    --plot

# Tau
python src/structure_prep/parse_iupred.py \
    --input data/iupred/tau_iupred.txt \
    --threshold 0.5 --padding 20 --max_length 200 \
    --fasta_out data/constructs/tau/tau_mtbr_trimmed.fasta \
    --plot
```
