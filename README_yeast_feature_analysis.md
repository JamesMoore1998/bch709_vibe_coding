# Yeast Feature Analysis Environment (micromamba)

This README documents how to recreate the Python environment and rerun the yeast GFF3 feature distribution analysis in this repo.

## What this analysis does

The script `scripts/analyze_yeast_features.py` analyzes a Saccharomyces cerevisiae GFF3 file and answers:

- How are genomic features (`gene`, `exon`, `tRNA`, `snoRNA`) distributed across chromosomes?
- Does feature density correlate with chromosome size?

Key counting rules:

- Exons are counted as unique intervals per chromosome using `(start, end, strand)` to prevent overcounting.
- Chromosome lengths come only from `data/chrom.sizes`.
- GFF3 seqids not found in `chrom.sizes` are excluded and logged to `results/dropped_seqids.txt`.

## Files used

- GFF3 input: `data/saccharomyces_cerevisiae.gff.gz`
- Chrom sizes: `data/chrom.sizes`
- Script: `scripts/analyze_yeast_features.py`

Note: The script reads `.gz` files directly. You do not need to unzip the GFF3.

## Create the micromamba environment

If the environment does not exist yet:

```bash
micromamba create -n bch709_vibe_coding -c conda-forge python=3.11 pandas numpy matplotlib seaborn biopython tqdm scipy
```

Why `numpy` and `scipy`?

- `numpy` is used by the script.
- `scipy` is used for Pearson/Spearman p-values when available (the script can partially fall back if needed).

## Activate the environment

```bash
micromamba activate bch709_vibe_coding
```

If activation is not configured in your shell yet (macOS + zsh), run once:

```bash
micromamba shell init -s zsh -r ~/micromamba
```

Then restart your terminal.

## Run the analysis

From the repo root:

```bash
python3 scripts/analyze_yeast_features.py \
  --gff3 data/saccharomyces_cerevisiae.gff.gz \
  --chrom-sizes data/chrom.sizes \
  --outdir results \
  --prefix saccharomyces
```

Alternative (no activation required):

```bash
micromamba run -n bch709_vibe_coding python3 scripts/analyze_yeast_features.py \
  --gff3 data/saccharomyces_cerevisiae.gff.gz \
  --chrom-sizes data/chrom.sizes \
  --outdir results \
  --prefix saccharomyces
```

## Outputs (written to `results/`)

- `saccharomyces_feature_distribution_by_chromosome.csv`
- `saccharomyces_feature_density_correlations.csv`
- `dropped_seqids.txt`
- `saccharomyces_counts_by_chromosome.png`
- `saccharomyces_densities_by_chromosome.png`
- `saccharomyces_density_vs_chromosome_size_scatter.png`

## Optional: unzip the GFF3 on macOS

You usually do not need this, but if you want an uncompressed copy:

```bash
gunzip -k data/saccharomyces_cerevisiae.gff.gz
```

This keeps the original `.gz` and creates:

- `data/saccharomyces_cerevisiae.gff`

Then run the script with `--gff3 data/saccharomyces_cerevisiae.gff`.

## Common micromamba issue (libmamba error)

If you see:

- `No prefix found at: .../envs/bch709_vibe_coding`
- `Environment must first be created ...`

That means the environment name does not exist yet (even if your shell prompt shows `bch709_vibe_coding git:(master)`, which is usually just the repo folder and Git branch).

Fix:

```bash
micromamba create -n bch709_vibe_coding -c conda-forge python=3.11 pandas numpy matplotlib seaborn biopython tqdm scipy
```

## Quick checks (later)

List environments:

```bash
micromamba env list
```

Check Python path inside the environment:

```bash
micromamba run -n bch709_vibe_coding python3 -c "import sys; print(sys.executable)"
```

Check required packages:

```bash
micromamba run -n bch709_vibe_coding python3 -c "import pandas, numpy, matplotlib, seaborn, Bio, tqdm, scipy; print('ok')"
```

