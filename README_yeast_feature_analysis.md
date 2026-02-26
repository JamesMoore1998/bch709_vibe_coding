# Yeast Feature Analysis Environment (micromamba)

This README documents how to recreate the Python environment and rerun the yeast GFF3 feature distribution analysis in this repo.

## What this analysis does

The script `scripts/analyze_yeast_features.py` analyzes a Saccharomyces cerevisiae GFF3 file and answers:

- How are genomic features (`gene`, `exon`, `tRNA`, `snoRNA`) distributed across chromosomes?
- Does feature density correlate with chromosome size?

Key counting rules:

- Exons are counted as unique intervals per chromosome using `(start, end, strand)` to prevent overcounting.
- `noncoding_exon` rows are included in the exon count (treated as exon-like features).
- Chromosome lengths come only from `data/chrom.sizes`.
- Mitochondrial seqid aliases (for example `chrmt` in GFF vs `chrM` in `chrom.sizes`) are normalized when possible.
- GFF3 seqids not found in `chrom.sizes` are excluded and logged to `results/dropped_seqids.txt`.

## Files used

- GFF3 input: `data/saccharomyces_cerevisiae.gff.gz`
- Chrom sizes: `data/chrom.sizes`
- Script: `scripts/analyze_yeast_features.py`

Note: The script reads `.gz` files directly. You do not need to unzip the GFF3.

## Create the micromamba environment

If the environment does not exist yet:

```bash
micromamba create -n bch709_vibe_coding -c conda-forge python=3.11 pandas numpy matplotlib seaborn biopython tqdm scipy gawk
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
micromamba create -n bch709_vibe_coding -c conda-forge python=3.11 pandas numpy matplotlib seaborn biopython tqdm scipy gawk
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

Check GNU awk:

```bash
micromamba run -n bch709_vibe_coding gawk --version | head -1
```

## Synthetic test data + GNU one-liner QC

Small fixture files are included for QC regression checks:

- `tests/fixtures/yeast_qc/chrom.sizes`
- `tests/fixtures/yeast_qc/test.gff3`
- `tests/fixtures/yeast_qc/expected_qc_feature_distribution_by_chromosome.csv`
- `tests/fixtures/yeast_qc/expected_dropped_seqids.txt`

This fixture intentionally includes:

- duplicate exon rows on the same chromosome/strand (should count once)
- same exon coordinates on opposite strands (should count twice)
- a `chrUn` feature not present in `chrom.sizes` (should be dropped)

Run the analysis on the fixture:

```bash
python3 scripts/analyze_yeast_features.py \
  --gff3 tests/fixtures/yeast_qc/test.gff3 \
  --chrom-sizes tests/fixtures/yeast_qc/chrom.sizes \
  --outdir tests/fixtures/yeast_qc/out \
  --prefix qc
```

GNU `awk` one-liner QC (validates counts and density columns in `qc_feature_distribution_by_chromosome.csv` against the fixture):

```bash
gawk -F'\t|,' '
BEGIN{tol=1e-6}
FILENAME==ARGV[1] && $1 !~ /^#/ && NF>=2 {len[$1]=$2+0; next}
FILENAME==ARGV[2] && $0 !~ /^#/ && NF==9 {
  c=$1; t=$3; s=$4+0; e=$5+0; st=$7;
  if(!(c in len)) next;
  if(t=="gene") g[c]++;
  else if(t=="tRNA") trna[c]++;
  else if(t=="snoRNA") sno[c]++;
  else if(t=="exon"){k=c SUBSEP s SUBSEP e SUBSEP st; if(!(k in seen)){seen[k]=1; ex[c]++}}
  next
}
FILENAME==ARGV[3] && FNR==1 {next}
FILENAME==ARGV[3] {
  c=$1; L=$2+0;
  if(!(c in len)) {print "FAIL unexpected chromosome in CSV:", c; bad=1; next}
  seen_csv[c]++;
  exp_g=(c in g)?g[c]:0; exp_ex=(c in ex)?ex[c]:0; exp_t=(c in trna)?trna[c]:0; exp_s=(c in sno)?sno[c]:0;
  if($3+0!=exp_g || $4+0!=exp_ex || $5+0!=exp_t || $6+0!=exp_s){print "FAIL count mismatch", c, "csv=" $3,$4,$5,$6, "exp=" exp_g,exp_ex,exp_t,exp_s; bad=1}
  mb=L/1000000.0;
  if(mb<=0){print "FAIL nonpositive length", c; bad=1}
  if((($7+0)-(exp_g/mb))^2>tol || (($8+0)-(exp_ex/mb))^2>tol || (($9+0)-(exp_t/mb))^2>tol || (($10+0)-(exp_s/mb))^2>tol){
    print "FAIL density mismatch", c; bad=1
  }
  next
}
END{
  for(c in len){
    if(!(c in seen_csv)){print "FAIL missing chromosome in CSV:", c; bad=1}
    if(seen_csv[c]>1){print "FAIL duplicate CSV row for chromosome:", c; bad=1}
  }
  if(!bad) print "PASS QC: counts and densities match fixture-derived expectations"
  exit bad
}' tests/fixtures/yeast_qc/chrom.sizes tests/fixtures/yeast_qc/test.gff3 tests/fixtures/yeast_qc/out/qc_feature_distribution_by_chromosome.csv
```

Optional exact-file checks after running:

```bash
diff -u tests/fixtures/yeast_qc/expected_qc_feature_distribution_by_chromosome.csv tests/fixtures/yeast_qc/out/qc_feature_distribution_by_chromosome.csv
diff -u tests/fixtures/yeast_qc/expected_dropped_seqids.txt tests/fixtures/yeast_qc/out/dropped_seqids.txt
```
