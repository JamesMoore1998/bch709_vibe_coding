#!/usr/bin/env python3
"""
BCH709 Homework 1: Yeast mRNA FASTA GC distribution analysis.

Interpretation hypotheses (for submission):
1) Distinct GC-content groups can reflect transcript classes with different codon usage biases,
   which in turn can relate to translational optimization and expression programs.
2) GC-content variation may also track gene-family/genomic-context differences (for example,
   local mutational pressure and selection on RNA stability), producing subpopulation-like structure.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from statistics import mean, median, pstdev
from typing import Iterable, Iterator, List, Optional, Tuple

INPUT_FASTA = "data/mrna.fa.gz"
OUT_TSV = "results/mrna_metrics.tsv"
OUT_PNG = "results/gc_content_distribution.png"


@dataclass
class Record:
    accession: str
    length: int
    gc_content: Optional[float]


def parse_fasta_gz(path: str) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    seq_parts: List[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def extract_accession(header: str) -> str:
    return header.split()[0]


def gc_metrics(seq: str) -> Tuple[int, Optional[float]]:
    cleaned = re.sub(r"\s+", "", seq).upper()
    # Keep canonical nucleotide symbols for length; drop any other characters.
    filtered = "".join(ch for ch in cleaned if ch in {"A", "C", "G", "T", "N"})
    length = len(filtered)
    if length == 0:
        return 0, None
    gc = (filtered.count("G") + filtered.count("C")) / length
    return length, gc


def ensure_results_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_metrics_tsv(path: str, records: List[Record]) -> None:
    sorted_records = sorted(
        records,
        key=lambda r: (-1 if r.gc_content is None else 0, -(r.gc_content or -1.0), r.accession),
    )
    with open(path, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["accession", "length", "gc_content"])
        for rec in sorted_records:
            gc_str = "" if rec.gc_content is None else f"{rec.gc_content:.4f}"
            writer.writerow([rec.accession, rec.length, gc_str])


def write_gc_values(path: str, gc_values: Iterable[float]) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for value in gc_values:
            out.write(f"{value:.10f}\n")


def plot_with_r(gc_values_path: str, out_png: str) -> None:
    r_code = f"""
vals <- scan('{gc_values_path}', quiet=TRUE)
png('{out_png}', width=1600, height=900, res=200)
par(mar=c(7,5,4,2)+0.1, mgp=c(4,1,0))
h <- hist(vals, breaks=40, probability=TRUE,
          col='#9ecae1', border='white',
          xlim=c(0,1), xlab='GC content',
          main='Yeast mRNA GC Content Distribution')
lines(density(vals), lwd=2, col='#08519c')
mu <- mean(vals)
med <- median(vals)
abline(v=mu, lty=2, lwd=2, col='#cb181d')
abline(v=med, lty=2, lwd=2, col='#238b45')
legend('topright', bty='n',
       legend=c(sprintf('Mean = %.4f', mu), sprintf('Median = %.4f', med)),
       lty=2, lwd=2, col=c('#cb181d','#238b45'))
caption <- sprintf('n=%d  mean=%.4f  median=%.4f  sd=%.4f', length(vals), mu, med, sd(vals))
mtext(caption, side=1, line=5.6, cex=0.9)
dev.off()
"""
    subprocess.run(["Rscript", "-e", r_code], check=True)


def main() -> int:
    if not os.path.exists(INPUT_FASTA):
        print(f"ERROR: Missing input file: {INPUT_FASTA}", file=sys.stderr)
        return 1

    ensure_results_dir("results")

    records: List[Record] = []
    total = 0
    invalid = 0

    for header, seq in parse_fasta_gz(INPUT_FASTA):
        total += 1
        accession = extract_accession(header)
        length, gc = gc_metrics(seq)
        if gc is None:
            invalid += 1
        records.append(Record(accession=accession, length=length, gc_content=gc))

    write_metrics_tsv(OUT_TSV, records)

    valid_gc = [r.gc_content for r in records if r.gc_content is not None]
    if not valid_gc:
        print("ERROR: No valid sequences for GC statistics.", file=sys.stderr)
        return 1

    gc_tmp = "results/.gc_values_tmp.txt"
    write_gc_values(gc_tmp, valid_gc)
    try:
        plot_with_r(gc_tmp, OUT_PNG)
    finally:
        if os.path.exists(gc_tmp):
            os.remove(gc_tmp)

    print(f"Total records parsed: {total}")
    print(f"Valid sequences used for GC stats: {len(valid_gc)}")
    print(f"Invalid/empty sequences: {invalid}")
    print(f"GC min: {min(valid_gc):.4f}")
    print(f"GC max: {max(valid_gc):.4f}")
    print(f"GC mean: {mean(valid_gc):.4f}")
    print(f"GC median: {median(valid_gc):.4f}")
    print(f"GC sd: {pstdev(valid_gc):.4f}")
    print(f"Wrote: {OUT_TSV}")
    print(f"Wrote: {OUT_PNG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
