#!/usr/bin/env Rscript
# Plot top-10 genes by CV from the Gasch 2000 dataset as a publication-quality heatmap.
# Usage:
#   Rscript plot_top10_heatmap.R /path/to/gasch2000.txt out.png

suppressPackageStartupMessages({
  required <- c("readr","dplyr","tidyr","ggplot2","RColorBrewer")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) stop("Install required packages first: ", paste(missing, collapse=", "))
})

args <- commandArgs(trailingOnly = TRUE)
infile <- if(length(args) >= 1) args[1] else "gasch2000.txt"
outfile <- if(length(args) >= 2) args[2] else "gasch_top10_heatmap.png"

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

if(!file.exists(infile)) stop("Input file not found: ", infile)

# Read everything as character to detect numeric columns robustly
df_raw <- read_tsv(infile, col_types = cols(.default = "c"))

gene_col <- names(df_raw)[1]

# Determine which columns look numeric (fraction of non-missing values parseable as numeric)
is_numeric_frac <- function(x){
  x <- trimws(x)
  x[x==""] <- NA
  n_non_na <- sum(!is.na(x))
  if(n_non_na == 0) return(0)
  n_num <- sum(!is.na(suppressWarnings(as.numeric(x))))
  n_num / n_non_na
}

frac_numeric <- sapply(df_raw[-1], is_numeric_frac)
# choose columns with >= 60% numeric entries
num_cols <- names(frac_numeric)[frac_numeric >= 0.6]
if(length(num_cols) < 2) stop("Could not detect enough numeric condition columns automatically.")

# Build numeric expression table: keep gene identifier plus detected numeric columns
expr <- df_raw %>%
  select(all_of(c(gene_col, num_cols))) %>%
  rename(Gene = !!sym(gene_col)) %>%
  mutate(across(all_of(num_cols), ~as.numeric(.)))

# Quality checks: count available values per row
expr$non_na_count <- apply(expr[num_cols], 1, function(r) sum(!is.na(r)))
min_non_na_required <- ceiling(length(num_cols) * 0.5)

# Compute mean, sd, CV (sd / |mean|), handling near-zero means with small epsilon
eps <- 1e-8
expr <- expr %>%
  rowwise() %>%
  mutate(
    mean_val = mean(c_across(all_of(num_cols)), na.rm = TRUE),
    sd_val   = sd(c_across(all_of(num_cols)), na.rm = TRUE),
    cv = ifelse(abs(mean_val) < 0.1, NA_real_, sd_val / (abs(mean_val) + eps))
  ) %>%
  ungroup()

# Filter out rows with too many missing values or undefined CV
expr_filt <- expr %>% filter(non_na_count >= min_non_na_required & !is.na(cv))
if(nrow(expr_filt) < 10) stop("Not enough genes passed filters to select top 10.")

# Select top 10 by CV
top_genes <- expr_filt %>% arrange(desc(cv)) %>% slice_head(n = 10) %>% pull(Gene)

# Subset expression matrix to top genes and reshape for ggplot
top_mat <- expr %>% filter(Gene %in% top_genes) %>% select(Gene, all_of(num_cols))
long <- top_mat %>% pivot_longer(-Gene, names_to = "Condition", values_to = "Value")

# Ensure ordering: genes on x-axis in descending CV order; conditions in file order
gene_levels <- top_genes
condition_levels <- num_cols
long$Gene <- factor(long$Gene, levels = gene_levels)
long$Condition <- factor(long$Condition, levels = condition_levels)

# Optional font handling: try showtext if available (ensures embedding on macOS)
if(requireNamespace("showtext", quietly=TRUE) && requireNamespace("sysfonts", quietly=TRUE)){
  library(showtext)
  library(sysfonts)
  # Try to register Times New Roman from common macOS path, otherwise rely on system name
  ttf_path <- "/Library/Fonts/Times New Roman.ttf"
  if(file.exists(ttf_path)) sysfonts::font_add("Times_New_Roman", ttf_path)
  showtext::showtext_auto()
  base_family <- if(file.exists(ttf_path)) "Times_New_Roman" else "Times New Roman"
} else {
  base_family <- "Times New Roman"
}

# Create heatmap
p <- ggplot(long, aes(x = Gene, y = Condition, fill = Value)) +
  geom_tile(color = "black", size = 0.25) +
  scale_fill_distiller(palette = "PuGn", na.value = "grey95") +
  labs(x = "Gene", y = "Condition", fill = "Expression (log)", title = "Top 10 genes by CV") +
  theme_bw(base_size = 11, base_family = base_family) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = c(1,1),
    legend.justification = c(1,1),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save outputs
cat("Top 10 genes (by CV):\n")
print(top_genes)
write.csv(data.frame(Gene = top_genes), file = "top10_genes.csv", row.names = FALSE)

ggsave(outfile, p, width = 5, height = 5, units = "in", dpi = 300, bg = "white")
cat("Saved heatmap to:", outfile, "\n")

invisible(list(top_genes = top_genes, plot = p))
