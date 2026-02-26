# ============================================================
# gasch_heatmap.R
# Heat map: top 10 genes by CV from Gasch 2000 dataset
# Output: results/gasch_top10_heatmap.png (5x5 in, 300 DPI)
# ============================================================

# ── 1. Packages ─────────────────────────────────────────────
required_pkgs <- c("readr", "dplyr", "tidyr", "ggplot2",
                   "RColorBrewer", "showtext", "sysfonts")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(showtext)
library(sysfonts)

# ── 2. Font setup (Times New Roman via showtext) ─────────────
font_add("Times New Roman",
         regular      = "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
         bold         = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
         italic       = "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf",
         bolditalic   = "/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf")
showtext_auto()

# ── 3. Load data ─────────────────────────────────────────────
# Columns: UID | NAME | GWEIGHT | <conditions...>
# Skip NAME and GWEIGHT; keep UID + all numeric condition columns.
raw <- read_tsv(
  "../data/gasch2000.txt",
  col_types = cols(.default = col_double(), UID = col_character(),
                   NAME = col_character(), GWEIGHT = col_double()),
  na = c("", "NA")
)

# Drop NAME and GWEIGHT; rename UID → gene
expr <- raw %>%
  select(-NAME, -GWEIGHT) %>%
  rename(gene = UID)

# ── 4. Values are already log2 ratios — confirm range ────────
vals <- unlist(select(expr, -gene), use.names = FALSE)
message(sprintf("Value range (ignoring NA): [%.2f, %.2f]",
                min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))

# ── 5. Compute CV per gene and select top 10 ─────────────────
# CV = sd / |mean|; exclude genes where |mean| < eps or
# too many NAs (require at least 50% non-missing).
n_cond <- ncol(expr) - 1   # number of condition columns
eps    <- 0.1              # minimum |mean| to avoid CV blowup

cv_df <- expr %>%
  rowwise() %>%
  mutate(
    n_obs  = sum(!is.na(c_across(-gene))),
    m      = mean(c_across(-gene), na.rm = TRUE),
    s      = sd(c_across(-gene),   na.rm = TRUE),
    cv     = if_else(abs(m) >= eps & n_obs >= n_cond * 0.5,
                     s / abs(m), NA_real_)
  ) %>%
  ungroup()

top10_genes <- cv_df %>%
  filter(!is.na(cv)) %>%
  slice_max(cv, n = 10) %>%
  pull(gene)

message("Top 10 genes by CV:")
message(paste(top10_genes, collapse = ", "))

# ── 6. Subset and reshape: rows = conditions, cols = genes ───
mat <- expr %>%
  filter(gene %in% top10_genes) %>%
  column_to_rownames("gene")        # genes × conditions

# Transpose so conditions become rows (Y-axis)
mat_t <- as.data.frame(t(mat))      # conditions × genes

# Select top 30 conditions by variance across the top-10 genes
cond_var <- apply(mat_t, 1, var, na.rm = TRUE)
top30_conds <- names(sort(cond_var, decreasing = TRUE))[1:30]
mat_t <- mat_t[top30_conds, , drop = FALSE]
message(sprintf("Conditions in plot: %d", nrow(mat_t)))

# Convert to long format for ggplot2
mat_t$condition <- rownames(mat_t)

long <- mat_t %>%
  pivot_longer(cols = -condition,
               names_to  = "gene",
               values_to = "expression")

# ── 7. Order genes by CV (highest first) on x-axis ───────────
gene_order <- cv_df %>%
  filter(gene %in% top10_genes) %>%
  arrange(desc(cv)) %>%
  pull(gene)

long <- long %>%
  mutate(gene = factor(gene, levels = gene_order))

# ── 8. Build the heat map ─────────────────────────────────────
p <- ggplot(long, aes(x = gene, y = condition, fill = expression)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_distiller(
    palette   = "PRGn",
    direction = 1,
    name      = "Gene expression\nLog scale",
    na.value  = "grey80"
  ) +
  labs(
    title = "Gene Top 10",
    x     = "Gene",
    y     = "Conditions"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    # Tile area
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    # Background
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white"),
    # Axis text
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1,
                                 size = 7, family = "Times New Roman"),
    axis.text.y  = element_text(size = 7, family = "Times New Roman"),
    axis.title   = element_text(size = 10, family = "Times New Roman"),
    plot.title   = element_text(size = 12, family = "Times New Roman",
                                hjust = 0.5),
    # Legend: top-right inside plot
    legend.position        = c(1, 1),
    legend.justification   = c(1, 1),
    legend.background      = element_rect(fill = alpha("white", 0.7),
                                           color = "grey60", linewidth = 0.3),
    legend.title           = element_text(size = 8,  family = "Times New Roman"),
    legend.text            = element_text(size = 7,  family = "Times New Roman"),
    legend.key.height      = unit(0.4, "cm"),
    legend.key.width       = unit(0.3, "cm")
  )

# ── 9. Export ─────────────────────────────────────────────────
dir.create("../results", showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = "../results/gasch_top10_heatmap.png",
  plot     = p,
  width    = 5,
  height   = 5,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

message("Saved: results/gasch_top10_heatmap.png")
