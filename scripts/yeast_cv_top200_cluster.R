#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
})

input_top200 <- "results/yeast_stress_cv_top200.tsv"
input_expr <- "data/gasch2000.txt"
out_heatmap <- "results/cv_top200_cluster_heatmap.pdf"
out_clusters <- "results/cluster_assignment.tsv"

if (!file.exists(input_expr)) {
  stop(sprintf("Missing required input: %s", input_expr))
}

read_top200_ids <- function(path) {
  top_df <- read_tsv(path, show_col_types = FALSE)
  if ("gene_id" %in% names(top_df)) {
    ids <- top_df$gene_id
  } else if ("UID" %in% names(top_df)) {
    ids <- top_df$UID
  } else {
    ids <- top_df[[1]]
  }
  ids <- unique(as.character(ids[!is.na(ids) & nzchar(ids)]))
  ids
}

sanitize_label <- function(x) {
  x <- gsub("\\.\\.\\.[0-9]+$", "", x)
  x <- gsub("\"", "", x, fixed = TRUE)
  x <- gsub("030inutes", "030 minutes", x, fixed = TRUE)
  x <- gsub("\\s+", " ", trimws(x))
  x
}

extract_time_min <- function(label) {
  m <- regmatches(
    tolower(label),
    regexpr("([0-9]{1,3})\\s*(min|minutes)", tolower(label), perl = TRUE)
  )
  if (length(m) == 0 || identical(m, character(0)) || m == "") {
    return(NA_integer_)
  }
  as.integer(gsub("[^0-9]", "", m))
}

extract_temp_pair <- function(label) {
  # Patterns like "29C to 33C" or "17 to 37"
  if (grepl("[0-9]+\\s*C.*to\\s*[0-9]+\\s*C", label, ignore.case = TRUE, perl = TRUE)) {
    parts <- regmatches(
      label,
      regexec("([0-9]+)\\s*C.*to\\s*([0-9]+)\\s*C", label, ignore.case = TRUE, perl = TRUE)
    )[[1]]
    return(c(as.integer(parts[2]), as.integer(parts[3])))
  }
  if (grepl("[0-9]+\\s*C\\s*to\\s*[0-9]+\\s*C", label, ignore.case = TRUE, perl = TRUE)) {
    parts <- regmatches(
      label,
      regexec("([0-9]+)\\s*C\\s*to\\s*([0-9]+)\\s*C", label, ignore.case = TRUE, perl = TRUE)
    )[[1]]
    return(c(as.integer(parts[2]), as.integer(parts[3])))
  }
  if (grepl("[0-9]+\\s*to\\s*[0-9]+", label, ignore.case = TRUE, perl = TRUE)) {
    parts <- regmatches(
      label,
      regexec("([0-9]+)\\s*to\\s*([0-9]+)", label, ignore.case = TRUE, perl = TRUE)
    )[[1]]
    return(c(as.integer(parts[2]), as.integer(parts[3])))
  }
  c(NA_integer_, NA_integer_)
}

parse_replicate <- function(label) {
  m <- regmatches(label, regexpr("hs-[0-9]+", label, ignore.case = TRUE, perl = TRUE))
  if (length(m) == 0 || identical(m, character(0)) || m == "") {
    return(NA_character_)
  }
  tolower(m)
}

parse_condition_meta <- function(label) {
  clean <- sanitize_label(label)
  lc <- tolower(clean)
  temps <- extract_temp_pair(clean)
  time_min <- extract_time_min(clean)
  replicate <- parse_replicate(clean)

  condition_type <- NA_character_
  start_temp_C <- temps[1]
  end_temp_C <- temps[2]

  if (grepl("^37c to 25c shock", lc)) {
    condition_type <- "downshift"
  } else if (grepl("^heat shock\\s+[0-9]+\\s*to\\s*[0-9]+", lc)) {
    condition_type <- "upshift"
  } else if (grepl("^heat shock", lc)) {
    condition_type <- "heat_shock"
    # Classic heat-shock time-course without explicit start temp in label.
    if (is.na(end_temp_C)) {
      end_temp_C <- 37L
    }
  } else if (grepl("sorbitol", lc)) {
    condition_type <- "osmotic_plus_shift"
  } else if (grepl("29c to 33c|33c vs\\. 30c", lc)) {
    condition_type <- "mild_shift"
    if (grepl("33c vs\\. 30c", lc) && is.na(start_temp_C) && is.na(end_temp_C)) {
      start_temp_C <- 30L
      end_temp_C <- 33L
    }
  }

  if (is.na(condition_type)) {
    condition_type <- "mild_shift"
  }

  rep_tag <- ""
  if (!is.na(replicate)) {
    rep_num <- sub("^hs-", "", replicate)
    rep_tag <- sprintf(" (rep%s)", rep_num)
  }

  display_label <- clean
  if (condition_type == "heat_shock") {
    display_label <- sprintf("HS %sC %sm%s", end_temp_C, time_min, rep_tag)
  } else if (condition_type == "downshift") {
    display_label <- sprintf("Downshift %s->%s %sm", start_temp_C, end_temp_C, time_min)
  } else if (condition_type == "upshift") {
    display_label <- sprintf("Upshift %s->%s %sm", start_temp_C, end_temp_C, time_min)
  } else if (condition_type == "mild_shift") {
    display_label <- sprintf("%s->%s %sm", start_temp_C, end_temp_C, time_min)
  } else if (condition_type == "osmotic_plus_shift") {
    display_label <- sprintf("%s->%s +1M sorbitol %sm", start_temp_C, end_temp_C, time_min)
  }

  tibble(
    raw_label = label,
    cleaned_label = clean,
    condition_type = condition_type,
    start_temp_C = start_temp_C,
    end_temp_C = end_temp_C,
    time_min = time_min,
    replicate = replicate,
    display_label = display_label
  )
}

raw <- read_tsv(
  input_expr,
  col_types = cols(.default = col_double(), UID = col_character(), NAME = col_character(), GWEIGHT = col_double()),
  na = c("", "NA"),
  show_col_types = FALSE
)

if (!("UID" %in% names(raw))) {
  stop("Expected UID column in data/gasch2000.txt")
}

if (file.exists(input_top200)) {
  gene_ids <- read_top200_ids(input_top200)
  message(sprintf("Using existing top-200 list: %s", input_top200))
} else {
  message(sprintf("Top-200 file not found; regenerating: %s", input_top200))
  expr_for_cv <- raw %>% select(-any_of(c("NAME", "GWEIGHT")))
  n_cond <- ncol(expr_for_cv) - 1
  eps <- 0.1
  cv_df <- expr_for_cv %>%
    rowwise() %>%
    mutate(
      n_obs = sum(!is.na(c_across(-UID))),
      m = mean(c_across(-UID), na.rm = TRUE),
      s = sd(c_across(-UID), na.rm = TRUE),
      cv = if_else(abs(m) >= eps & n_obs >= n_cond * 0.5, s / abs(m), NA_real_)
    ) %>%
    ungroup() %>%
    filter(!is.na(cv)) %>%
    arrange(desc(cv)) %>%
    select(gene_id = UID, cv)

  top200_df <- cv_df %>% slice_head(n = 200)
  dir.create(dirname(input_top200), showWarnings = FALSE, recursive = TRUE)
  write_tsv(top200_df, input_top200)
  gene_ids <- top200_df$gene_id
  message(sprintf("Wrote regenerated top-200 list: %s", input_top200))
}

gene_ids <- unique(gene_ids)
message(sprintf("Requested genes: %d", length(gene_ids)))

expr <- raw %>%
  filter(UID %in% gene_ids) %>%
  select(-any_of(c("NAME", "GWEIGHT")))

message(sprintf("Matched genes in gasch2000: %d", nrow(expr)))
if (nrow(expr) == 0) {
  stop("No genes matched between top-200 list and gasch2000 UID column.")
}

expr_mat <- expr %>%
  column_to_rownames("UID") %>%
  as.matrix()

mode(expr_mat) <- "numeric"

message(sprintf("Matrix dims before 30-col subset: %d genes x %d conditions", nrow(expr_mat), ncol(expr_mat)))
if (ncol(expr_mat) < 30) {
  stop(sprintf("Need at least 30 condition columns, found %d", ncol(expr_mat)))
}

expr_30 <- expr_mat[, seq_len(30), drop = FALSE]
message(sprintf("Matrix dims after 30-col subset: %d genes x %d conditions", nrow(expr_30), ncol(expr_30)))

# Build explicit column metadata and standardized display labels (keep existing column order).
col_meta <- bind_rows(lapply(colnames(expr_30), parse_condition_meta)) %>%
  mutate(col_index = row_number())

# Avoid ambiguous duplicate display labels by adding explicit run index.
dup_idx <- ave(seq_len(nrow(col_meta)), col_meta$display_label, FUN = seq_along)
dup_n <- ave(seq_len(nrow(col_meta)), col_meta$display_label, FUN = length)
col_meta <- col_meta %>%
  mutate(
    display_label = if_else(
      dup_n > 1,
      sprintf("%s [run %d/%d]", display_label, dup_idx, dup_n),
      display_label
    )
  )

# Apply cleaned display labels to plotting matrix without reordering columns.
colnames(expr_30) <- col_meta$display_label

# Row-wise z-score with NA-aware statistics.
# Keep informative values, and only impute missing cells to 0 (row mean after z-scoring).
row_means <- rowMeans(expr_30, na.rm = TRUE)
row_sds <- apply(expr_30, 1, sd, na.rm = TRUE)
row_z <- sweep(expr_30, 1, row_means, "-")
row_z <- sweep(row_z, 1, row_sds, "/")

# Rows with undefined/zero SD cannot be scaled meaningfully; set whole row to 0.
bad_sd <- which(!is.finite(row_sds) | row_sds == 0)
if (length(bad_sd) > 0) {
  row_z[bad_sd, ] <- 0
  message(sprintf("Rows with undefined/zero SD set to 0: %d", length(bad_sd)))
}

# For partially missing rows, fill NA/Inf cells with 0 after z-scoring (neutral color).
bad_cells <- !is.finite(row_z)
if (any(bad_cells)) {
  row_z[bad_cells] <- 0
  message(sprintf("Non-finite z-score cells imputed to 0: %d", sum(bad_cells)))
}

row_dist <- dist(row_z, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")
clusters <- cutree(row_hclust, k = 4)

cluster_df <- tibble(
  gene_id = names(clusters),
  cluster = as.integer(clusters)
) %>% arrange(cluster, gene_id)

write_tsv(cluster_df, out_clusters)
message(sprintf("Saved cluster assignments: %s", out_clusters))

annotation_row <- data.frame(cluster = factor(clusters, levels = 1:4))
rownames(annotation_row) <- names(clusters)
annotation_colors <- list(cluster = setNames(brewer.pal(4, "Set2"), as.character(1:4)))

pdf(out_heatmap, width = 8, height = 12)
ph <- pheatmap(
  mat = row_z,
  cluster_rows = row_hclust,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  annotation_names_row = FALSE,
  show_rownames = FALSE,
  fontsize_col = 6,
  angle_col = 45,
  main = "Top 200 CV Yeast Genes: Row-wise Z-score Clustering",
  silent = TRUE
)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
grid::grid.text(
  "Condition / timepoint",
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.015, "npc"),
  gp = grid::gpar(fontsize = 10)
)
dev.off()
message(sprintf("Saved clustered heatmap: %s", out_heatmap))

message("Final x-axis labels in order:")
for (i in seq_len(nrow(col_meta))) {
  message(sprintf("%02d\t%s", i, col_meta$display_label[i]))
}

message("Done.")
