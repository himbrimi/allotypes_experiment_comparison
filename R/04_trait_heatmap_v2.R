# 05_trait_heatmap.R
# Heatmap of glycan trait means per antibody (HC × LC), per experiment.
# Rows = traits, columns = 8 antibodies, annotated by HC and LC.
# High values: warm (red/orange); low values: cold (blue).
#
# Input:  long-format data frames with columns:
#         Sample, peptide, type, LC, HC, glycan, narea, tmpid
#         One row per replicate × trait observation.
#         Adjust object names in DATA list below if yours differ.

library(tidyverse)
library(pheatmap)

## load data
load("./data/processed/02-X_dt_EXP02.Rdata") # loads traits_exp02
traits_exp02 <- X_dt

load("./data/processed/02-X_dt_EXP03.RData") # loads traits_exp03
traits_exp03 <- X_dt



# ── 0. Config ────────────────────────────────────────────────────────────────

TRAITS <- c("Agalactosylation", "Galactosylation", "Sialylation", "High Mannose", "Bisection", "Antennary fucosylation",
            "Monoantennary", "Hybrid")

TRAIT_LABELS <- setNames(TRAITS, TRAITS)

ANCHOR_COL <- "Y_NO"

anchor_callback <- function(hc, ...) {
  if (ANCHOR_COL %in% hc$labels) {                 # only the column tree has this label
    wts <- ifelse(hc$labels == ANCHOR_COL, -Inf, seq_along(hc$labels))
    d <- reorder(as.dendrogram(hc), wts, agglo.FUN = mean)
    return(as.hclust(d))
  }
  hc  
}

# Traits with known inter-experiment glycoform discrepancies — flagged with *
AFFECTED_TRAITS <- c("Sialylation", "Bisection", "Antennary fucosylation", "Monoantennary")

# Adjust these to your actual Step 2 output object names
DATA <- list(
  EXP02 = traits_exp02,
  EXP03 = traits_exp03
)

OUT_DIR <- "output/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# Clustering settings (applied to the row z-scored matrix)
CLUST_DIST   <- "euclidean"   # "euclidean" or "correlation"
CLUST_METHOD <- "ward.D2"     # e.g. "ward.D2", "complete", "average"


# ── 1. Colour palette: cold (blue) → white → warm (red) ─────────────────────

heatmap_colours <- colorRampPalette(
  c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#FDAE61", "#D6604D", "#B2182B")
)(100)

# ── 2. Annotation colours ────────────────────────────────────────────────────

# hc_colours <- c(Y = "#4DAF4A", YA = "#377EB8", YF = "#FF7F00", YI = "#E41A1C")
# lc_colours <- c(WT = "#999999", NO = "#333333")

# ── 3. Function: build matrix and plot ───────────────────────────────────────

make_heatmap <- function(df, exp_label) {
  
  # Compute means per trait per HC × LC condition from long format
  means <- df |>
    filter(glycan %in% TRAITS) |>
    group_by(HC, LC, glycan) |>
    summarise(mean_narea = mean(narea, na.rm = TRUE), .groups = "drop") |>
    mutate(antibody = paste0(HC, "_", LC)) |>
    arrange(HC, LC)
  
  # Pivot to wide: rows = traits, columns = antibodies
  mat <- means |>
    select(glycan, antibody, mean_narea) |>
    pivot_wider(names_from = antibody, values_from = mean_narea) |>
    # Enforce trait order
    mutate(glycan = factor(glycan, levels = TRAITS)) |>
    arrange(glycan) |>
    column_to_rownames("glycan") |>
    as.matrix()
  
  
  
  # Drop glycans that are all-NA or constant (can't be z-scored / clustered)
  keep <- apply(mat, 1, function(r) sum(!is.na(r)) >= 2 && sd(r, na.rm = TRUE) > 0)
  mat  <- mat[keep, , drop = FALSE]

  # Keep raw means matrix (same dims) for cell labels — formatted to 1 decimal place
  mat_raw <- mat
  mat_labels <- matrix(
    sprintf("%.1f", mat_raw),
    nrow = nrow(mat_raw), ncol = ncol(mat_raw),
    dimnames = dimnames(mat_raw)
  )

  # Z-score each row so colour encodes deviation from trait mean
  mat <- t(scale(t(mat)))
  
  # Replace any residual NA/NaN with 0 (= row mean in z-space) so clustering runs
  mat[is.na(mat)] <- 0
  mat_labels[is.na(mat_labels)] <- ""

  # Symmetric breaks centred on 0
  max_abs <- max(abs(mat), na.rm = TRUE)
  breaks  <- seq(-max_abs, max_abs, length.out = 101)
  
  # # Row labels: full names, * for affected traits
  # rownames(mat) <- ifelse(
  #   rownames(mat) %in% AFFECTED_TRAITS,
  #   paste0(rownames(mat), "*"),
  #   rownames(mat)
  # )
  rownames(mat_labels) <- rownames(mat)

  pheatmap(
    mat,
    color                    = heatmap_colours,
    breaks                   = breaks,
    scale                    = "none",          # already z-scored above
    cluster_rows             = TRUE,            # cluster glycans
    cluster_cols             = TRUE,            # cluster sample types
    clustering_distance_rows = CLUST_DIST,
    clustering_distance_cols = CLUST_DIST,
    clustering_method        = CLUST_METHOD,
    clustering_callback      = anchor_callback,
    annotation_col           = NA,              # no HC/LC colour bands
    show_colnames            = TRUE,
    show_rownames            = TRUE,
    display_numbers          = mat_labels,      # raw (non-z-scored) means in each cell
    number_color             = "black",
    fontsize                 = 10,
    fontsize_row             = 9,
    fontsize_col             = 9,
    fontsize_number          = 7,
    angle_col                = 90,
    border_color             = "white",
    cellwidth                = 42,             # slightly wider to fit the number
    cellheight               = 22,
    legend                   = TRUE,
    #main                    = paste0(exp_label, "— trait z-scores per antibody (row-scaled)"),
    filename                 = file.path(OUT_DIR, paste0("04_02_heatmap_", tolower(exp_label), ".pdf")),
    width                    = 9,
    height                   = 5.5
  )
  
  message("Saved: ", file.path(OUT_DIR, paste0("04_0_heatmap_", tolower(exp_label), ".pdf")))
}

# ── 4. Run for both experiments ───────────────────────────────────────────────

iwalk(DATA, make_heatmap)
