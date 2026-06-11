# 05_glycan_heatmap.R
# Heatmap of directly measured glycan means per antibody (HC × LC), per experiment.
# Rows = individual measured glycans, columns = 8 antibodies (sample types).
# BOTH rows (glycans) and columns (sample types) are clustered by similarity.
# Columns are annotated by HC and LC so the column dendrogram can be read
# against antibody identity.
# High values: warm (red/orange); low values: cold (blue).
#
# Input:  long-format data frames with columns:
#         Sample, peptide, type, LC, HC, glycan, narea, tmpid
#         One row per replicate × glycan observation.
#         Adjust object names / paths in the DATA block below if yours differ.

library(tidyverse)
library(pheatmap)

## load data ------------------------------------------------------------------
## NOTE: these should be the *directly measured glycan* objects, NOT the derived
## traits used in 04_trait_heatmap.R. Rename the files / objects below to match
## your actual Step 2 glycan output (e.g. X_glycans, X_gl, X_norm, ...).
load("./data/processed/02-X_glycans_EXP02.Rdata") # adjust if your object differs
glycans_exp02 <- X_glycans

load("./data/processed/02-X_glycans_EXP03.RData") # adjust if your object differs
glycans_exp03 <- X_glycans


# ── 0. Config ────────────────────────────────────────────────────────────────

# Derived traits to EXCLUDE: these are computed summaries, not measured glycans.
# Any row whose `glycan` matches one of these is dropped, so the heatmap shows
# only directly measured glycan species. (Safe even if your object is glycans-only.)
DERIVED_TRAITS <- c("Agalactosylation", "Galactosylation", "Sialylation",
                    "High Mannose", "Bisection", "Antennary fucosylation",
                    "Monoantennary", "Hybrid")

# Clustering settings (applied to the row z-scored matrix)
CLUST_DIST   <- "euclidean"   # "euclidean" or "correlation"
CLUST_METHOD <- "ward.D2"     # e.g. "ward.D2", "complete", "average"

DATA <- list(
  EXP02 = glycans_exp02,
  EXP03 = glycans_exp03
)

OUT_DIR <- "output/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Colour palette: cold (blue) → white → warm (red) ─────────────────────

heatmap_colours <- colorRampPalette(
  c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#FDAE61", "#D6604D", "#B2182B")
)(100)

# ── 2. Function: build matrix and plot ───────────────────────────────────────

make_glycan_heatmap <- function(df, exp_label) {

  # Mean signal per measured glycan per HC × LC condition (from long format)
  means <- df |>
    filter(!glycan %in% DERIVED_TRAITS) |>   # keep only directly measured glycans
    group_by(HC, LC, glycan) |>
    summarise(mean_narea = mean(narea, na.rm = TRUE), .groups = "drop") |>
    mutate(antibody = paste0(HC, "_", LC))

  # Wide matrix: rows = glycans, columns = antibodies (order set by clustering)
  mat <- means |>
    select(glycan, antibody, mean_narea) |>
    pivot_wider(names_from = antibody, values_from = mean_narea) |>
    column_to_rownames("glycan") |>
    as.matrix()

  # Column annotation (HC / LC) — helps interpret the column dendrogram
  annotation_col <- means |>
    distinct(antibody, HC, LC) |>
    column_to_rownames("antibody")
  annotation_col <- annotation_col[colnames(mat), , drop = FALSE]

  # Drop glycans that are all-NA or constant (can't be z-scored / clustered)
  keep <- apply(mat, 1, function(r) sum(!is.na(r)) >= 2 && sd(r, na.rm = TRUE) > 0)
  mat  <- mat[keep, , drop = FALSE]

  # Z-score each row so colour AND clustering encode pattern, not magnitude
  mat <- t(scale(t(mat)))

  # Replace any residual NA/NaN with 0 (= row mean in z-space) so clustering runs
  mat[is.na(mat)] <- 0

  # Symmetric breaks centred on 0 so the diverging palette is balanced
  max_abs <- max(abs(mat), na.rm = TRUE)
  breaks  <- seq(-max_abs, max_abs, length.out = 101)

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
    annotation_col           = annotation_col,
    annotation_names_col     = TRUE,
    show_colnames            = TRUE,
    show_rownames            = TRUE,
    fontsize                 = 10,
    fontsize_row             = 7,
    fontsize_col             = 9,
    angle_col                = 90,
    border_color             = "white",
    cellwidth                = 16,
    cellheight               = 10,
    treeheight_row           = 30,
    treeheight_col           = 25,
    legend                   = TRUE,
    main                     = paste0(exp_label, " — measured glycan z-scores per antibody (row-scaled, clustered)"),
    filename                 = file.path(OUT_DIR, paste0("05_glycan_heatmap_", tolower(exp_label), ".pdf")),
    width                    = 8,
    height                   = max(5.5, 0.16 * nrow(mat) + 2)  # grow with glycan count
  )

  message("Saved: ", file.path(OUT_DIR, paste0("05_glycan_heatmap_", tolower(exp_label), ".pdf")))
}

# ── 3. Run for both experiments ───────────────────────────────────────────────

iwalk(DATA, make_glycan_heatmap)
