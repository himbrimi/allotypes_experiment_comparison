# ── 05_glycan_averages.R ────────────────────────────────────────────────────
# Calculate per-antibody mean normalised areas for individual glycans,
# separately for EXP02 and EXP03. Mirrors the derived-trait averaging in
# 02-derived_traits_stats.R.
# Output: output/tables/06-EXP02_glycan_averages.csv
#         output/tables/06-EXP03_glycan_averages.csv
# ────────────────────────────────────────────────────────────────────────────

library(data.table)
library(tidyverse)

# ── Load renormalised data ───────────────────────────────────────────────────
load("./data/processed/00-X-exp02_without_stands.RData")  # → X (or adjust name)
exp02 <- X
load("./data/processed/00-X-exp03_without_stands.RData")
exp03 <- X

# ── Helper function ──────────────────────────────────────────────────────────
glycan_averages <- function(df, out_path) {
  aver <- df %>%
    group_by(type, glycan) %>%
    summarise(mean_narea = mean(narea), .groups = "drop")
  
  aver_w <- aver %>%
    spread(type, mean_narea)
  
  fwrite(aver_w, file = out_path)
  invisible(aver_w)
}

# ── Run ──────────────────────────────────────────────────────────────────────
glycan_averages(exp02, "output/tables/06-EXP02_glycan_averages.csv")
glycan_averages(exp03, "output/tables/06-EXP03_glycan_averages.csv")

message("Done. Glycan averages written to output/tables/")