# =============================================================================
# EXP02 — Analyte curation
# Rules derived from DT sheet and 'check sum spectra' sheet of EXP02_LCMS_data.xlsx
#
# The input data has one row per sample × charge state (Charge = 2 or 3).
# Each peptide × glycan combination is a column, e.g. IgGI1H4N2.
#
# Curation strategy:
#   (a) Zero-out entire glycan columns that are not integrated for ANY allotype.
#   (b) Zero-out glycan columns that are absent / not integrated for a specific
#       allotype (HC type), using the 'type' column to identify rows.
#   (c) Zero-out 3+ rows for glycans where only 2+ is valid (globally or per
#       allotype).
#
# "Zero-out" means setting the value to NA (not zero), to distinguish
# genuinely absent signal from missing/excluded data in downstream steps.
#
# HC → peptide mapping:
#   Y   → IgGI1
#   YA  → IgGIA1
#   YF  → IgGIF1
#   YI  → IgGIILE1
# =============================================================================

library(dplyr)
library(stringr)
library(readxl)

# Helper: extract HC from 'type' column (e.g. "YA-WT " → "YA")
get_hc <- function(type_col) str_extract(str_trim(type_col), "^[^-]+")

# Helper: zero-out (→ NA) a glycan for specific peptides and/or charge states
# peptides: character vector of peptide prefixes, e.g. c("IgGI1", "IgGIA1")
#           NULL = all peptides
# charges:  integer vector, e.g. c(3L)  NULL = both
# hc_types: character vector of HC types to restrict to, e.g. c("YA", "YI")
#           NULL = all types
nullify <- function(df, glycan, peptides = NULL, charges = NULL, hc_types = NULL) {

  all_peptides <- c("IgGI1", "IgGIA1", "IgGIF1", "IgGIILE1")
  if (is.null(peptides)) peptides <- all_peptides

  cols <- paste0(peptides, glycan)
  cols <- cols[cols %in% names(df)]   # guard against absent columns

  if (length(cols) == 0) return(df)

  row_mask <- rep(TRUE, nrow(df))
  if (!is.null(charges))  row_mask <- row_mask & (df$Charge %in% charges)
  if (!is.null(hc_types)) row_mask <- row_mask & (get_hc(df$type) %in% hc_types)

  df[row_mask, cols] <- NA
  df
}

# ---- Load data ---------------------------------------------------------------
# Replace with your actual path / import step
df <- read_excel("../EXP02_contamination_flagged_cleaned.xlsx")

# Ensure Charge is integer
df$Charge <- as.integer(df$Charge)

# =============================================================================
# RULE BLOCK 1: Glycans not integrated in ANY allotype → zero out entirely
# Source: DT sheet 'decision' column = blank / "did not integrate" / "-" for all
# =============================================================================

not_integrated_globally <- c(
  "H3N4",         # row 2:  all "-"
  "H3N6F3",       # row 5:  not integrating (possibly present in YI but excluded)
  "H4N3",         # row 7:  no clear decision to integrate; FA1G1 assignment unclear
  "H4N3F2",       # row 9:  "did not integrate but present"
  "H4N4",         # row 10: all "-"
  "H4N5F2S1",     # row 16: all "-"
  "H5N4",         # row 18: all "-"
  "H5N4F3",       # row 24: x / absent
  "H5N5F1S1",     # row 26: "not integrating, possibly present in Y"
  "H5N5F2",       # row 27: all "-"
  "H5N5F3",       # row 28: all "-"
  "H6N3F1G1",     # row 30 equiv: all "-"
  "H6N3F2",       # row 30: all "-"
  "H6N5F1",       # row 31: all "-"
  "H5N3F2",       # row 34: "most probably present, but would not integrate"
  "H6N3F3"        # DT row 36: "see below" — no integration decision; globally excluded
)

for (g in not_integrated_globally) {
  df <- nullify(df, glycan = g)
}

# =============================================================================
# RULE BLOCK 2: Glycans not integrated for specific allotypes
# Source: DT sheet per-allotype columns ("+" vs "-" vs notes)
# =============================================================================

# H5N4F1S2: absent in YA; low but retained for YI, Y, and YF
# DT row 21: YA = "-", YI = "low" → exclude YA only
df <- nullify(df, glycan = "H5N4F1S2", hc_types = "YA")

# H5N4F2: absent in YA and YI
# DT row 22: YA = "-", YI = "-"
df <- nullify(df, glycan = "H5N4F2", hc_types = c("YA", "YI"))

# H5N4F2S1: absent in YA and YI
# DT row 23: YA = "-", YI = "-"
df <- nullify(df, glycan = "H5N4F2S1", hc_types = c("YA", "YI"))

# H5N5F1: absent in YA and YI
# DT row 25: YA = "-", YI = "-"
df <- nullify(df, glycan = "H5N5F1", hc_types = c("YA", "YI"))

# =============================================================================
# RULE BLOCK 3: 3+ charge states excluded globally for specific glycans
# Source: DT 'decision' column "integrating only 2+" + 'check sum spectra' sheet
# =============================================================================

only_2plus_globally <- c(
  "H4N2",    # DT row 6 + check sum spectra: "discarded 3+ in all allotypes"
  "H4N4F2",  # DT row 13 + check sum spectra: "discard all 3+"
  "H5N2",    # DT row 32: "integrating only 2+"
  "H6N2"     # DT row 33: integrate (YF special window), but 3+ excluded globally
)

for (g in only_2plus_globally) {
  df <- nullify(df, glycan = g, charges = 3L)
}

# H4N4F2: additionally absent in YA and YI entirely (both charge states)
# DT row 13: YA = "-", YI = "-"
df <- nullify(df, glycan = "H4N4F2", hc_types = c("YA", "YI"))

# =============================================================================
# RULE BLOCK 4: 3+ charge states excluded for specific allotype × glycan
# Source: 'check sum spectra' sheet individual decisions
# =============================================================================

# H4N3F2: 3+ bad in all; already excluded from integration entirely (Block 1),
# but if any residual values remain, also strip 3+.
df <- nullify(df, glycan = "H4N3F2", charges = 3L)

# H5N4F2: 3+ excluded for Y and YF
# check sum spectra: Y 3+ "low but ok" is superseded by your correction — exclude both
df <- nullify(df, glycan = "H5N4F2", charges = 3L, hc_types = c("Y", "YF"))

# H5N4F2S1: 3+ check for YF — "a bit dirty but we can try" → retain
# No additional nullification needed beyond Block 2 (YA/YI already excluded).

# =============================================================================
# RULE BLOCK 5: H6N2 special integration window for YF
# DT row 33: "IgG_YF - integrate starting 5.6 (dirt)"
# This is a retention-time / integration window issue that must be handled
# at the raw integration stage, not in this post-processing script.
# Flag here for documentation; no code action possible at this stage.
# =============================================================================
# NOTE: IgGIF1H6N2 integration window correction (start 5.6) must be applied
# in the LC-MS software before export. Values currently in the table may
# include dirt signal for YF samples. Review IgGIF1H6N2 manually.
message("MANUAL CHECK REQUIRED: IgGIF1H6N2 — integration window starts at 5.6 for YF samples.")

# =============================================================================
# RULE BLOCK 6: Structures flagged "to check at some point" (DT rows 38–55)
# These were not integrated and have no values; no action needed.
# Listed here for documentation only.
# =============================================================================
# H3N5, H4N5F1S2, H5N6F1S1, H5N5F2S1, H6N4F1, H3N4F2, H7N3F1, H3N6F1,
# H3N4F1S1, H3N5F2, H3N5F1S1, H3N6F2, H4N5F2, H3N6F1S1, H6N5F2, H4N5F3,
# H3N6F2S1, H5N5F4

# =============================================================================
# DONE — save curated data
# =============================================================================

write.csv(df, "./output/EXP02_curation/EXP02_analytes_curated.csv", row.names = FALSE)
message("Analyte curation complete. Output: data/processed/EXP02_analytes_curated.csv")
