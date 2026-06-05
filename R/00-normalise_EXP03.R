# ==============================================================================
#
# 00-normalise_EXP03.R
#
# EXP03 — Normalisation, contamination removal, and unification
#
# Input:  data/raw/2026-03-18-1359Z_Summary_corrected_20260603.csv
#           Corrected LaCyTools summary (charge-state errors fixed 2026-06-03).
#           Single wide-format file; each row is one sample × one charge state.
#
# Processing:
#   1. Sum 2+ and 3+ ion areas per glycopeptide per sample
#   2. Pivot to long format; parse peptide and glycan from column names
#   3. Replace NAs with 0 (glycoforms not detected in a sample)
#   4. Parse sample metadata (HC, LC, type) from sample names
#      EXP03 sample names follow the pattern HC-LC_replicate_... (e.g. Y-NO_1_...)
#      so no external metadata file is required (unlike EXP02)
#   5. Remove blanks and standards (analysed separately)
#   6. Normalise: each glycoform expressed as % of total area per peptide
#      per sample (sum = 100 within each peptide × sample combination)
#   7. Flag contamination: glycopeptide signal present in a sample whose
#      HC type does not match that peptide (carry-over / calibration error)
#   8. Remove contaminated measurements
#   9. Complete peptide × glycan grid: insert missing rows as 0 so every
#      sample has an entry for every glycan observed anywhere in the dataset
#  10. Remove H5N4F2S1 and H5N4F2S2 globally (all-zero after source correction)
#  11. Replace remaining NAs with 0; explicitly zero H4N4F2 for IgGIA1 (YA)
#      and IgGIILE1 (YI) samples (below quantification threshold in both)
#  12. Unification: drop bookkeeping columns, finalise HC/LC/type
#  13. QC plot on unified data
#  14. Save final standards-free dataset
#
# Note: EXP02 normalisation follows a different path because sample identity
#       is provided via an external metadata file rather than encoded in
#       sample names.
#
# Outputs:
#   data/processed/00-X-exp03_v1_1_without_stands.RData  — final clean dataset
#   output/figures/00-QC-EXP03.pdf                       — QC plot
#
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Olga's MacBook Pro/work/glycomics/GENOS LAB/PROJECTS/2024_IgG_glycans_with_Anika/20260513_EXP02_vs_EXP03")

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

path_csv   <- "./data/raw/2026-03-18-1359Z_Summary_corrected_20260603.csv"
path_rdata <- "./data/processed/00-X-exp03_v1_1_without_stands.RData"
path_pdf   <- "./output/figures/00-QC-EXP03.pdf"


# ------------------------------------------------------------------------------
# 1. Load corrected source CSV
# ------------------------------------------------------------------------------

raw <- fread(path_csv)

# Remove rows with NA Charge — can result from manual editing artefacts in the
# source file (empty rows, header remnants, etc.)
raw <- raw %>% filter(!is.na(Charge))

# Expect columns: Sample, Charge, then one column per glycopeptide (e.g. IgGI1H3N3F1)


# ------------------------------------------------------------------------------
# 2. Sum charge states (2+ and 3+) per glycopeptide per sample
# ------------------------------------------------------------------------------

X <- raw %>%
  pivot_longer(
    cols      = -c(Sample, Charge),
    names_to  = "glycopeptide",
    values_to = "area"
  ) %>%
  mutate(area = replace_na(area, 0)) %>%
  group_by(Sample, glycopeptide) %>%
  summarise(area = sum(area, na.rm = TRUE), .groups = "drop")


# ------------------------------------------------------------------------------
# 3. Parse peptide type and glycan from glycopeptide column name
#    Column names follow the pattern: IgGI1H3N4F1, IgGIA1H4N4F1, etc.
#    Peptide = prefix up to and including the trailing "1" (e.g. IgGI1, IgGIA1)
#    Glycan  = everything after that prefix        (e.g. H3N4F1)
# ------------------------------------------------------------------------------

X <- X %>%
  mutate(
    peptide = str_extract(glycopeptide, "^IgGI[A-Z]*1"),
    glycan  = str_remove( glycopeptide, "^IgGI[A-Z]*1")
  ) %>%
  select(-glycopeptide)


# ------------------------------------------------------------------------------
# 4. Parse sample metadata from sample names
#    Sample names follow the pattern: HC-LC_replicate_... (e.g. Y-NO_1_1-D,...)
#    HC  = everything before the first hyphen  (Y, YA, YF, YI)
#    LC  = everything between the first hyphen and the first underscore (NO, WT)
#    type = HC-LC (e.g. Y-NO, YA-WT)
# ------------------------------------------------------------------------------

X <- X %>%
  mutate(
    HC   = str_extract(Sample, "^[^-]+"),
    LC   = str_extract(Sample, "(?<=-)[^_]+"),
    type = paste(HC, LC, sep = "-")
  )


# ------------------------------------------------------------------------------
# 5. Remove blanks and standards
#    Standards are retained for independent analysis; blanks discarded entirely.
#    Identified by "stand" or "blank" at the start of the sample name.
# ------------------------------------------------------------------------------

X <- X %>% filter(!grepl("^stand|^blank", Sample, ignore.case = TRUE))


# ------------------------------------------------------------------------------
# 6. Normalise: % area per peptide per sample (sum = 100)
# ------------------------------------------------------------------------------

X <- X %>%
  group_by(Sample, peptide) %>%
  mutate(narea = area / sum(area, na.rm = TRUE) * 100) %>%
  ungroup() %>%
  select(-area)


# ------------------------------------------------------------------------------
# 7. Flag contamination
#    A measurement is flagged as contaminated (carry-over or calibration error)
#    when the glycopeptide does not match the sample's HC type.
#    EXP03 HC values (Y, YA, YF, YI) are mapped to their IgG peptide prefixes
#    for comparison with the parsed peptide column.
#    cond == 20  → clean measurement (retained)
#    cond ==  8  → contamination flag (removed in step 8)
# ------------------------------------------------------------------------------

hc_to_peptide <- c(Y = "IgGI", YA = "IgGIA", YF = "IgGIF", YI = "IgGIILE")

X <- X %>%
  mutate(
    cond = ifelse(str_remove(peptide, "1$") != hc_to_peptide[HC], 8, 20)
  )


# ------------------------------------------------------------------------------
# 8. Remove contaminated measurements
# ------------------------------------------------------------------------------

X <- X %>% filter(cond == 20) %>% select(-cond)


# ------------------------------------------------------------------------------
# 9. Complete peptide × glycan grid per sample
#    Some glycoform × allotype combinations have no column in the source CSV
#    (e.g. H4N4F2 for YA and YI) and are therefore absent from the long-format
#    data entirely after the pivot. complete() inserts the missing rows so that
#    every sample has a row for every peptide × glycan combination observed
#    anywhere in the dataset. Missing values are filled with 0.
# ------------------------------------------------------------------------------

X <- X %>%
  complete(
    nesting(Sample, HC, LC, type, peptide),
    glycan,
    fill = list(narea = 0)
  )


# ------------------------------------------------------------------------------
# 10. Remove H5N4F2S1 and H5N4F2S2 globally
#     H5N4F2S1: all charge states removed during source correction (2026-06-03);
#     column still exists in the CSV with all-zero values.
#     H5N4F2S2: zero or NA across all sample types.
#     Both excluded to avoid spurious zeros in normalised totals.
# ------------------------------------------------------------------------------

X <- X %>% filter(!glycan %in% c("H5N4F2S1", "H5N4F2S2"))


# ------------------------------------------------------------------------------
# 11. Replace remaining NAs with 0; explicitly zero H4N4F2 for YA and YI
#     H4N4F2 is below quantification threshold in IgGIA1 (YA) and IgGIILE1
#     (YI) in EXP03; no column exists for these allotypes in the source CSV.
#     Set to 0 explicitly for consistency in downstream analysis and plotting.
# ------------------------------------------------------------------------------

X <- X %>%
  mutate(narea = replace_na(narea, 0)) %>%
  mutate(narea = ifelse(
    HC %in% c("YA", "YI") & glycan == "H4N4F2",
    0, narea
  ))


# ------------------------------------------------------------------------------
# 12. Unification
#     HC, LC, and type are already derived from sample names above.
#     No further column renaming needed; no bookkeeping columns to drop.
# ------------------------------------------------------------------------------

# (no additional transforms needed for EXP03)


# ------------------------------------------------------------------------------
# 13. QC plot
# ------------------------------------------------------------------------------

pdf(path_pdf, width = 14, height = 20)
p <- ggplot(X, aes(x = type, y = narea))
print(
  p
  + facet_wrap(~ glycan, scales = "free", ncol = 3)
  + geom_jitter(aes(shape = LC, colour = peptide))
  + ylab("% Normalised Area")
  + theme_bw()
  + theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text   = element_text(size = 7),
    legend.position = "bottom"
  )
)
dev.off()


# ------------------------------------------------------------------------------
# 14. Save final standards-free dataset
# ------------------------------------------------------------------------------

save(X, file = path_rdata)

message("Done.")
message("Final dataset: ", path_rdata)
message("QC plot:       ", path_pdf)
