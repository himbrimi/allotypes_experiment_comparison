# ==============================================================================
#
# 00-normalise_EXP02.R
#
# EXP02 — Normalisation, sample annotation, and contamination removal
#
# Input:  data/raw/exp02-all-data-raw_EXCEL_MAC-e3dited_corrected_20260603.csv
#           Corrected LaCyTools summary (charge-state errors fixed 2026-06-03).
#           Single wide-format file; each row is one sample × one charge state.
#         data/raw/EXP02_samples.csv
#           Sample metadata: genos_id, allotype, type, and other descriptors.
#
# Processing:
#   1. Sum 2+ and 3+ ion areas per glycopeptide per sample
#   2. Pivot to long format; parse peptide and glycan from column names
#   3. Replace NAs with 0 (glycoforms not detected in a sample)
#   4. Remove blanks and standards (analysed separately)
#   5. Normalise: each glycoform expressed as % of total area per peptide
#      per sample (sum = 100 within each peptide × sample combination)
#   6. Merge sample metadata
#   7. Flag contamination: glycopeptide signal present in a sample whose
#      HC type does not match that peptide (carry-over / calibration error)
#   8. Remove contaminated measurements
#   9. Complete peptide × glycan grid: insert missing rows as 0 so every
#      sample has an entry for every glycan observed anywhere in the dataset
#  10. Remove H5N4F2S2 and H5N4F2S1 globally (zero/NA across all samples)
#  11. Replace remaining NAs with 0; explicitly zero H4N4F2, H5N4F2, H5N5F1
#      for IgGIA1 (YA) samples (below quantification threshold)
#  12. Unification: rename construct → LC, derive HC, drop bookkeeping
#      columns; save intermediate dataset (_v1)
#  13. QC plot on unified data
#  14. Reconstruct type as HC-LC; save final standards-free dataset
#
# Note: EXP03 normalisation follows a different path because sample identity
#       is encoded in sample names rather than in a separate metadata file.
#
# Outputs:
#   data/processed/00-X-exp02_v1.RData                   — intermediate dataset
#   data/processed/00-X-exp02_v1_1_without_stands.RData  — final clean dataset
#   output/figures/00-QC-EXP02.pdf                       — QC plot
#
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Olga’s MacBook Pro/work/glycomics/GENOS LAB/PROJECTS/2024_IgG_glycans_with_Anika/20260513_EXP02_vs_EXP03")

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

path_csv          <- "./data/raw/exp02-all-data-raw_EXCEL_MAC-e3dited_corrected_20260603.csv"
path_samples      <- "./data/raw/EXP02_samples.csv"
path_rdata        <- "./data/processed/00-X-exp02_without_stands.RData"
path_pdf          <- "./output/figures/00-QC-EXP02.pdf"


# ------------------------------------------------------------------------------
# 1. Load corrected source CSV
# ------------------------------------------------------------------------------

raw <- fread(path_csv)

# Expect columns: Sample, Charge, then one column per glycopeptide (e.g. IgGIA1H3N3F1)


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
# 4. Remove blanks and standards
#    Blanks are discarded entirely.
#    Sample names for blanks/standards are expected to contain "blank" or
#    "stand" (case-insensitive); adjust the pattern below if naming differs.
# ------------------------------------------------------------------------------

X <- X %>% filter(!grepl("stand|blank", Sample, ignore.case = TRUE))


# ------------------------------------------------------------------------------
# 5. Normalise: % area per peptide per sample (sum = 100)
# ------------------------------------------------------------------------------

X <- X %>%
  group_by(Sample, peptide) %>%
  mutate(narea = area / sum(area, na.rm = TRUE) * 100) %>%
  ungroup() %>%
  select(-area)


# ------------------------------------------------------------------------------
# 6. Add sample metadata
#    EXP02 sample names are anonymised run IDs; allotype and group information
#    comes from EXP02_samples.csv via a genos_id key.
#    genos_id is extracted from Sample by removing the replicate/injection
#    suffix (e.g. "id_01_2-E,1_01_55390.raw" -> "id_01").
# ------------------------------------------------------------------------------

samples <- fread(path_samples, header = TRUE, sep = ",", stringsAsFactors = FALSE)

X <- X %>%
  mutate(genos_id = str_extract(Sample, "^id_[0-9]+"))

X <- merge(X, samples, by = "genos_id")


# ------------------------------------------------------------------------------
# 7. Flag contamination
#    A measurement is flagged as contaminated (carry-over or calibration error)
#    when the peptide column does not match the sample's allotype.
#    cond == 20  → clean measurement (retained)
#    cond ==  8  → contamination flag (plotted as ✱, then removed)
# ------------------------------------------------------------------------------

X <- X %>%
  mutate(
    cond = ifelse(str_remove(peptide, "1$") != allotype, 8, 20)
  )

# ------------------------------------------------------------------------------
# 8. Reemove contaminated measurements
# ------------------------------------------------------------------------------

X <- X %>% filter(cond == 20) %>% select(-cond)



# ------------------------------------------------------------------------------
# 9. Complete peptide × glycan grid per sample
#    Some glycoform × allotype combinations have no column in the source CSV
#    (e.g. H4N4F2 for IgGIA1) and are therefore absent from the long-format
#    data entirely after the pivot. complete() inserts the missing rows so that
#    every sample has a row for every peptide × glycan combination observed
#    anywhere in the dataset. Missing values are filled with 0.
# ------------------------------------------------------------------------------

# X <- X %>%
#   group_by(Sample) %>%
#   complete(peptide, glycan, fill = list(narea = 0)) %>%
#   fill(genos_id, type, num, allotype, construct, .direction = "downup") %>%
#   ungroup()

X <- X %>%
  complete(
    nesting(Sample, genos_id, type, num, allotype, construct, peptide),
    glycan,
    fill = list(narea = 0)
  )


# ------------------------------------------------------------------------------
# 10. Remove H5N4F2S2 globally
#    Signal is zero or NA across all sample types; excluded before downstream
#    analysis to avoid contributing spurious zeros to normalised totals.
# ------------------------------------------------------------------------------

X <- X %>% filter(!(glycan %in% c("H5N4F2S2", "H5N4F2S1")))


# ------------------------------------------------------------------------------
# 11. Replace NAs with 0; explicitly zero H4N4F2, H5N4F2, H5N5F1 for YA
#     These structures are absent from IgGIA1 (YA) samples — below
#     quantification threshold or not detected. Set to 0 explicitly rather
#     than left as NA for consistency in downstream analysis and plotting.
#     All remaining NAs are also replaced with 0.
# ------------------------------------------------------------------------------

X <- X %>%
  mutate(narea = replace_na(narea, 0)) %>%
  mutate(narea = ifelse(
    allotype == "IgGIA" & glycan %in% c("H4N4F2", "H5N4F2", "H5N5F1"),
    0, narea
  ))


# ------------------------------------------------------------------------------
# 12. Unification
#     Rename construct to LC; derive HC from the type column (everything before
#     the first hyphen). Drop internal bookkeeping columns not needed downstream.
#     Save intermediate dataset (used for QC plot below).
# ------------------------------------------------------------------------------

X <- X %>%
  rename(LC = construct) %>%
  mutate(HC = gsub("-.*$", "", type)) %>%
  mutate(type = paste(HC, LC, sep = "-")) %>%
  select(-allotype, -num, -genos_id)
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
# 14. Reconstruct type as HC-LC and save final standards-free dataset
# ------------------------------------------------------------------------------


save(X, file = path_rdata)

message("Done.")
message("Intermediate dataset: ", path_rdata_v1)
message("Final dataset:        ", path_rdata)
message("QC plot:              ", path_pdf)
