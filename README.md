# IgG Glycome Sequence Study
**Project:** 2024_IgG_glycans_with_Anika  
**Analysis folder:** 20260513_EXP02_vs_EXP03

## Overview

This project investigates how antibody heavy chain (HC) and light chain (LC) amino acid sequence influences IgG glycosylation profiles. Eight antibodies (4 HC × 2 LC, full factorial design) were expressed in the same cell culture system and analysed by LC-MS of tryptic glycopeptides. Data are expressed as % normalised areas per glycoform per sample.

Two independent experiments were included: **EXP02** and **EXP03**, performed approximately one year apart. Some variation in detected glycoforms between experiments is expected due to instrument servicing and differences in cell metabolic state.

## Experimental Design

| Factor | Levels |
|---|---|
| Heavy chain (HC) | Y, YA, YF, YI |
| Light chain (LC) | WT, NO |
| Antibodies | 8 (full factorial) |
| Replicates | 3 per antibody per experiment |

HC types map to analysed glycopeptides as follows: Y → IgGI1, YA → IgGIA1, YF → IgGIF1, YI → IgGIILE1.

## Repository Structure

```
20260513_EXP02_vs_EXP03/
├── 20260611_EXP02_vs_EXP03_OlgaZaitseva.pptx  # Summary presentation
├── Note_IMPORTANT.docx                          # Project notes
├── project_instructions.md                      # Full analysis plan and decisions log
├── README.md
├── archive/                                     # Superseded scripts (retained for reference)
│   ├── 00-EXP02-combining-data.R
│   ├── 01-EXP02-format.R
│   ├── 02_analyte_curation.R
│   ├── 03-normalization.R
│   ├── 03-normalization_v2.R
│   ├── 04-detection_rates.R
│   └── 05-QC.R
├── data/
│   ├── raw/                                     # Corrected LaCyTools source CSVs (read-only)
│   │   ├── 00-X-exp03.RData                     # EXP03 raw combined data (original)
│   │   ├── 01-X-exp02.RData                     # EXP02 raw combined data (original)
│   │   ├── 2026-03-18-1359Z_Summary_corrected_20260603.csv
│   │   ├── exp02-all-data-raw_EXCEL_MAC-e3dited_corrected_20260603.csv
│   │   └── EXP02_samples.csv                    # EXP02 sample metadata
│   └── processed/
│       ├── 00-X-exp02_without_stands.RData      # EXP02 renormalised, clean
│       ├── 00-X-exp03_without_stands.RData      # EXP03 renormalised, clean
│       ├── 01-X-EXP02.RData                     # EXP02 with derived traits
│       ├── 01-X-EXP03.RData                     # EXP03 with derived traits
│       ├── 02-X_dt_EXP02.RData                  # EXP02 traits, stats-ready
│       └── 02-X_dt_EXP03.RData                  # EXP03 traits, stats-ready
├── docs/
│   ├── 20260528comparing_EXP02_EXP03.csv        # Inter-experiment comparison (flat)
│   ├── EXP02_charge_state_summary.xlsx          # Charge-state audit EXP02
│   ├── EXP02_LCMS_data.xlsx                     # Integration decisions and QC notes EXP02
│   ├── EXP02_vs_EXP03_integration_differences.xlsx
│   └── EXP03_charge_state_summary_20260603.xlsx # Charge-state audit EXP03
├── output/
│   ├── figures/
│   │   ├── 00-QC-EXP02.pdf                      # QC plot EXP02
│   │   ├── 00-QC-EXP03.pdf                      # QC plot EXP03
│   │   ├── 01-EXP02-dt.pdf                      # Derived trait distributions EXP02
│   │   ├── 01-EXP03-dt.pdf                      # Derived trait distributions EXP03
│   │   ├── 02-EXP02_dt_bw.pdf                   # Trait QC plot EXP02, b/w
│   │   ├── 02-EXP02-HC_trajectories.pdf         # HC marginal emmeans EXP02
│   │   ├── 02-EXP02-LC_effect.pdf               # LC marginal emmeans EXP02
│   │   ├── 02-EXP02-LC_effect_per_HC.pdf        # HC×LC cell means EXP02
│   │   ├── 02-EXP02-LC_slope_graph.pdf          # HC×LC slope graph EXP02
│   │   ├── 02-EXP03_dt_bw.pdf                   # Trait QC plot EXP03, b/w
│   │   ├── 02-EXP03-HC_trajectories.pdf         # HC marginal emmeans EXP03
│   │   ├── 02-EXP03-LC_effect.pdf               # LC marginal emmeans EXP03
│   │   ├── 02-EXP03-LC_effect_per_HC.pdf        # HC×LC cell means EXP03
│   │   ├── 02-EXP03-LC_slope_graph.pdf          # HC×LC slope graph EXP03
│   │   ├── 03-combined-HC-effect.pdf            # HC marginal emmeans, both experiments
│   │   ├── 03-combined-LC-effect.pdf            # LC marginal emmeans, both experiments
│   │   ├── 03-EXP02_vs_EXP03_traits.pdf         # EXP02 vs EXP03 scatter, all traits
│   │   ├── 03-EXP02_vs_EXP03_traits_per_HC.pdf  # EXP02 vs EXP03 scatter, faceted by HC
│   │   ├── 04_heatmap_exp02.pdf                 # Derived trait z-score heatmap EXP02
│   │   ├── 04_heatmap_exp03.pdf                 # Derived trait z-score heatmap EXP03
│   │   ├── 05_glycan_heatmap_exp02.pdf          # Raw glycan heatmap EXP02
│   │   └── 05_glycan_heatmap_exp03.pdf          # Raw glycan heatmap EXP03
│   └── tables/
│       ├── 02-EXP02-art-anova.csv               # ART-ANOVA results EXP02 (with BH-FDR)
│       ├── 02-EXP02_data_averages.csv           # HC×LC trait means EXP02
│       ├── 02-EXP03-art-anova.csv               # ART-ANOVA results EXP03 (with BH-FDR)
│       ├── 02-EXP03_data_averages.csv           # HC×LC trait means EXP03
│       ├── combined-art-anova-replication.xlsx  # Consolidated replication table
│       └── combined-art-anova-replication_PRESENT.xlsx
└── R/
    ├── 00-normalise_EXP02.R                     # Normalisation, QC, unification EXP02
    ├── 00-normalise_EXP03.R                     # Normalisation, QC, unification EXP03
    ├── 01-derived_traits_exp02.R                # Trait derivation EXP02
    ├── 01-derived_traits_exp03.R                # Trait derivation EXP03
    ├── 02-derived_traits_stats.R                # ART-ANOVA, BH-FDR, emmeans
    ├── 03-combined_effects_plot.R               # Combined emmeans and scatter plots
    ├── 04_trait_heatmap.R                       # Derived trait heatmap
    └── 05_glycan_heatmap.R                      # Raw glycan heatmap
```

## Analysis Pipeline

### Step 0 — Normalisation, QC, contamination removal, and unification (`00-normalise_EXP02.R`, `00-normalise_EXP03.R`)

Both LaCyTools source files were manually corrected for charge-state integration errors (see Source Corrections below). Each corrected CSV is read; charge states summed per glycopeptide per sample; blanks and standards removed; absolute areas renormalised to sum = 100 per sample per peptide; contamination flagged and removed; peptide × glycan grid completed; globally absent structures removed; HC and LC columns derived. QC plot saved as PDF. Output: `00-X-exp02_without_stands.RData` and `00-X-exp03_without_stands.RData`.

These scripts incorporate and supersede the previously separate `01-QC_unification_EXP02.R` and `01-QC_unification_EXP03.R`.

### Step 1 — Analyte curation and derived traits (`01-derived_traits_exp02.R`, `01-derived_traits.R`)

Remaining analyte exclusions applied per experiment:

- **Global exclusions:** H5N4F2S1 (present in one HC type only per experiment)
- **Allotype-specific exclusions:** structures absent or below QC threshold in specific HC types (see Allotype-Specific Exclusions below)

Individual glycoforms aggregated into 8 derived traits. Fucosylation excluded — all non-high-mannose structures carry core fucose, making it perfectly collinear with 100 − M.

### Step 2 — Statistics and inter-experiment comparison (`02-derived_traits_stats.R`)

ART-ANOVA (`art(trait ~ HC * LC)`) run independently per experiment for all 8 traits (24 tests per experiment). BH-FDR correction applied within each experiment. Results exported to `output/tables/`. Descriptive inter-experiment comparison of HC×LC trait means computed and exported as CSV and annotated XLSX.

### Step 3 — Visualisation (`03-combined_effects_plot.R`, `04_trait_heatmap.R`)

Emmeans from ART models plotted for both experiments on the same graph. Experiment encoded as line type (EXP02 solid, EXP03 dashed) and colour lightness. HC marginal means and LC marginal means produced as combined figures. HC×LC interaction plots omitted — no interactions were significant after BH-FDR correction in either experiment. Trait heatmaps produced separately.

## Source Corrections (2026-06-03)

Both source CSVs were audited by systematic charge-state consistency checking across all glycoform × HC type combinations. Corrections involve removal of erroneously included ions only; no new integrations were added. Full audit tables in `docs/`.

### EXP03 corrections

| Glycoform | Correction | Scope | Rationale |
|---|---|---|---|
| H3N3F1, H3N5F1, H4N2, H4N3F1, H4N5F1, H5N2, H6N2,  H6N3F1 | Removed 3+ ions |All HC types | In some allotypes 3+ ions did not pass the QC,so 3+ ions were removed for all HC types for consistency |
| H4N4F2 | Removed all ions | IgGIA1 (YA) and IgGIILE1 (YI) | S/N below threshold of 9 |
| H5N4F1S1 | Removed all ions | IgGIA1 (YA)  | Did not pass the QC |
| H5N4F2S1 | Removed all ions | was integrated for IgGI1 (Y) only | removed  because S/N below threshold of 9 |
| H5N4F1S2 | Removed all ions | IgGIA1 (YA) and IgGIILE1 (YI) | Did not pass the QC |
| H5N5F1| Removed all ions | IgGIA1 (YA) and IgGIILE1 (YI) | Did not pass the QC |


### EXP02 corrections

| Glycoform | Correction | Scope | Rationale |
|---|---|---|---|
| H4N4F2 | Removed all ions | IgGIA1 (YA) and IgGIILE1 (YI)| S/N below threshold of 9 |
| H5N4F2 | Retained 2+ only | IgGI1 (Y) and IgGIF1 (YF) | 3+ excluded; absent for YA and YI |
| H5N4F2S1 | Removed all ions | All HC types |Did not pass the QC |
| H3N3F1, H3N5F1, H4N3F1, H4N5F1, H5N2, H6N2, H6N3F1, H6N3F1S1 | Removed 3+ ions | All HC types | 3+ excluded globally |


## Derived Trait Definitions

| Label | Name | Glycoforms (both experiments) | EXP02 difference |
|---|---|---|---|
| A1 | Monoantennary | H3N3F1 | — |
| G0 | Agalactosylation | H3N3F1 + H3N4F1 + H3N5F1 | — |
| G | Galactosylation | H4N4F1 + H4N5F1 + H5N4F1 + H5N4F2 + H5N5F1 + H6N3F1 | — |
| S | Sialylation | H4N4F1S1 + H5N4F1S2 + H6N3F1S1 + H5N4F1S1 | ⚠ See note below |
| M | High mannose | H4N2 + H5N2 + H6N2 | − H4N2 (integration error) |
| B | Bisection | H3N5F1 + H4N5F1 + H5N5F1 | — |
| AntennaryF | Antennary fucosylation | H5N4F2 | — |
| H | Hybrid | H5N3F1 + H6N3F1 + H6N3F1S1 | ⚠ See note below |

**Inter-experiment charge-state differences (verified 2026-06-17 against audit tables):**

Two glycoforms contributing to S and H are integrated with different charge states across experiments. Glycoform formulas are identical but the integrated signal differs:

- **H5N4F1S1** — EXP02: 2+ and 3+ (Y, YF, YI); EXP03: 2+ only (Y, YF, YI). EXP02 S values are slightly higher for these three HC types due to the additional 3+ contribution.
- **H6N3F1S1** — EXP02: 2+ only (all HC types); EXP03: 2+ and 3+ (all HC types). EXP03 S and H values are slightly higher for all HC types due to the additional 3+ contribution.

These two biases partially offset each other in the S trait but do not cancel cleanly (they affect different HC types). Inter-experiment comparisons of S and H should be interpreted with this caveat. All other traits are unaffected by charge-state differences.



## Dependencies

```r
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(ARTool)
library(lme4)
library(BlandAltmanLeh)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(rstatix)
```

Reproducibility managed with `renv`. Run `renv::restore()` to install the correct package versions.

## Key References

- Bland & Altman (1986). Statistical methods for assessing agreement. *The Lancet.*
- Wobbrock et al. (2011). The Aligned Rank Transform. *CHI 2011.*
- Elkin et al. (2021). ART for Multifactor Contrast Tests. *UIST 2021.*
- Benjamini & Hochberg (1995). Controlling the False Discovery Rate. *JRSS-B.*

## Status

- [x] Analysis plan defined
- [x] Trait definitions finalised (8 traits; Fuc dropped; AntennaryF updated — H5N4F2S1 removed)
- [x] EXP03 source corrections applied and verified (2026-06-03)
- [x] EXP02 source corrections applied and verified (2026-06-03)
- [x] H5N4F2S1 removed from both experiments (present in one HC type only per experiment)
- [x] Charge-state integration audit tables produced (`docs/EXP03_charge_state_summary_v2.xlsx`, `docs/EXP02_charge_state_summary.xlsx`)
- [x] Inter-experiment charge-state differences identified (2026-06-17): S and H traits affected; all other traits unaffected
- [x] Inter-experiment glycoform discrepancies documented (Step 1 prerequisite)
- [x] Step 0 scripts written and executed — `R/00-normalise_EXP02.R`, `R/00-normalise_EXP03.R`; new RData files generated; supersede `01-QC_unification_*.R`
- [x] Step 1: Trait derivation — rerun for both experiments (`R/01-derived_traits.R`, `R/01-derived_traits_exp02.R`); outputs: `data/processed/01-X-EXP02.RData`, `data/processed/01-X-EXP03.RData`
- [x] Step 2: ART-ANOVA, BH-FDR, emmeans — rerun (`R/02-derived_traits_stats.R`); results: `output/tables/02-EXP02-art-anova.csv`, `output/tables/02-EXP02_data_averages.csv`, `output/tables/02-EXP03-art-anova.csv`, `output/tables/02-EXP03_data_averages.csv`
- [x] Step 3: Visualisation — `R/03-combined_effects_plot.R` and `R/04_trait_heatmap.R` run; combined emmeans plots, inter-experiment scatter, and trait heatmaps produced
- [ ] Synthesis
