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
| Replicates | 3 per antibody per experiment (6 when pooled) |

HC types map to analysed glycopeptides as follows: Y → IgGI1, YA → IgGIA1, YF → IgGIF1, YI → IgGIILE1.

## Repository Structure

```
20260513_EXP02_vs_EXP03/
├── data/
│   ├── raw/                          # Original LC-MS output files (read-only)
│   │   ├── 2026-03-18-1359Z_Summary  # EXP03 input
│   │   ├── 2025-04-11-1513Z_Summary  # EXP02 primary run
│   │   ├── 2025-04-11-1541Z_Summary  # EXP02 secondary run
│   │   ├── 2025-04-28-1117Z_Summary  # EXP02 third run
│   │   └── EXP02_samples             # EXP02 sample metadata
│   └── processed/
│       ├── EXP02_curation/           # EXP02 intermediate and curated data
│       └── EXP03_curation/           # EXP03 intermediate and curated data
├── R/
│   ├── 01_data_formatting.R          # Combining and formatting raw LC-MS output
│   ├── 02_analyte_curation.R         # QC-based analyte exclusions and charge state rules
│   └── 03_derived_traits.R           # Derived glycan trait calculations (in progress)
├── output/
│   ├── figures/
│   └── tables/
├── docs/
│   └── EXP02_LCMS_data.xlsx          # Integration decisions and QC notes for EXP02
└── README.md
```

## Analysis Pipeline

### Step 1 — Data formatting (`01_data_formatting.R`)
Raw LC-MS summary files are combined and reformatted into a single wide-format table with one row per sample × charge state.

### Step 2 — Analyte curation (`02_analyte_curation.R`)
Three types of exclusions are applied based on QC review documented in `docs/EXP02_LCMS_data.xlsx`:

- **Contamination flagging:** values in peptide columns that do not correspond to the sample's HC type are removed (carry-over and calibration errors)
- **Global exclusions:** 16 glycan structures not integrated in any allotype are set to NA
- **Allotype-specific exclusions:** structures absent or failed QC in specific HC types
- **Charge state exclusions:** 3+ ions excluded globally or per allotype where signal quality was insufficient

### Step 3 — Trait harmonisation and derived traits (`03_derived_traits.R`)
Individual glycoforms are aggregated into 9 biologically meaningful derived traits. Definitions are harmonised between EXP02 and EXP03 before pooling. See trait definitions below.

### Steps 4–6 — Statistics and visualisation (planned)
- Inter-experiment reproducibility (Spearman ρ, Bland-Altman)
- ART-ANOVA with experiment as random blocking factor
- BH-FDR correction; heatmaps, bar charts, PCA

## Derived Trait Definitions

| Label | Name | Key glycoforms |
|---|---|---|
| G0 | Agalactosylation | H3N3F1, H3N4F1, H3N5F1, H5N3F1 |
| G | Galactosylation | H4N4F1, H4N5F1, H5N4F1, H5N4F2, H5N5F1, H6N3F1 |
| S | Sialylation | H4N4F1S1, H5N4F1S1, H5N4F1S2, H6N3F1S1 |
| M | High mannose | H4N2, H5N2, H6N2 |
| B | Bisection | H3N5F1, H4N5F1, H5N5F1 |
| Fuc | Fucosylation | All fucosylated structures |
| AntennaryF | Antennary fucosylation | H5N4F2 |
| A1 | Monoantennary | H3N3F1 |
| H | Hybrid | H5N3F1, H6N3F1, H6N3F1S1 |

Trait definitions are harmonised between experiments (H5N4F2S1 excluded from S and AntennaryF in both; H4N5F1S1 exclusion from B pending domain decision).

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
- [x] Trait definitions finalised
- [x] EXP02 data formatted and curated
- [ ] EXP03 data formatting and curation
- [ ] Trait harmonisation implemented
- [ ] Inter-experiment reproducibility (Step 2)
- [ ] ART-ANOVA and FDR correction (Steps 3–4)
- [ ] Visualisation (Step 5)
