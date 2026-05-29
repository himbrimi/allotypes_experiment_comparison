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
├── data/
│   ├── raw/                                    # Original LC-MS output files (read-only)
│   │   ├── 00-X-exp02.RData                    # EXP02 raw combined data
│   │   └── 00-X-exp03.RData                    # EXP03 raw combined data
│   └── processed/                              # Intermediate and curated RData files
│       ├── 00-X-exp02_v1.RData                 # EXP02 formatted
│       ├── 00-X-exp02_v1_1_without_stands.RData# EXP02 standards removed
│       ├── 00-X-exp03_v1.RData                 # EXP03 formatted
│       ├── 00-X-exp03_v1_1_without_stands.RData# EXP03 standards removed
│       ├── 02-X-EXP02.RData                    # EXP02 curated analytes
│       ├── 02-X-EXP03.RData                    # EXP03 curated analytes
│       ├── 03-X_dt                             # Derived traits EXP03
│       └── 03-X_dt_EXP02                       # Derived traits EXP02
├── docs/
│   └── EXP02_LCMS_data.xlsx                    # Integration decisions and QC notes for EXP02
├── output/
│   ├── figures/
│   │   ├── 01-QC-EXP02.pdf                     # QC plots EXP02
│   │   ├── 01-QC-EXP03.pdf                     # QC plots EXP03
│   │   ├── 02-derived-traits_EXP02.pdf         # Derived trait distributions EXP02
│   │   ├── 02-derived-traits_EXP03.pdf         # Derived trait distributions EXP03
│   │   ├── 03-dt-EXP02.pdf                     # Trait summary EXP02
│   │   ├── 03-dt-EXP03.pdf                     # Trait summary EXP03
│   │   ├── 03-EXP02-HC-effect.pdf              # HC effect plots EXP02
│   │   ├── 03-EXP02-LC-effect.pdf              # LC effect plots EXP02
│   │   ├── 03-EXP02-raw-data.pdf               # Raw data overview EXP02
│   │   ├── 03-EXP03-HC-effect.pdf              # HC effect plots EXP03
│   │   ├── 03-EXP03-LC-effect.pdf              # LC effect plots EXP03
│   │   ├── 03-EXP03-raw-data.pdf               # Raw data overview EXP03
│   │   ├── 04-combined-HC-effect.pdf           # HC marginal emmeans, both experiments
│   │   └── 04-combined-LC-effect.pdf           # LC marginal emmeans, both experiments
│   └── tables/
│       ├── 03-art-anova-EXP02.csv              # ART-ANOVA results EXP02 (with BH-FDR)
│       ├── 03-art-anova-EXP03.csv              # ART-ANOVA results EXP03 (with BH-FDR)
│       ├── 03-EXP02_data_averages.csv          # HC×LC trait means EXP02
│       ├── 03-EXP03_data_averages.csv          # HC×LC trait means EXP03
│       ├── 20260528-comparing_EXP02_EXP03.csv  # Inter-experiment comparison (flat)
│       └── 20260528-comparing_EXP02_EXP03.xlsx # Inter-experiment comparison (with difference formulas)
├── R/
│   ├── 01-QC_unification_EXP02.R               # QC and formatting EXP02
│   ├── 01-QC_unification_EXP03.R               # QC and formatting EXP03
│   ├── 02-derived_traits_exp02.R               # Analyte curation and trait derivation EXP02
│   ├── 02-derived_traits.R                     # Analyte curation and trait derivation EXP03
│   ├── 03-derived_traits_stats.R               # ART-ANOVA, BH-FDR, inter-experiment comparison
│   └── 04-combined_effects_plot.R              # Combined HC and LC emmeans plots (both experiments)
├── project_instructions.md                     # Full analysis plan and decisions log
└── README.md
```

## Analysis Pipeline

### Step 1 — QC and data formatting (`01-QC_unification_EXP02.R`, `01-QC_unification_EXP03.R`)
Each experiment is processed independently. Raw LC-MS summary files are reformatted into a wide-format table per experiment. Contamination flags are applied: values in peptide columns that do not correspond to the sample's HC type are removed (carry-over and calibration errors). Charge state harmonisation is applied based on QC review documented in `docs/EXP02_LCMS_data.xlsx`: 3+ ions excluded globally or per allotype where signal quality was insufficient.

### Step 2 — Analyte curation and derived traits (`02-derived_traits_exp02.R`, `02-derived_traits.R`)
Remaining analyte exclusions are applied per experiment:

- **Global exclusions:** glycan structures not integrated in any allotype are set to NA
- **Allotype-specific exclusions:** structures absent or failed QC in specific HC types

Individual glycoforms are then aggregated into 8 biologically meaningful derived traits. Fucosylation (Fuc) was excluded — all non-high-mannose structures carry core fucose, making it perfectly collinear with 100 − M. Trait definitions are experiment-specific where glycoform availability differs (see Derived Trait Definitions below).

### Step 3 — Statistics and inter-experiment comparison (`03-derived_traits_stats.R`)
ART-ANOVA (`art(trait ~ HC * LC)`) run independently per experiment for all 8 traits (24 tests per experiment). BH-FDR correction applied within each experiment. Results exported to `output/tables/`. Descriptive inter-experiment comparison of HC×LC trait means computed and exported as CSV and annotated XLSX.

### Step 4 — Visualisation (`04-combined_effects_plot.R`)
Emmeans from ART models plotted for both experiments on the same graph. Experiment encoded as line type (EXP02 solid, EXP03 dashed) and colour lightness. HC marginal means and LC marginal means produced as combined figures. HC×LC interaction plots omitted — no interactions were significant after BH-FDR correction in either experiment.

## Derived Trait Definitions

| Label | Name | EXP03 | EXP02 additions |
|---|---|---|---|
| A1 | Monoantennary | H3N3F1 | — |
| G0 | Agalactosylation | H3N3F1 + H3N4F1 + H3N5F1 | — |
| G | Galactosylation | H4N4F1 + H4N5F1 + H5N4F1 + H5N4F2 + H5N5F1 + H6N3F1 | — |
| S | Sialylation | H4N4F1S1 + H5N4F1S1 + H5N4F1S2 + H6N3F1S1 | + H4N5F1S1 + H5N4F2S1 |
| M | High mannose | H4N2 + H5N2 + H6N2 | − H4N2 (integration error) |
| B | Bisection | H3N5F1 + H4N5F1 + H5N5F1 | + H4N5F1S1 |
| AntennaryF | Antennary fucosylation | H5N4F2 | + H5N4F2S1 |
| H | Hybrid | H5N3F1 + H6N3F1 + H6N3F1S1 | — |

H4N5F1S1 could not be reliably quantified in EXP03 and is absent from that dataset. H5N4F2S1 was excluded from EXP03 (S/N < 9).

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
- [x] EXP03 data formatted and curated
- [x] Trait harmonisation implemented
- [x] Inter-experiment comparison — descriptive means per HC×LC (`20260528comparing_EXP02_EXP03.csv`)
- [x] ART-ANOVA per experiment (`03artanovaEXP02.csv`, `03artanovaEXP03.csv`)
- [x] BH-FDR correction (Step 5; embedded in ART-ANOVA output files as `p_adj`)
- [x] Visualisation (Step 6) — combined HC and LC emmeans plots (`04-combined_effects_plot.R`); HC×LC plots omitted — no significant interactions in either experiment
- [ ] Synthesis (Step 7)
