# IgG Glycome Sequence Study — Project Briefing

## Overview

This project analyses the impact of antibody heavy chain (HC) and light chain (LC) amino acid sequence on IgG glycosylation profiles. Eight antibodies were expressed in the same cell culture system and analysed by LC-MS of tryptic glycopeptides. Data are expressed as % normalized areas per glycoform per sample.

Two independent experiments were performed approximately one year apart (EXP02 and EXP03). The LC-MS instrument underwent service procedures between experiments and cells may have been in a different metabolic state. Differences in detected glycoforms between experiments are documented below and preclude pooling; each experiment is analysed independently.

---

## Experimental Design

- **Layout:** 4 Heavy Chains × 2 Light Chains = 8 antibodies (full factorial design)
- **Replicates:** 3 measurements per antibody per experiment
- **Peptides analysed:** IgGI1 (HC=Y), IgGIA1 (HC=YA), IgGIF1 (HC=YF), IgGIILE1 (HC=YI)
- **Primary scientific questions:**
  - Do HC sequence, LC sequence, or their interaction drive specific glycan traits?
  - Which glycan trait effects replicate consistently across both experiments?

---

## Raw Data Files

The corrected LaCyTools summary CSVs are the authoritative source inputs for this analysis. Both were manually corrected on 2026-06-03 and verified by systematic charge-state audit. Charge-state integration summaries are documented in `docs/`.

| File | Experiment | Location |
|---|---|---|
| `2026-03-18-1359Z_Summary_corrected_20260603.csv` | EXP03 | `data/raw/` |
| `exp02-all-data-raw_EXCEL_MAC-e3dited_corrected_20260603.csv` | EXP02 | `data/raw/` |

These replace the previously used `00-X-exp03_without_stands.RData` and `00-X-exp02_without_stands.RData` as starting points. Renormalised RData files will be regenerated from these CSVs by `R/00-normalise_EXP03.R` and `R/00-normalise_EXP02.R` (Step 0).

### EXP03 source corrections (2026-06-03)

The original LaCyTools output contained charge-state integration errors. All corrections involve removal of erroneously included ions; no new integrations were added.

| Glycoform | Correction | Scope | Rationale |
|---|---|---|---|
| H3N5F1 | Removed 3+ ions | IgGIA1 (YA) only | Inconsistent with all other HC types where 3+ was not integrated |
| H5N2 | Removed 3+ ions | All HC types | 3+ ions not used for this structure in any sample type |
| H6N3F1 | Removed 3+ ions | IgGIILE1 (YI) only | Inconsistent with all other HC types where 3+ was excluded |
| H4N4F2 | Removed 2+ ions | IgGIA1 (YA) and IgGIILE1 (YI) | S/N below threshold of 9; retained for IgGI1 (Y) and IgGIF1 (YF) |
| H5N4F1S1 | Removed all ions | All HC types | Signal determined to be noise upon re-inspection; structure absent from EXP03 entirely |
| H5N4F2S1 | Removed all ions | IgGI1 (Y) only | Present in Y only; removed for cross-experiment comparability (absent from EXP02 in all HC types) |

**Verification:** all corrections confirmed by systematic charge-state audit (Python; 2026-06-03). Charge-state integration summary: `docs/EXP03_charge_state_summary_v2.xlsx`.

### EXP02 source corrections (2026-06-03)

The original LaCyTools output also contained charge-state integration errors. All corrections involve removal of erroneously included ions.

| Glycoform | Correction | Scope | Rationale |
|---|---|---|---|
| H5N4F1S1 | Removed all ions | IgGIA1 (YA) only | S/N below threshold; retained for Y, YF, YI |
| H4N4F2 | Removed all ions | IgGIILE1 (YI) only | S/N below threshold of 9; retained for Y and YF (YA has no column) |
| H5N4F2 | 2+ only retained | IgGI1 (Y) and IgGIF1 (YF) | 3+ excluded; absent for YA and YI |
| H6N3F1S1 | Removed 3+ ions | All HC types | 3+ excluded globally |
| H5N4F2S1 | Absent | IgGIILE1 (YI) only, zeroed | Present in YI only; removed for cross-experiment comparability |

**Verification:** all corrections confirmed by systematic charge-state audit (Python; 2026-06-03). Charge-state integration summary: `docs/EXP02_charge_state_summary.xlsx`.

### H5N4F2S1 — removed from both experiments

H5N4F2S1 was present in only one HC type per experiment (Y in EXP03; YI in EXP02, zeroed after correction). It contributes to the S and AntennaryF traits but cannot be compared across experiments. It has been removed from both corrected CSVs and will be excluded during normalisation (Step 0). This supersedes the earlier entry in the inter-experiment discrepancy table.

---

## Analysis Plan (8 steps)

### Step 0 — Renormalisation from corrected source CSVs (new)

**Inputs:**
- `data/raw/2026-03-18-1359Z_Summary_corrected_20260603.csv` (EXP03)
- `data/raw/exp02-all-data-raw_EXCEL_MAC-e3dited_corrected_20260603.csv` (EXP02)

For each experiment: read corrected LaCyTools summary CSV, remove standards, exclude H5N4F2S1 columns, and renormalise absolute areas to sum = 100 per sample per peptide. Export renormalised standards-free datasets as:
- `data/processed/00-X-exp03_without_stands.RData`
- `data/processed/00-X-exp02_without_stands.RData`

Scripts: `R/00-normalise_EXP03.R` and `R/00-normalise_EXP02.R`.

### Step 1 — Inter-experiment glycoform audit (prerequisite, not a statistical step)

Compare the glycoforms present in the two curated datasets and document all differences in representation between EXP02 and EXP03. This serves three purposes: (1) justifies analysing experiments separately rather than pooling; (2) identifies which derived traits are affected and in which direction; (3) provides a factual report for cell culture teams on inter-experiment variation.

**Output:** a table listing, for each glycoform that differs: its status in each experiment, which traits are affected, and the expected directional bias on those traits. Classify traits as unaffected (directly comparable) or affected (inter-experiment comparison requires explicit bias notation).

**Confirmed discrepancies (post-correction):**

| Glycoform | EXP02 | EXP03 | Affected traits | Directional bias |
|---|---|---|---|---|
| H4N5F1S1 | Present | Absent (could not be reliably quantified) | B, S | B and S slightly higher in EXP02 |
| H5N4F1S1 | Present (Y, YF, YI) | Absent (noise; source correction 2026-06-03) | S | S slightly higher in EXP02 |
| H4N2 | Absent (integration error; removed globally) | Present | M | M slightly lower in EXP02 |

H5N4F2S1 removed from both experiments — no longer a discrepancy.

**Traits unaffected by discrepancies:** G0, G, A1, H (directly comparable across experiments).
**Traits affected:** S, B, M (inter-experiment comparison requires explicit bias notation). AntennaryF now unaffected (H5N4F2S1 removed from both).

### Step 2 — Derive traits per experiment separately

Calculate all 8 traits independently for EXP02 and EXP03 using experiment-specific glycoform availability. Do not apply a common restricted definition — each experiment uses the best available data for that dataset.

**Note for EXP03:** H5N4F1S1 absent from EXP03 (noise; source correction) — excluded from S trait. H4N4F2 excluded for YA and YI (allotype-specific; already excluded from G and AntennaryF globally).

**Note for EXP02:** H5N4F1S1 absent for YA (allotype-specific exclusion). H4N4F2 absent for YI. H5N4F2S1 excluded from both experiments.

### Step 3 — ART-ANOVA per experiment separately

Aligned Rank Transform ANOVA using the `ARTool` package, run independently for each experiment:

```r
art(trait ~ HC * LC, data = df)
```

Tests HC main effect, LC main effect, and HC×LC interaction for each trait (3 tests × 8 traits = 24 tests per experiment).

**Known caveat:** ART interaction tests can have inflated type I error when main effects are very large (Elkin et al. 2021). If a dominant main effect is found, cross-validate the interaction using standard parametric two-way ANOVA.

**Power note:** With only 3 replicates per HC×LC cell, power is limited — especially for the interaction term. Effects must replicate across both experiments to be considered robust.

### Step 4 — Inter-experiment comparison (descriptive)

For traits unaffected by glycoform discrepancies (G0, G, A1, H, AntennaryF), compute Spearman ρ between EXP02 and EXP03 antibody means across the 8 antibodies as a reproducibility check. No pass/fail threshold; ρ reported as descriptive evidence of consistency.

For traits affected by missing glycoforms (S, B, M), note the directional bias explicitly — do not compute ρ as a comparability measure.

Spearman ρ results are intended for internal use and supplementary material; not a primary paper result.

### Step 5 — FDR correction

Apply Benjamini-Hochberg FDR separately within each experiment across all 24 tests (8 traits × 3 effects: HC, LC, HC×LC). Significance threshold: q < 0.05.

### Step 6 — Visualisation

Implemented in `04-combined_effects_plot.R`. Emmeans from ART models plotted for both experiments on the same graph, with experiment encoded as line type (EXP02 solid, EXP03 dashed) and colour lightness. Two plot types produced:

- **HC marginal means** (pooled across LC) — `04-combined-HC-effect.pdf`
- **LC marginal means** (pooled across HC) — `04-combined-LC-effect.pdf`

HC×LC cell means code is present in the script but commented out — no HC×LC interactions were significant after BH-FDR correction in either experiment.

### Step 7 — Synthesis

Biological conclusions drawn only from effects that are:
- Statistically significant (q < 0.05) in both experiments independently, **and**
- Consistent in direction in both experiments.

Effects significant in only one experiment are reported as exploratory.

---

## Trait Definitions

Eight derived traits, calculated independently per experiment using available glycoforms.

| Label | Name | EXP03 | EXP02 additions |
|---|---|---|---|
| A1 | Monoantennary | H3N3F1 | — |
| G0 | Agalactosylation | H3N3F1 + H3N4F1 + H3N5F1 | — |
| G | Galactosylation | H4N4F1 + H4N5F1 + H5N4F1 + H5N4F2 + H5N5F1 + H6N3F1 | — |
| S | Sialylation | H4N4F1S1 + H5N4F1S2 + H6N3F1S1 | + H4N5F1S1 + H5N4F1S1 (Y, YF, YI only) |
| M | High Mannose | H4N2 + H5N2 + H6N2 | − H4N2 (integration error) |
| B | Bisection | H3N5F1 + H4N5F1 + H5N5F1 | + H4N5F1S1 |
| AntennaryF | Antennary fucosylation | H5N4F2 | — (H5N4F2S1 removed from both experiments) |
| H | Hybrid | H5N3F1 + H6N3F1 + H6N3F1S1 | — |

Fucosylation trait dropped: all non-high-mannose structures in this dataset carry core fucose, making Fuc perfectly collinear with 100 − M.

H5N4F1S1: absent from EXP03 (noise; source correction 2026-06-03). Present in EXP02 for Y, YF, YI; absent for YA (below QC threshold).
H5N4F2S1: removed from both experiments — present in only one HC type per experiment, precluding meaningful comparison.

---

## Glycan Structure Exclusions

### Common to both experiments

| Glycan | Excluded from | Reason |
|---|---|---|
| H4N3F1 | G0, G, B | Co-eluting isomer mixture |
| H4N4F2 | G, AntennaryF | Antennary fucose position unconfirmed |
| H5N4F2S1 | All traits | Present in only one HC type per experiment; removed from both CSVs before normalisation (2026-06-03) |

### Experiment-specific exclusions

| Glycoform | EXP02 | EXP03 | Reason |
|---|---|---|---|
| H4N5F1S1 | Present | Absent | Could not be reliably quantified in EXP03 |
| H5N4F1S1 | Present (Y, YF, YI); absent (YA) | Absent (all HC types) | EXP03: noise upon re-inspection (source correction 2026-06-03). EXP02: below QC threshold for YA |
| H4N2 | Absent | Present | Integration error in EXP02 (detected in IgGIF1 only); removed globally |

### Allotype-specific exclusions

| Glycoform | Experiment | Excluded for | Retained for | Reason |
|---|---|---|---|---|
| H4N4F2 | EXP03 | YA, YI | Y, YF | S/N below threshold in YA and YI |
| H4N4F2 | EXP02 | YI | Y, YF | S/N below threshold of 9; no column for YA |
| H5N4F1S1 | EXP02 | YA | Y, YF, YI | S/N below threshold |
| H5N4F2 | EXP02 | YA, YI | Y, YF | Not detected; 2+ only retained for Y and YF |

---

## R Environment

- Language: R
- Key packages: `ARTool`, `ggplot2`, `pheatmap`, `tidyverse`, `rstatix`, `data.table`
- Reproducibility: `renv`
- Version control: GitHub

---

## Key References

- Bland & Altman (1986). Statistical methods for assessing agreement. *The Lancet.*
- Wobbrock et al. (2011). The Aligned Rank Transform. *CHI 2011.*
- Elkin et al. (2021). ART for Multifactor Contrast Tests. *UIST 2021.*
- Benjamini & Hochberg (1995). Controlling the False Discovery Rate. *JRSS-B.*

---

## Current Status

- [x] Analysis plan defined
- [x] Trait definitions finalised per experiment (8 traits; Fuc dropped)
- [x] EXP03 source corrections applied and verified (2026-06-03)
- [x] EXP02 source corrections applied and verified (2026-06-03)
- [x] H5N4F2S1 removed from both experiments (present in one HC type only per experiment)
- [x] Charge-state integration summaries produced (`docs/EXP03_charge_state_summary_v2.xlsx`, `docs/EXP02_charge_state_summary.xlsx`)
- [x] Glycoform discrepancies between experiments documented and updated (Step 1)
- [ ] **Step 0: Renormalisation scripts** — `R/00-normalise_EXP03.R` and `R/00-normalise_EXP02.R` to be written; new RData files to be generated
- [ ] Step 2: Derive traits — to be rerun for both experiments after Step 0
- [ ] Step 3: ART-ANOVA — to be rerun for both experiments after Step 0
- [ ] Step 4: Inter-experiment comparison — to be rerun after Step 0
- [ ] Step 5: BH-FDR correction — to be rerun after Step 0
- [ ] Step 6: Visualisation — to be rerun after Step 0
- [ ] Step 7: Synthesis
