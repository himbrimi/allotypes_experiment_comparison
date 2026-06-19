# IgG Glycome Sequence Study — Project Briefing

## ⚠️ Action required

- [ ] **Verify H6N3F1S1 EXP02 charge-state integration** — 3+ ions absent in EXP02 raw CSV but present in EXP03. No documented rationale found in audit table. For large sialylated structures, 3+ is typically the dominant ion. Origin of the absence (non-detection, instrument settings, or data processing artefact) is unknown. **Raw EXP02 data must be inspected before S and H inter-experiment comparisons can be interpreted.**

## Overview

Impact of antibody heavy chain (HC) and light chain (LC) amino acid sequence on IgG glycosylation profiles. Eight antibodies (4 HC × 2 LC, full factorial) expressed in the same cell culture system; analysed by LC-MS of tryptic glycopeptides. Data: % normalised areas per glycoform per sample.

Two independent experiments (EXP02, EXP03) performed ~1 year apart. Instrument servicing and cell metabolic state differ between experiments. Glycoform discrepancies and charge-state integration differences preclude pooling; each experiment is analysed independently.

---

## Experimental Design

- **Layout:** Y × {WT, NO}, YA × {WT, NO}, YF × {WT, NO}, YI × {WT, NO} = 8 antibodies
- **Replicates:** 3 per antibody per experiment
- **Peptides:** IgGI1 (HC=Y), IgGIA1 (HC=YA), IgGIF1 (HC=YF), IgGIILE1 (HC=YI)
- **Questions:** Do HC, LC, or HC×LC interactions drive glycan traits? Which effects replicate?

---

## Raw Data Files

| File | Experiment | Location |
|---|---|---|
| `2026-03-18-1359Z_Summary_corrected_20260603.csv` | EXP03 | `data/raw/` |
| `exp02-all-data-raw_EXCEL_MAC-e3dited_corrected_20260603.csv` | EXP02 | `data/raw/` |
| `EXP02_samples.csv` | EXP02 metadata | `data/raw/` |

Both CSVs corrected 2026-06-03 (charge-state integration errors removed; no new integrations added). Audit tables: `docs/EXP02_charge_state_summary.xlsx`, `docs/EXP03_charge_state_summary_20260603.xlsx`.

---

## Analysis Pipeline

### Step 0 — Normalisation (`00-normalise_EXP02.R`, `00-normalise_EXP03.R`)

Read corrected source CSV → sum charge states per glycopeptide per sample → remove blanks/standards → normalise to sum = 100 per peptide per sample → flag and remove contamination (peptide–HC mismatch) → complete peptide × glycan grid → remove H5N4F2S1 and H5N4F2S2 globally → zero allotype-specific below-threshold structures (see below) → derive HC/LC/type columns → QC plot.

**EXP02** (`00-normalise_EXP02.R`): sample metadata from `EXP02_samples.csv` via `genos_id` key. Zeros H4N4F2, H5N4F2, H5N5F1 for YA.
**EXP03** (`00-normalise_EXP03.R`): metadata parsed from sample names (pattern `HC-LC_rep_...`). Includes NA Charge filter. Zeros H4N4F2 for YA and YI.

Outputs: `data/processed/00-X-exp02_without_stands.RData`, `data/processed/00-X-exp03_without_stands.RData`.
Supersedes: `01-QC_unification_EXP02.R`, `01-QC_unification_EXP03.R` (archived).

### Step 1 — Derived traits (`01-derived_traits_exp02.R`, `01-derived_traits_exp03.R`)

Load Step 0 output → spread to wide format → apply `derived_traits()` function → gather back to long → filter to 8 trait labels → set zeros to NA for plotting → save.

Exclusions applied implicitly via zeroing in Step 0: H4N4F2 excluded from G and AntennaryF for relevant HC types; H4N3F1 not included in any trait formula (co-eluting isomers).

**EXP02 vs EXP03 trait differences:**
- M: EXP02 uses `H5N2 + H6N2` (H4N2 absent — integration error); EXP03 uses `H4N2 + H5N2 + H6N2`
- S and H: same formula both experiments, but charge-state differences affect integrated signal (see Charge-State Bias below)

Outputs: `data/processed/01-X-EXP02.RData`, `data/processed/01-X-EXP03.RData`.

### Step 2 — ART-ANOVA + BH-FDR + emmeans (`02-derived_traits_stats.R`)

Load Step 1 output → filter to 8 traits → `art(narea ~ HC * LC)` per trait (24 tests per experiment) → BH-FDR across all 24 p-values → extract LC marginal emmeans, HC marginal emmeans, HC×LC cell means → QC plot, LC effect plot, slope graph, HC trajectory plot → save results.

**Caveat:** ART interaction tests can have inflated type I error when main effects are very large (Elkin et al. 2021). With 3 replicates per cell, power for the interaction term is limited.

Outputs: `output/tables/02-EXP02-art-anova.csv`, `output/tables/02-EXP02_data_averages.csv`, `output/tables/02-EXP03-art-anova.csv`, `output/tables/02-EXP03_data_averages.csv`.

### Step 3 — Visualisation (`03-combined_effects_plot.R`, `04_trait_heatmap.R`, `05_glycan_heatmap.R`)

**`03-combined_effects_plot.R`:** Re-fits ART models on both experiments; overlays LC and HC marginal emmeans (EXP02 solid/dark, EXP03 dashed/light); produces EXP02 vs EXP03 scatter plots (identity line, colour = trait, shape = HC, fill = LC). HC×LC plots omitted — no interactions significant after BH-FDR in either experiment.

**`04_trait_heatmap.R`:** Row-z-scored heatmap of derived trait means per antibody (HC × LC). Traits affected by inter-experiment discrepancies asterisked. Bidirectional clustering; Y_NO anchored leftmost.

**`05_glycan_heatmap.R`:** Same approach for directly measured glycans; both rows (glycans) and columns (antibodies) clustered. Load paths require adjustment to match Step 1 glycan output object name.

Outputs: `output/figures/03-*`, `output/figures/04_heatmap_exp0*.pdf`, `output/figures/05_glycan_heatmap_exp0*.pdf`.

### Step 4 — Synthesis

Biological conclusions drawn only from effects significant (q < 0.05) in both experiments independently and consistent in direction. Effects significant in one experiment only are reported as exploratory.

---

## Trait Definitions

| Label | Name | Glycoforms (both experiments) | EXP02 difference |
|---|---|---|---|
| A1 | Monoantennary | H3N3F1 | — |
| G0 | Agalactosylation | H3N3F1 + H3N4F1 + H3N5F1 | — |
| G | Galactosylation | H4N4F1 + H4N5F1 + H5N4F1 + H5N4F2 + H5N5F1 + H6N3F1 | — |
| S | Sialylation | H4N4F1S1 + H5N4F1S1 + H5N4F1S2 + H6N3F1S1 | ⚠ See charge-state bias note |
| M | High Mannose | H4N2 + H5N2 + H6N2 | − H4N2 (integration error) |
| B | Bisection | H3N5F1 + H4N5F1 + H5N5F1 | — |
| AntennaryF | Antennary fucosylation | H5N4F2 | — |
| H | Hybrid | H5N3F1 + H6N3F1 + H6N3F1S1 | ⚠ See charge-state bias note |
| S_A2 | Monoantennary sialylation | H4N4F1S1 + H5N4F1S1 + H5N4F1S2 | ⚠ Same charge-state bias as S; H5N4F1S1 absent for YA in both experiments |

S_A2 is an exploratory trait motivated by elevated sialylation in YI driven by monoantennary sialylated structures. Computed in `01-derived_traits_exp02_v2.R` and `01-derived_traits_exp03_v2.R` (kept separate from the main pipeline for now).

Fucosylation dropped: all non-HM structures carry core fucose → perfectly collinear with 100 − M.

**Charge-state bias (verified 2026-06-17):**
- **H5N4F1S1** — EXP02: 2+ and 3+ (Y, YF, YI); EXP03: 2+ only. EXP02 S slightly higher for Y, YF, YI.
- **H6N3F1S1** — EXP02: 2+ only; EXP03: 2+ and 3+ (all HC). EXP03 S and H slightly higher for all HC types.

> ⚠️ **H6N3F1S1 EXP02 — origin of 3+ absence unknown:** 3+ ions are absent from the EXP02 raw CSV with no documented rationale in the audit table. For large sialylated structures, 3+ is typically the dominant ion. Whether this reflects genuine non-detection, instrument settings, or a data processing artefact is unclear. Raw EXP02 data must be inspected before S and H inter-experiment comparisons can be interpreted.
- The two S biases partially offset but do not cancel (different HC types affected). Inter-experiment comparison of S and H requires this caveat. All other traits unaffected.

---

## Glycan Exclusions

**Global (both experiments):**

| Glycan | Excluded from | Reason |
|---|---|---|
| H4N3F1 | All traits | Co-eluting isomer mixture |
| H4N4F2 | G, AntennaryF | Antennary fucose position unconfirmed |
| H5N4F2S1 | All traits | Present in one HC type only per experiment; removed from both CSVs before normalisation |

**Allotype-specific:**

| Glycoform | Experiment | Excluded for | Retained for | Reason |
|---|---|---|---|---|
| H4N4F2 | EXP03 | YA, YI | Y, YF | S/N below threshold |
| H4N4F2 | EXP02 | YA, YI | Y, YF | S/N below threshold; no column for YA in source |
| H5N4F2 | EXP02 | YA, YI | Y, YF | Not detected; 2+ only for Y and YF |
| H5N5F1 | EXP02 | YA | Y, YF, YI | Not detected |

---

## Inter-Experiment Glycoform Discrepancies

| Glycoform | EXP02 | EXP03 | Affected traits | Bias |
|---|---|---|---|---|
| H4N2 | Absent (integration error) | Present | M | M slightly lower in EXP02 |

H5N4F2S1 removed from both experiments — no longer a discrepancy. H5N4F1S1 and H6N3F1S1 present in both but with different charge states integrated (see Charge-State Bias above).

**Traits unaffected (directly comparable):** G0, G, A1, B, AntennaryF.
**Traits requiring caveat:** M (missing glycoform), S and H (charge-state bias).

---

## Key References

- Wobbrock et al. (2011). The Aligned Rank Transform. *CHI 2011.*
- Elkin et al. (2021). ART for Multifactor Contrast Tests. *UIST 2021.*
- Benjamini & Hochberg (1995). Controlling the False Discovery Rate. *JRSS-B.*

---

## Current Status

- [x] Source corrections applied and verified (2026-06-03); audit tables in `docs/` (`EXP02_charge_state_summary.xlsx`, `EXP03_charge_state_summary_20260603.xlsx`)
- [x] Inter-experiment glycoform discrepancies and charge-state differences documented (2026-06-17)
- [x] Step 0: Normalisation scripts written and executed
- [x] Step 1: Derived traits — both experiments; outputs: `01-X-EXP02.RData`, `01-X-EXP03.RData`
- [x] Step 2: ART-ANOVA + BH-FDR + emmeans — both experiments; results in `output/tables/`
- [x] Step 3: Visualisation — combined emmeans plots, scatter comparisons, trait heatmaps produced
- [ ] Step 4: Synthesis
