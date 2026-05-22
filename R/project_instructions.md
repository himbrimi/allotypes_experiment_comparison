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

Both input files are pre-curated long-format datasets (total area normalised, contamination-filtered), produced during individual analysis of each experiment.

| File | Experiment |
|---|---|
| `data/processed/00-X-exp02_v1_1_without_stands.RData` | EXP02 |
| `data/processed/00-X-exp03_v1_1_without_stands.RData` | EXP03 |

---

## Analysis Plan (7 steps)

### Step 1 — Inter-experiment glycoform audit (prerequisite, not a statistical step)

Compare the glycoforms present in the two curated datasets and document all differences in representation between EXP02 and EXP03. This serves three purposes: (1) justifies analysing experiments separately rather than pooling; (2) identifies which derived traits are affected and in which direction; (3) provides a factual report for cell culture teams on inter-experiment variation.

**Output:** a table listing, for each glycoform that differs: its status in each experiment, which traits are affected, and the expected directional bias on those traits. Classify traits as unaffected (directly comparable) or affected (inter-experiment comparison requires explicit bias notation).

**Confirmed discrepancies (from RData inspection):**

| Glycoform | EXP02 | EXP03 | Affected traits | Directional bias |
|---|---|---|---|---|
| H4N5F1S1 | Present | Absent | B, S | B and S slightly higher in EXP02 |
| H5N4F2S1 | Present | Absent (excluded, S/N < 9) | S, AntennaryF | S and AntennaryF slightly higher in EXP02 |

Both datasets renormalized to sum = 100 per sample per peptide. EXP03 contains 19 glycoforms; EXP02 contains 21.

**Traits unaffected by discrepancies:** G0, G, M, A1, H (directly comparable across experiments).
**Traits affected:** S, B, AntennaryF (EXP02 values systematically slightly higher; inter-experiment comparison requires explicit bias notation).

### Step 2 — Derive traits per experiment separately

Calculate all 8 traits independently for EXP02 and EXP03 using experiment-specific glycoform availability. Do not apply a common restricted definition — each experiment uses the best available data for that dataset.

### Step 3 — ART-ANOVA per experiment separately

Aligned Rank Transform ANOVA using the `ARTool` package, run independently for each experiment:

```r
art(trait ~ HC * LC, data = df)
```

Tests HC main effect, LC main effect, and HC×LC interaction for each trait (3 tests × 8 traits = 24 tests per experiment).

**Known caveat:** ART interaction tests can have inflated type I error when main effects are very large (Elkin et al. 2021). If a dominant main effect is found, cross-validate the interaction using standard parametric two-way ANOVA.

**Power note:** With only 3 replicates per HC×LC cell, power is limited — especially for the interaction term. Effects must replicate across both experiments to be considered robust.

### Step 4 — Inter-experiment comparison (descriptive)

For traits unaffected by glycoform discrepancies (G0, G, M, A1, H), compute Spearman ρ between EXP02 and EXP03 antibody means across the 8 antibodies as a reproducibility check. No pass/fail threshold; ρ reported as descriptive evidence of consistency.

For traits affected by missing glycoforms (S, B, AntennaryF), note the directional bias explicitly — do not compute ρ as a comparability measure.

Spearman ρ results are intended for internal use and supplementary material; not a primary paper result.

### Step 5 — FDR correction

Apply Benjamini-Hochberg FDR separately within each experiment across all 24 tests (8 traits × 3 effects: HC, LC, HC×LC). Significance threshold: q < 0.05.

### Step 6 — Visualisation

To be planned after Steps 2–5 are implemented. Anticipated outputs:
- Heatmap of trait means per experiment (rows = traits, columns = 8 antibodies; annotated by HC and LC)
- Per-trait bar charts with individual data points for significant findings, experiments shown separately
- PCA on full trait matrix coloured by HC and LC
- Supplementary reproducibility scatter plots (EXP02 vs EXP03 per trait with ρ) for unaffected traits

### Step 7 — Synthesis

Biological conclusions drawn only from effects that are:
- Statistically significant (q < 0.05) in both experiments independently, **and**
- Consistent in direction in both experiments.

Effects significant in only one experiment are reported as exploratory.

---

## Trait Definitions

Eight derived traits, calculated independently per experiment using available glycoforms.

| Label | Name | Formula | EXP02 additions |
|---|---|---|---|
| G0 | Agalactosylation | H3N3F1 + H3N4F1 + H3N5F1 | — |
| G | Galactosylation | H4N4F1 + H4N5F1 + H5N4F1 + H5N4F2 + H5N5F1 + H6N3F1 | — |
| S | Sialylation | H4N4F1S1 + H5N4F1S1 + H5N4F1S2 + H6N3F1S1 | + H4N5F1S1 + H5N4F2S1 |
| M | High Mannose | H4N2 + H5N2 + H6N2 | — |
| B | Bisection | H3N5F1 + H4N5F1 + H5N5F1 | + H4N5F1S1 |
| AntennaryF | Antennary fucosylation | H5N4F2 | + H5N4F2S1 |
| A1 | Monoantennary | H3N3F1 | — |
| H | Hybrid | H5N3F1 + H6N3F1 + H6N3F1S1 | — |

Fucosylation trait dropped: all non-high-mannose structures in this dataset carry core fucose, making Fuc perfectly collinear with 100 − M.

---

## Glycan Structure Exclusions

### Common to both experiments

| Glycan | Excluded from | Reason |
|---|---|---|
| H4N3F1 | G0, G, B | Co-eluting isomer mixture |
| H4N4F2 | G, AntennaryF | Antennary fucose position unconfirmed |

### Experiment-specific exclusions

| Glycoform | EXP02 | EXP03 | Reason |
|---|---|---|---|
| H4N5F1S1 | Present | Absent | Not detected |
| H5N4F2S1 | Present | Excluded | S/N = 5.5, below threshold of 9 |

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
- [x] Glycoform discrepancies between experiments documented (Step 1 complete)
- [x] Pre-curated RData inputs identified, cleaned, and verified
- [x] Step 2: Derive traits per experiment (script ready)
- [ ] Step 3: ART-ANOVA per experiment
- [ ] Step 4: Inter-experiment comparison
- [ ] Step 5: BH-FDR correction
- [ ] Step 6: Visualisation
- [ ] Step 7: Synthesis
