# ==============================================================================
# Combined EXP02 + EXP03 — LC and HC effect plots (emmeans from ART models)
# ==============================================================================
# Produces three combined plots analogous to 03-derived_traits_stats.R:
#   1. LC marginal means (pooled across HC) — both experiments overlaid
#   2. HC x LC cell means (slope graph)    — both experiments overlaid
#   3. HC marginal means (pooled across LC) — both experiments overlaid
#
# Experiment is encoded as line type (EXP02 = solid, EXP03 = dashed) and
# colour lightness (EXP02 = dark, EXP03 = light), so HC/LC aesthetics are
# preserved from the original script.
# ==============================================================================



# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ARTool)
library(emmeans)



# ------------------------------------------------------------------------------
# Shared helpers
# ------------------------------------------------------------------------------

glycan_levels  <- c("G0", "G", "S", "B", "M", "AntennaryF", "A1", "H")
glycan_labels  <- c("Agalactosylation", "Galactosylation", "Sialylation",
                    "Bisection", "High Mannose", "Antennary fucosylation",
                    "Monoantennary", "Hybrid")

prepare_data <- function(path) {
  load(path)                              # loads object X into environment
  X %>%
    filter(glycan %in% glycan_levels,
           type   != "stand") %>%
    mutate(
      tmpid = paste(peptide, HC, sep = "_"),
      glycan = factor(glycan, levels = glycan_levels, labels = glycan_labels),
      LC = as.factor(LC),
      HC = as.factor(HC)
    ) %>%
    filter(tmpid %in% c("IgGI_Y", "IgGIA_YA", "IgGIF_YF", "IgGIILE_YI")) %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0)))
}

# ART emmeans extractors (same logic as original script)
get_em_LC <- function(d) {
  m <- art(narea ~ HC * LC, data = d)
  as.data.frame(emmeans(artlm(m, "LC"), ~ LC))
}

get_em_HC_LC <- function(d) {
  m <- art(narea ~ HC * LC, data = d)
  as.data.frame(emmeans(artlm(m, "HC:LC"), ~ HC * LC))
}

get_em_HC <- function(d) {
  m <- art(narea ~ HC * LC, data = d)
  as.data.frame(emmeans(artlm(m, "HC"), ~ HC))
}

run_emmeans <- function(X_dt, getter) {
  X_dt %>%
    group_by(glycan) %>%
    do(getter(.)) %>%
    ungroup()
}



# ------------------------------------------------------------------------------
# Load and prepare both experiments
# ------------------------------------------------------------------------------

X02 <- prepare_data("./data/processed/02-X-EXP02.Rdata")
X03 <- prepare_data("./data/processed/02-X-EXP03.Rdata")



# ------------------------------------------------------------------------------
# Compute emmeans per experiment
# ------------------------------------------------------------------------------

em_LC_02   <- run_emmeans(X02, get_em_LC)    %>% mutate(experiment = "EXP02")
em_LC_03   <- run_emmeans(X03, get_em_LC)    %>% mutate(experiment = "EXP03")

em_HC_LC_02 <- run_emmeans(X02, get_em_HC_LC) %>% mutate(experiment = "EXP02")
em_HC_LC_03 <- run_emmeans(X03, get_em_HC_LC) %>% mutate(experiment = "EXP03")

em_HC_02   <- run_emmeans(X02, get_em_HC)    %>% mutate(experiment = "EXP02")
em_HC_03   <- run_emmeans(X03, get_em_HC)    %>% mutate(experiment = "EXP03")

em_LC    <- bind_rows(em_LC_02,    em_LC_03)
em_HC_LC <- bind_rows(em_HC_LC_02, em_HC_LC_03)
em_HC    <- bind_rows(em_HC_02,    em_HC_03)



# ------------------------------------------------------------------------------
# Shared theme and scales
# ------------------------------------------------------------------------------

exp_linetypes <- c("EXP02" = "solid", "EXP03" = "dashed")
exp_alphas    <- c("EXP02" = 0.20,    "EXP03" = 0.10)   # ribbon transparency

base_theme <- theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text   = element_text(size = 16),
    legend.position = "bottom"
  )



# ------------------------------------------------------------------------------
# Plot 1 — LC marginal means, both experiments
# ------------------------------------------------------------------------------

pdf("./output/figures/04-combined-LC-effect.pdf", width = 12, height = 14)

ggplot(em_LC, aes(x = LC, y = emmean, group = experiment,
                  linetype = experiment, color = experiment)) +
  facet_wrap(~ glycan, scales = "free", ncol = 3) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL,
                  fill = experiment, alpha = experiment),
              color = NA) +
  geom_point(size = 3) +
  geom_line(linewidth = 1.2) +
  scale_linetype_manual(values = exp_linetypes) +
  scale_alpha_manual(values = exp_alphas, guide = "none") +
  scale_color_manual(values = c("EXP02" = "grey20", "EXP03" = "grey60")) +
  scale_fill_manual( values = c("EXP02" = "grey20", "EXP03" = "grey60")) +
  ylab("Aligned rank transformed normalized area") +
  base_theme

dev.off()



# ------------------------------------------------------------------------------
# Plot 2 — HC x LC cell means (slope graph), both experiments
# ------------------------------------------------------------------------------
# 
# pdf("./output/figures/04-combined-LC_HC-effect.pdf", width = 12, height = 14)
# 
# ggplot(em_HC_LC,
#        aes(x = LC, y = emmean, group = interaction(HC, experiment),
#            color = HC, shape = HC, linetype = experiment)) +
#   facet_wrap(~ glycan, scales = "free", ncol = 3) +
#   geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL,
#                   fill = HC, alpha = experiment),
#               color = NA,
#               position = position_dodge(width = 0.3)) +
#   geom_point(size = 3,
#              position = position_dodge(width = 0.3)) +
#   geom_line(linewidth = 0.8,
#             position = position_dodge(width = 0.3)) +
#   scale_shape_manual(values = c(21, 24, 23, 22)) +
#   scale_linetype_manual(values = exp_linetypes) +
#   scale_alpha_manual(values = exp_alphas, guide = "none") +
#   ylab("Aligned rank transformed normalized area") +
#   base_theme
# 
# dev.off()



# ------------------------------------------------------------------------------
# Plot 3 — HC marginal means, both experiments
# ------------------------------------------------------------------------------

pdf("./output/figures/04-combined-HC-effect.pdf", width = 12, height = 14)

ggplot(em_HC, aes(x = HC, y = emmean, group = experiment,
                  linetype = experiment, color = experiment)) +
  facet_wrap(~ glycan, scales = "free", ncol = 3) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL,
                  fill = experiment, alpha = experiment),
              color = NA) +
  geom_point(size = 3) +
  geom_line(linewidth = 1.2) +
  scale_linetype_manual(values = exp_linetypes) +
  scale_alpha_manual(values = exp_alphas, guide = "none") +
  scale_color_manual(values = c("EXP02" = "grey20", "EXP03" = "grey60")) +
  scale_fill_manual( values = c("EXP02" = "grey20", "EXP03" = "grey60")) +
  ylab("Aligned rank transformed normalized area") +
  base_theme

dev.off()
