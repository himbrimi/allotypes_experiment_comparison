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

glycan_levels  <- c("Agalactosylation", "Galactosylation", "Sialylation", "Sialylation (A2)",
                    "Bisection", "High Mannose", "Antennary fucosylation",
                    "Monoantennary", "Hybrid")
glycan_labels  <- c("Agalactosylation", "Galactosylation", "Sialylation","Sialylation (A2)",
                    "Bisection", "High Mannose", "Antennary fucosylation",
                    "Monoantennary", "Hybrid")

prepare_data <- function(path) {
  load(path)                              # loads object X into environment
  X_dt %>%
    filter(glycan %in% glycan_levels,
           type   != "stand") %>%
    mutate(
      tmpid = paste(peptide, HC, sep = "_"),
      glycan = factor(glycan, levels = glycan_levels, labels = glycan_labels),
      LC = as.factor(LC),
      HC = as.factor(HC)
    ) %>%
    filter(tmpid %in% c("IgGI1_Y", "IgGIA1_YA", "IgGIF1_YF", "IgGIILE1_YI")) %>%
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

X02 <- prepare_data("./data/processed/02-X_dt_EXP02_S_A2.RData")
X03 <- prepare_data("./data/processed/02-X_dt_EXP03_S_A2.RData")



# ------------------------------------------------------------------------------
# Compute emmeans per experiment
# ------------------------------------------------------------------------------

em_LC_02   <- run_emmeans(X02, get_em_LC)    %>% mutate(experiment = "EXP 1")
em_LC_03   <- run_emmeans(X03, get_em_LC)    %>% mutate(experiment = "EXP 2")

em_HC_LC_02 <- run_emmeans(X02, get_em_HC_LC) %>% mutate(experiment = "EXP 1")
em_HC_LC_03 <- run_emmeans(X03, get_em_HC_LC) %>% mutate(experiment = "EXP 2")

em_HC_02   <- run_emmeans(X02, get_em_HC)    %>% mutate(experiment = "EXP 1")
em_HC_03   <- run_emmeans(X03, get_em_HC)    %>% mutate(experiment = "EXP 2")

em_LC    <- bind_rows(em_LC_02,    em_LC_03)
em_HC_LC <- bind_rows(em_HC_LC_02, em_HC_LC_03)
em_HC    <- bind_rows(em_HC_02,    em_HC_03)



# ------------------------------------------------------------------------------
# Shared theme and scales
# ------------------------------------------------------------------------------

exp_linetypes <- c("EXP 1" = "solid", "EXP 2" = "dashed")
exp_alphas    <- c("EXP 1" = 0.20,    "EXP 2" = 0.10)   # ribbon transparency

base_theme <- theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text   = element_text(size = 16),
    legend.position = "bottom"
  )



# ------------------------------------------------------------------------------
# Plot 1 — LC marginal means, both experiments
# ------------------------------------------------------------------------------

pdf("./output/figures/03-combined-LC-effect.pdf", width = 12, height = 14)

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
  scale_color_manual(values = c("EXP 1" = "grey20", "EXP 2" = "grey60")) +
  scale_fill_manual( values = c("EXP 1" = "grey20", "EXP 2" = "grey60")) +
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

pdf("./output/figures/03-combined-HC-effect.pdf", width = 12, height = 14)

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
  scale_color_manual(values = c("EXP 1" = "grey20", "EXP 2" = "grey60")) +
  scale_fill_manual( values = c("EXP 1" = "grey20", "EXP 2" = "grey60")) +
  ylab("Aligned rank transformed normalized area") +
  base_theme

dev.off()


# ------------------------------------------------------------------------------
# Plot 4 — correlation between EXP02 and EXP02
# ------------------------------------------------------------------------------

X02_m <- X02 %>%
  summarise(across(where(is.numeric), mean), .by = c(glycan, LC, HC)) %>%
  mutate(experiment = "EXP1")

X03_m <- X03 %>%
  summarise(across(where(is.numeric), mean), .by = c(glycan, LC, HC)) %>%
  mutate(experiment = "EXP2")


X_m <- bind_rows(X02_m, X03_m)

X_mw <- X_m %>%
  pivot_wider(names_from = experiment, values_from = narea)

X_mw <- X_mw %>%
  mutate(genotype = paste(HC, LC, sep = "-"))

#pdf("./output/figures/04-combined-HC-effect.pdf", width = 12, height = 14)

glycan_colours <- c(
  "Agalactosylation"       = "#E69F00",
  "Galactosylation"        = "#56B4E9",
  "Sialylation"            = "#009E73",
  "Sialylation (A2)" = "#6D5A9C",
  "Bisection"              = "#F0E442",
  "High Mannose"           = "#0072B2",
  "Antennary fucosylation" = "#D55E00",
  "Monoantennary"          = "#CC79A7",
  "Hybrid"                 = "#000000"
)


hc_shapes <- c(Y = 21, YA = 22, YF = 23, YI = 24)

pdf("./output/figures/03-EXP02_vs_EXP03_traits.pdf", width = 12, height = 14)
p <- ggplot(X_mw, aes(x = EXP1, y = EXP2))

print(
  p 
    # WT: filled symbols
    + geom_point(
      data = subset(X_mw, LC == "WT"),
      aes(color = glycan, fill = glycan, shape = HC),
      size = 3, stroke = 0.8
    )
  # NO: hollow symbols (same border colour, fill = NA)
    + geom_point(
      data = subset(X_mw, LC == "NO"),
      aes(color = glycan, shape = HC),
      fill = NA,
      size = 3, stroke = 0.8
    )
  # Invisible points to inject an LC legend
   +  geom_point(aes(alpha = LC), shape = NA, size = 0)
  
  # Scales
    + scale_shape_manual(name = "HC", values = hc_shapes)
    + scale_color_manual(
      name   = "Trait",
      values = glycan_colours,
      guide  = guide_legend(
        override.aes = list(shape = 21, size = 3,
                            fill  = unname(glycan_colours))))
    + scale_fill_manual(values = glycan_colours, guide = "none")
    + scale_alpha_manual(
      name   = "LC",
      values = c(WT = 1, NO = 1),       # alpha has no visual effect
      guide  = guide_legend(
        override.aes = list(
          shape    = c(21, 21),
          color    = c("grey40", "grey40"),
          fill     = c("grey40", NA),   # WT filled, NO hollow
          size     = 3,
          stroke   = 0.8,
          alpha    = 1
        )
      )
    )
    + labs(
      x = "Mean normalised area EXP 1 (%)",
      y = "Mean normalised area EXP 2 (%)"
    )
      + theme_bw(base_size = 11) 
      + theme(legend.key = element_rect(fill = NA)) 
  ## identity line
   + geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey60", linewidth = 0.5)
)

dev.off()

#### separately per IgG heavy chain

pdf("./output/figures/03-EXP02_vs_EXP03_traits_per_HC.pdf", width = 12, height = 12)


p <- ggplot(X_mw, aes(x = EXP1, y = EXP2))

print(
  p 
  # WT: filled symbols
  + geom_point(
    data = subset(X_mw, LC == "WT"),
    aes(color = glycan, fill = glycan, shape = LC),
    size = 3, stroke = 0.8
  )
  # NO: hollow symbols (same border colour, fill = NA)
  + geom_point(
    data = subset(X_mw, LC == "NO"),
    aes(color = glycan, shape = LC),
    fill = NA,
    size = 3, stroke = 0.8
  )
  + facet_wrap(~HC, scales="free", ncol=2)
  + scale_color_manual(values = glycan_colours)
  # Invisible points to inject an LC legend
  #+  geom_point(aes(alpha = LC), shape = NA, size = 0)
  
  # Scales
  # + scale_shape_manual(name = "HC", values = hc_shapes)
  # + scale_color_manual(
  #   name   = "Trait",
  #   values = glycan_colours,
  #   guide  = guide_legend(
  #     override.aes = list(shape = 21, size = 3,
  #                         fill  = unname(glycan_colours))))
  # + scale_fill_manual(values = glycan_colours, guide = "none")
  # + scale_alpha_manual(
  #   name   = "LC",
  #   values = c(WT = 1, NO = 1),       # alpha has no visual effect
  #   guide  = guide_legend(
  #     override.aes = list(
  #       shape    = c(21, 21),
  #       color    = c("grey40", "grey40"),
  #       fill     = c("grey40", NA),   # WT filled, NO hollow
  #       size     = 3,
  #       stroke   = 0.8,
  #       alpha    = 1
  #     )
  #   )
  # )
  + labs(
    x = "Mean normalised area EXP 1 (%)",
    y = "Mean normalised area EXP 2 (%)"
  )
  + theme_bw(base_size = 11) 
  #+ theme(legend.key = element_rect(fill = NA)) 
  ## identity line
  + geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey60", linewidth = 0.5)
)
dev.off()




