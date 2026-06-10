# ==============================================================================
# Derived glycan traits — ART-ANOVA and LC effect visualisation
# ==============================================================================



# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(ARTool)
library(emmeans)



# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

load("./data/processed/01-X-EXP03.Rdata")



# ------------------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------------------

X_dt <- X %>%
  filter(glycan %in% c("A1", "G0", "G", "S","M","B", "AntennaryF", "H"))



X_dt$glycan <- factor(X_dt$glycan, levels =c( "G0", "G", "S","B","M", "AntennaryF", "A1","H"),
                      labels = c("Agalactosylation", "Galactosylation", "Sialylation",  "Bisection", 
                                 "High Mannose", "Antennary fucosylation", "Monoantennary", "Hybrid"))

## wo standards

X_dt <- X_dt %>%
  filter(type != "stand")


X_dt$LC <- as.factor(X_dt$LC)
X_dt$HC <- as.factor(X_dt$HC)
X_dt[is.na(X_dt)] <- 0



# ------------------------------------------------------------------------------
# QC plot
# ------------------------------------------------------------------------------

p <- ggplot(data=X_dt, aes(x=type, y=narea))

print(
  p
  + facet_wrap(~glycan, scales="free", ncol=3)
  + geom_jitter(aes(shape = LC, size = 2, fill = HC))
  + ylab("% Normalized Area")
  + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5))
  + scale_size(guide = "none")
  + theme(strip.text=element_text(size=20))
  + scale_shape_manual(values = c(21,24,8))
  + scale_fill_manual(values= c("black", "grey52", "grey78", "grey90", "white" ))
)




# ------------------------------------------------------------------------------
# ART-ANOVA: function definition
# ------------------------------------------------------------------------------

anova_art_glycans <- function(mydata) {
  art_model <- art(narea~ HC * LC, data = mydata)
  anova_art <- anova(art_model)
  art_data <-  data.frame(anova_art$Term, anova_art$Df, anova_art$Df.res, anova_art$`F value`, anova_art$`Pr(>F)`)
  colnames(art_data) <- c("Term", "Df", "Df.res", "F.value", "Pr")
  art_data
}



# ------------------------------------------------------------------------------
# Emmeans: function definitions
# ------------------------------------------------------------------------------

# LC marginal means (pooled across HC)
get_em_LC <- function(mydata) {
  art_model <- art(narea~ HC * LC, data = mydata)
  as.data.frame(emmeans(artlm(art_model, "LC"), ~ LC))
}

# HC x LC cell means
get_em_HC_LC <- function(mydata) {
  art_model <- art(narea~ HC * LC, data = mydata)
  as.data.frame(emmeans(artlm(art_model, "HC:LC"), ~ HC * LC))
}

# HC marginal means (pooled across LC)
get_em_HC <- function(mydata) {
  art_model <- art(narea~ HC * LC, data = mydata)
  as.data.frame(emmeans(artlm(art_model, "HC"), ~ HC))
}


# ------------------------------------------------------------------------------
# ART-ANOVA: apply across all glycans, FDR correction, save
# ------------------------------------------------------------------------------

anova_results <- X_dt %>%
  group_by(glycan) %>%
  do(anova_art_glycans(.)) %>%
  ungroup()

anova_results$p_adj <- p.adjust(anova_results$Pr, method = 'fdr' )

anova_results <- anova_results %>%
  mutate(sign = ifelse(p_adj < 0.05, "sign", "no"))

write.csv(anova_results, file= "./output/tables/02-EXP03-art-anova.csv")



# ------------------------------------------------------------------------------
# Emmeans: compute for all glycans
# ------------------------------------------------------------------------------

em_LC_dt <- X_dt %>%
  group_by(glycan) %>%
  do(get_em_LC(.)) %>%
  ungroup()

em_LC_HC <- X_dt %>%
  group_by(glycan) %>%
  do(get_em_HC_LC(.)) %>%
  ungroup()


em_HC <- X_dt %>%
  group_by(glycan) %>%
  do(get_em_HC(.)) %>%
  ungroup()


# ------------------------------------------------------------------------------
# Plot: raw data with per-group means (original scale)
# ------------------------------------------------------------------------------

X_summary <- X_dt %>%
  group_by(glycan, HC, LC) %>%
  summarise(mean_narea = mean(narea), .groups="drop")

p <- ggplot(data=X_summary, aes(x=LC, y=mean_narea))

print(
  p
  + facet_wrap(~glycan, scales="free", ncol=3)
  
  # Individual replicates — small, faded
  + geom_point(data=X_dt,
               aes(x=LC, y=narea, shape=HC, fill=LC),
               size=2, alpha=0.4,
               position=position_dodge(width=0.3))
  
  # Means per HC x LC — larger, solid
  + geom_point(aes(shape=HC, fill=LC),
               size=4,
               position=position_dodge(width=0.3))
  
  # Lines connecting means across LC, one line per HC
  + geom_line(aes(group=HC, color=HC),
              linewidth=0.8, alpha=0.8,
              position=position_dodge(width=0.3))
  
  # Grand mean across HC levels — thick dashed black line
  + stat_summary(aes(group=1),
                 fun=mean, geom="line",
                 linewidth=1.5, color="black", linetype="dashed")
  
  + scale_shape_manual(values=c(21,24,23,22))
  + scale_fill_manual(values=c("grey20","grey85"))
  + ylab("% Normalized Area")
  + theme_bw()
  + theme(
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
      strip.text=element_text(size=20))
)



# ------------------------------------------------------------------------------
# Plot: LC marginal means (pooled across HC)
# ------------------------------------------------------------------------------


p <- ggplot(data=em_LC_dt,
            aes(x=LC, y=emmean, group=1))

print(
  p
  + facet_wrap(~glycan, scales="free", ncol=3)
  + geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL),
              alpha=0.15, color=NA, fill="grey40")
  + geom_point(size=4)
  + geom_line(linewidth=1.2)
  + ylab("Aligned rank transformed normalized area")
  + theme_bw()
  + theme(
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
      strip.text=element_text(size=20))
)



# ------------------------------------------------------------------------------
# Plot: HC x LC cell means (slope graph)
# ------------------------------------------------------------------------------

p1 <- ggplot(data = em_LC_HC, aes(x=LC, y=emmean, group=HC, color=HC, shape=HC))

print(
  p1
  + facet_wrap(~glycan, scales="free", ncol=3)
  + geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL, fill=HC),
                alpha=0.15, color=NA,
                position=position_dodge(width=0.3))
  + geom_point(size=4,
               position=position_dodge(width=0.3))
  + geom_line(linewidth=0.8,
              position=position_dodge(width=0.3))
  + scale_shape_manual(values=c(21,24,23,22))
  + ylab("Aligned rank transformed normalized area")
  + theme_bw()
  + theme(
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
      strip.text=element_text(size=20))
)



# ------------------------------------------------------------------------------
# Plot: HC marginal means (pooled across LC)
# ------------------------------------------------------------------------------


p <- ggplot(data=em_HC,
            aes(x=HC, y=emmean, group=1))

print(
  p
  + facet_wrap(~glycan, scales="free", ncol=3)
  + geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL),
                alpha=0.15, color=NA, fill="grey40")
  + geom_point(size=4)
  + geom_line(linewidth=1.2)
  + ylab("Aligned rank transformed normalized area")
  + theme_bw()
  + theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    strip.text=element_text(size=20))
)






# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

save(X_dt, file="./data/processed/02-X_dt_EXP03.RData")

### calculate averages

X_aver <- X_dt %>%
  group_by(type, glycan)%>%
  summarise(mean_narea = mean(narea)) %>%
  ungroup()

Xa_w <- X_aver %>%
  spread(type, mean_narea)

fwrite(Xa_w, file="./output/tables/02-EXP03_data_averages.csv")
