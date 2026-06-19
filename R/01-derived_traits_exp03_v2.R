# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
#
# Derived trait calculation
#
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)



#
# load the data
#

load("./data/processed/00-X-exp03_without_stands.RData")

Xw <- X %>%
  spread(glycan, narea)

Xw[is.na(Xw)] <- 0

# H4N3F1 excluded from G0, G, B: ambiguous mixture of two co-eluting structures;
# core fucosylation is unambiguous regardless of isomer, so included in Fuc only.

# H4N4F2 excluded from AntennaryF and G : structure was not confirmed
# by fragmentation; included in Fuc only as core fucosylation is unambiguous.

#### fucosylation is just everythig that is not HM so i propose removing it

derived_traits <- function(data){
  data <-data %>%
    mutate(A1 = H3N3F1) %>%
    mutate(G0 = H3N3F1 + H3N4F1 + H3N5F1) %>%
    mutate(G = H4N4F1 + H4N5F1+H5N4F1+H5N4F2+H5N5F1+H6N3F1) %>%
    mutate(S = H4N4F1S1 + H5N4F1S1 + H5N4F1S2 + H6N3F1S1) %>%
    mutate(S_A2 = H4N4F1S1 + H5N4F1S1 + H5N4F1S2) %>%
    mutate(M = H4N2+H5N2+H6N2) %>%
    mutate(B = H3N5F1+H4N5F1+H5N5F1) %>%
    mutate(AntennaryF = H5N4F2) %>%
    mutate(H = H5N3F1 + H6N3F1 + H6N3F1S1)
  
}


X_all <- Xw %>%
  group_by(peptide) %>%
  do(derived_traits(.)) %>%
  ungroup()

X_l <- X_all %>%
  gather(glycan, narea, H3N3F1:H)

X_dt <- X_l %>%
  filter(glycan %in% c("A1", "G0", "G", "S", "S_A2", "M","B", "AntennaryF", "H"))




X_dt$glycan <- factor(X_dt$glycan, levels =c( "G0", "G", "S", "S_A2", "B","M", "AntennaryF", "A1","H"),
                      labels = c("Agalactosylation", "Galactosylation", "Sialylation", "Sialylation (A2)",  "Bisection", 
                                 "High Mannose", "Antennary fucosylation", "Monoantennary", "Hybrid"))


X_dt[X_dt$narea == 0,]$narea <- NA

p <- ggplot(data=X_dt, aes(x=type, y=narea))

print(
  p
  + facet_wrap(~glycan, scales="free", ncol=3)
  + geom_jitter(aes(colour = peptide, shape = factor(LC)))
  #+ scale_shape_manual(values = c(8,20))
  + ylab("% Normalized Area")
  + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5))
  
)



X <- X_l

save(X, file="./data/processed/01-X-EXP03_V2.RData")




