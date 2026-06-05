# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
#
# quality control
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


load("./data/raw/01-X-exp02.RData")

# depends on which EXP we are processing I might create HC and LC columns in a different way

X <- X %>%
  rename(LC = construct) %>%
  mutate(HC = gsub("-.*$", "", type))

X[X$narea == 0,]$narea <- NA


# plot according to the sample type

p <- ggplot(data=X, aes(x=type, y=narea))

print(
  p
  + facet_wrap(~glycan, scales="free", ncol=3)
  + geom_jitter(aes(shape = LC, colour = peptide))
  + ylab("% Normalized Area")
  + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5))
  
)

X <- X %>%
  select(-allotype, -num, -genos_id, -cond)
save(X, file = "./data/processed/00-X-exp02_v1.RData")


X <- X %>%
  filter(!(grepl("stand", type))) 

X <- X %>%
  filter(!(grepl("blank", type))) 

X <- X %>%
  mutate(type = paste(HC, LC, sep = "-"))



save(X, file = "./data/processed/00-X-exp02_v1_1_without_stands.RData")
