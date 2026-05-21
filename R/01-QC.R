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


load("./data/raw/00-X-exp03.RData")

# depends on which EXP we are processing I might create HC and LC columns in a different way

# X <- X %>%
#   rename(LC = construct) %>%
#   mutate(HC = gsub("-.*$", "", type))

# X[X$narea == 0,]$narea <- NA

X <- X %>%
  mutate(HC = ifelse(grepl("stand", Sample), "stand", 
                     ifelse(grepl("blank", Sample), "blank", gsub("-.*$","", Sample))))
X <- X %>%
  mutate(genos_id = ifelse(grepl("stand", Sample), "stand", 
                     ifelse(grepl("blank", Sample), "blank", gsub("_.*$","", Sample))))

X <- X %>%
  mutate(LC = ifelse(grepl("stand", genos_id), "stand", 
                     ifelse(grepl("blank", genos_id), "blank", gsub("^.*-","", genos_id))))
  


# plot according to the sample type

p <- ggplot(data=X, aes(x=genos_id, y=narea))

print(
  p
  + facet_wrap(~glycan, scales="free", ncol=3)
  + geom_jitter(aes(shape = LC, colour = peptide))
  + ylab("% Normalized Area")
  + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5))
  
)


