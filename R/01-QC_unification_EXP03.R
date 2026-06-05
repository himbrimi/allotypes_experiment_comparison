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

save(X, file = "./data/processed/00-X-exp03_v1.RData")

### unification to structure the data in EXP02 and EXP03 similarly



X <- X %>%
  filter(!(grepl("stand", genos_id))) %>%
  filter(!(grepl("blank", genos_id))) %>%
  rename(type = genos_id)

X <- X %>%
  mutate(tmp = HC)


X$tmp <- factor(X$tmp, levels = c("Y",  "YA", "YF" ,"YI"), labels =c("IgGI",  "IgGIA",   "IgGIF",  "IgGIILE"))

X <- X %>%
  mutate(cond = ifelse( tmp != peptide, 8, 20 ))

X <- X %>%
  filter(!(cond == 8))

X <- X %>%
  select(-tmp, -cond)

save(X, file = "./data/processed/00-X-exp03_v1_1_without_stands.RData")


