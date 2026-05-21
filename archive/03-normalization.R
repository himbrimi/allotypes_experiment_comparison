library(dplyr)
library(tidyr)
library(data.table)



X <- read.csv("./data/processed/EXP02_curation/EXP02_analytes_curated_v3.csv", stringsAsFactors = F, sep = ";")


Xl <- X %>%
  gather(glycopeptide, area, 9:(ncol(X)))

Xl <- Xl %>%
  group_by(Sample, glycopeptide, genos_id, allotype, LC, type, num) %>%
  summarise(area_23 = sum(area, na.rm = T)) %>%
  ungroup()

Xl_ <- Xl %>%
  mutate(peptide = gsub("1+.*$", "", Xl$glycopeptide)) %>%
  mutate(glycan = gsub("^IgGI[[:alpha:]]*1", "", Xl$glycopeptide)) %>%
  dplyr::select(-glycopeptide)



X <- Xl_

save(X, file = "./data/processed/EXP02_curation/03-X.RData")