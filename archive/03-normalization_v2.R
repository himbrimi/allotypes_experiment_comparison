library(dplyr)
library(tidyr)
library(data.table)



X <- fread("./data/raw/2026-03-18-1359Z_EXP03.csv", dec=",")


Xl <- X %>%
  gather(glycopeptide, area, 3:(ncol(X)))

Xl <- Xl %>%
  group_by(Sample, glycopeptide ) %>%
  summarise(area_23 = sum(area, na.rm = T)) %>%
  ungroup()

Xl_ <- Xl %>%
  mutate(peptide = gsub("1+.*$", "", Xl$glycopeptide)) %>%
  mutate(glycan = gsub("^IgGI[[:alpha:]]*1", "", Xl$glycopeptide)) %>%
  dplyr::select(-glycopeptide)



X <- Xl_

X <- X %>%
  group_by(Sample, peptide) %>% 
  mutate(narea=area_23/sum(area_23, na.rm = T)*100) %>% 
  ungroup()
  

## add sample information, drop total area
X <- X %>%
  select(-area_23) %>%
  mutate(
    HC = str_match(Sample, "_(Y[^-]*)-")[, 2],
    LC = str_match(Sample, "-([A-Z]+)_\\d")[, 2])


save(X, file = "./data/processed/EXP03_curation/03-X.RData")