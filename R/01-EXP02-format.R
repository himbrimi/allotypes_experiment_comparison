
library(dplyr)
library(tidyr)
library(data.table)

X <- read.csv("./output/EXP02_curation/20260515_EXP02_combined.csv")

### order columns
glycans <- sort(names(X)[grepl("IgGI", names(X))])

X_ <- X[, c( "Sample", "Charge", "genos_id", "allotype", "construct", "type", "num", glycans)]

X <- X_ %>% rename(LC = construct)

write.csv(X, file = "./data/processed/EXP02_curation/20260520_EXP02_combined_formatted.csv")