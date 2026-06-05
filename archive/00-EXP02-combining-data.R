library(readxl)
library(dplyr)
library(tidyr)
library(data.table)

# set working dir if needed 
# setwd("C:/Users/Olga/iCloudDrive/Documents/Documents - Olga’s MacBook Pro/work/glycomics/GENOS LAB/PROJECTS/2024_IgG_glycans_with_Anika/20260513_EXP02_vs_EXP03")

# load all three data sources

X <- fread(file = "./data/raw/2025-04-11-1513Z_Summary.csv")
X1 <- fread("./data/raw/2025-04-11-1541Z_Summary.txt")

colnames(X1)[1] <- "Sample"
colnames(X1)[2] <- "Charge"

## check if sample list is the same

setdiff(unique(X$Sample), unique(X1$Sample))


setdiff(unique(X1$Sample), unique(X$Sample))

### for now remove the standards because there is an issue with stand_03 integration

X <- X %>%
  filter(!(grepl("stand", Sample)))

X1 <- X1 %>%
  filter(!(grepl("stand", Sample)))


## check which colnames overlap between the 2 documents

overlap <- colnames(X1)[colnames(X1) %in% colnames(X)]

overlap_glycans <- overlap[grepl("IgG", overlap)]

## Columns with the same name in X need to be replaced with the same columns from X1
## X is a datatabe, hence I have to subset it in a stupid way

X <- X[,!(names(X) %in% overlap_glycans),  with = FALSE]

X_ <- merge(X, X1)

## load X2

X2 <- fread("./data/raw/2025-04-28-1117Z_Summary.txt")
## check if sample list is the same

setdiff(unique(X$Sample), unique(X2$Sample))


setdiff(unique(X2$Sample), unique(X$Sample))

### for now remove the standards because there is an issue with stand_03 integration


X2 <- X2 %>%
  filter(!(grepl("stand", Sample)))


## check which colnames overlap between the 2 documents

overlap <- colnames(X2)[colnames(X2) %in% colnames(X_)]

###no overlapping columns with glycans and rename it to X

X <- merge(X_, X2)


# drop columns with IgGII1 
igg2 <- names(X)[grepl("IgGII1", names(X))]
X_ <- X[,!(names(X) %in% igg2),  with = FALSE]

X <- X_

X <- X %>%
  relocate(Sample)

## add sample info
samples <- fread("./data/raw/EXP02_samples.csv")

X <- X %>%
  mutate(genos_id = gsub("_2-.*$", "", Sample))

X <- merge(X, samples,by = "genos_id")

X <- X %>%
  relocate(Sample, Charge,genos_id,allotype,construct)


fwrite(X, file = "./data/processed/EXP02_curation/20260515_EXP02_combined.csv")