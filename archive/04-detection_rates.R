
library(dplyr)
library(tidyr)
library(stringr)

load("./data/processed/EXP03_curation/03-X.RData")

X <- X%>% 
  mutate(
    HC = str_match(Sample, "_(Y[^-]*)-")[, 2],
    LC = str_match(Sample, "-([A-Z]+)_\\d")[, 2]
  )


detection_rates_03 <- X %>%
  group_by(HC, peptide, glycan) %>%
  summarise(
    n_total    = n(),
    n_detected = sum(!is.na(area_23) & area_23 > 0),
    det_rate   = n_detected / n_total,
    .groups = "drop"
  )

load("./data/processed/EXP02_curation/03-X.RData")

X <- X %>%
  mutate(HC = gsub("-.*$", "", type))

detection_rates_02 <- X %>%
  group_by(HC, peptide, glycan) %>%
  summarise(
    n_total    = n(),
    n_detected = sum(!is.na(area_23) & area_23 > 0),
    det_rate   = n_detected / n_total,
    .groups = "drop"
  )


dr_2 <- detection_rates_02 %>%
  rename(det_rate_2 = det_rate) %>%
  select(-n_total, -n_detected)

dr_3 <- detection_rates_03 %>%
  rename(det_rate_3 = det_rate) %>%
  select(-n_total, -n_detected)


dr <- merge(dr_2, dr_3)

dr <- dr %>%
  mutate(flag = ifelse(det_rate_2 == det_rate_3, 0, 1))


write.csv(dr, file = "./output/tables/04-detection_rates.csv")

