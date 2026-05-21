library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)




#
# load the data
#


load("./data/processed/EXP03_curation/03-X.RData")

## set missing values to NA
 X[X == 0]<- NA

 
X <- X %>%
  mutate(type = paste(HC, LC,sep ="_"))
# plot according to the sample type

p <- ggplot(data=X, aes(x=type, y=narea))

print(
  p
  + facet_wrap(~glycan, scales="free", ncol=3)
  + geom_jitter(aes(shape = LC, colour = peptide))
  + ylab("% Normalized Area")
  + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5))
  
)

## calculate averages

X_aver <- tmp %>%
  group_by(type, glycan)%>%
  summarise(mean_narea = mean(narea, na.rm = T)) %>%
  ungroup()

Xa_w <- X_aver %>%
  spread(type, mean_narea)

write.csv(Xa_w, file = "./out/data/01-glycans_averages_EXP03.csv")