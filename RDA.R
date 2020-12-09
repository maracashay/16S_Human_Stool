decorana(asv.stool.20)
# transform with hellinger
asv.filt.hell <- decostand(asv.stool.20, "hell")

# remove samples with NA
meta.stool.20.no.na <- meta.stool.20 %>% drop_na(c("Death", "Age", "severity", "month", "Location", "APACHE", "Vent", "Day", "Gender"))

#library(dplyr)
meta.stool.filt.rda <- meta.stool.20 
asv.filt.hell.rda <- asv.filt.hell
asv.filt.hell.rda$Samples <- meta.stool.filt.rda$Sample
asv.filt.hell.3 <- asv.filt.hell.rda %>% filter(Samples %in% meta.stool.filt.no.na$Sample)
                          
asv.filt.hell.3$Samples <- NULL

dbrda = rda(asv.filt.hell.3 ~ Death + Age + Vent + APACHE, meta.stool.filt.no.na, distance = 'bray', scale = TRUE, na = na.omit)
# vector fit
vf = envfit(dbrda, meta.stool.filt.no.na, na.rm = TRUE, perm = 999)

vf

summary(dbrda)
