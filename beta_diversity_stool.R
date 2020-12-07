#### set up ####
library(phyloseq)
library(vegan)

#### PERMANOVA ####

adon1s <- adonis(stoolOTU ~ Death, data = stoolmet)
adon2s <- adonis(stoolOTU ~ Day, data=stoolmet)
adon3s <- adonis(stoolOTU ~ month, data=stoolmet)
adon3sb <- adonis(stoolOTU ~ year, data=stoolmet)
adon4s <- adonis(stoolOTU ~ Age, data=stoolmet)
adon5s <- adonis(stoolOTU ~ ARDS, data=stoolmet)
adon5b <- adonis(stoolOTU[which(!is.na(stoolmet$severity)),] ~ severity, data=stoolmet)
stoolmet$severity3 <- stoolmet$severity==3
adon5c <- adonis(stoolOTU[which(!is.na(stoolmet$severity)),] ~ severity3, data=stoolmet)
adon6s <- adonis(stoolOTU ~ Vent, data=stoolmet)
stoolmet$Vent14 <- stoolmet$Vent>14
adon6s.b <- adonis(stoolOTU ~ Vent14, data=stoolmet)
adon7s <- adonis(stoolOTU ~ Gender, data=stoolmet)
adon8s <- adonis(stoolOTU ~ Location, data = stoolmet)

#### PERMDISP ####

otu.bray <- vegdist(stoolOTU, method="bray", na.rm=T)

mod1 <- betadisper(otu.bray, stoolmet$severity)
anova(mod1)

mod3 <- betadisper(otu.bray, stoolmet$Day)
anova(mod3)

mod4 <- betadisper(otu.bray, stoolmet$month)
anova(mod4)

mod5 <- betadisper(otu.bray, stoolmet$Gender)
anova(mod5)

stoolmet$AgeGroup <- stoolmet$Age>50
mod6 <- betadisper(otu.bray, stoolmet$AgeGroup)
anova(mod6.p)

mod7 <- betadisper(otu.bray, stoolmet$ARDS)
anova(mod7)

mod8 <- betadisper(otu.bray, stoolmet$severity==3)
anova(mod8)

mod9 <- betadisper(otu.bray, stoolmet$year)
anova(mod9)

mod10 <- betadisper(otu.bray, stoolmet$Location)
anova(mod10)

mod11 <- betadisper(otu.bray, stoolmet$Vent14)
anova(mod11)

# bootsrap samping for severity

stoolOTU.sev3 <- stoolOTU.20[which(stoolmet.20$severity==3),]
stoolmet.sev3 <- stoolmet.20[which(stoolmet.20$severity==3),]
stoolmet.20$severity3 <- stoolmet.20$severity==3

randPermdispMod <- list()
randPermdispAnova <- 1:1000 
randPermanova <- 1:1000 #1000 is how many times to subsample

for(i in randPermanova) {
  sevRand <- floor(runif(21, min=1, max=108)) #21 is how many samples to include in the random sample
  sev3Rand <- floor(runif(21, min=1, max=21))
  stoolOTU.sevRand <- stoolOTU.20[which(stoolmet.20$severity!=3&!is.na(stoolmet.20$severity))[sevRand],]
  stoolmet.sevRand <- stoolmet.20[which(stoolmet.20$severity!=3&!is.na(stoolmet.20$severity))[sevRand],]
  stoolOTU.sev3Rand <- stoolOTU.sev3[sev3Rand,]
  stoolmet.sev3Rand <- stoolmet.sev3[sev3Rand,]
  
  stoolOTU.Rand <- rbind(stoolOTU.sev3Rand,stoolOTU.sevRand)
  stoolmet.Rand <- rbind(stoolmet.sev3Rand,stoolmet.sevRand)
  
  otu.bray.sev3 <- vegdist(stoolOTU.Rand, method="bray", na.rm=T)
  randPermdispMod[[i]] <- betadisper(otu.bray.sev3, stoolmet.Rand$severity==3)
  randPermdispAnova[i] <- anova(randPermdispMod[[i]])[1,5]
  randPermanova[i] <- adonis(stoolOTU.Rand ~ severity3, data = stoolmet.Rand)$aov.tab[1,6]
}

