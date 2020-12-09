# mvabund

spp.mva <- mvabund(ps.20@otu_table)
# check mean variance plot
meanvar.plot(spp.mva)

mod1 <- manyglm(spp.mva, meta$Vent)
plot(mod1)

anova(mod1, p.uni = "adjusted)
