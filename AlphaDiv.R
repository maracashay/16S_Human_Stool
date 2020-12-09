# Alpha-diversity code
set.seed(500)
ps.rare.full <-rarefy_even_depth(ps.stool, sample.size = 10000)

#### Shannon Diversity ####
shannon <- diversity(rareOTU.stool, index="shannon")
shannon.t <- t(shannon)

#### Richness ####
richness <- rowSums(t(rareOTU.stool>0))

#### Evenness (microbiome package) ####
evenness <- evenness(rareOTU.stool)

meta.stool.rare <- data.frame(sample_data(ps.rare.full))

meta.stool.rare$Shannon <- shannon.t
meta.stool.rare$Richness <- richness
meta.stool.rare$Evenness <- evenness$pielou

#t.tests
t.test(shannon_tuk ~ Day, data = met.stool.rare, var.equal = FALSE)

# kruskall-wallis test
m6 <- kruskal.test(meta.stool.rare$Shannon ~ meta.stool.rare$Year)
m6

# linear model
m5 <- lm(shannon_tuk ~ Vent, data = met.trach.df.no.na)
plot(m5)
