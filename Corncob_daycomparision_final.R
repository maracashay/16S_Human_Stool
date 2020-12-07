#Code first drafted 10.2.20 by Suzanne Fleishman (suf44@psu.edu)
#Finalized 12/3/20

#### Setup ####

### Clear workspace ###
rm(list=ls())

setwd("CorncobSept2020/Outputs")


### Packages ###
#Corncob Tutorial and more info: https://rdrr.io/github/bryandmartin/corncob/f/vignettes/corncob-intro.Rmd

### Load Packages ###
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("nlme","lme4","corncob","phyloseq", "ggplot2","DESeq2","microbiomeSeq","metagMisc","tidyverse","car","emmeans")
ipak(packages)

### Prep Phyloseq Objects ###
ps<-readRDS("Full-Phyloseq.rds")

#Separate by sample type
ps.stool<-subset_samples(ps, Type =="stool")
ps.trach<-subset_samples(ps, Type =="tracheal")

#upload new metadata indicating if samples have pairs (i.e. two dates of samples)
new.meta<-read.csv("metadata_updatesept.csv")

#before merging in new metadata, double check IDs match
ps.stool@sam_data$Basespace_ID==new.meta$Basespace_ID

row.names(new.meta)<-row.names(ps.stool@otu_table)
ps.stool@sam_data<-sample_data(new.meta)

#subset to only paired data with two dates
ps.stool.pair<-subset_samples(ps.stool, Par == "1")


####For separate Genera test, replicate entire code, but change to genus level####

#ps.stool.gen<-taxa_level(ps.stool.pair,"Genus") #microbiomeSeq package
#ps.stool.pair<-ps.stool.gen

####subset for taxa in 20% of samples####

quartz()
phyloseq_prevalence_plot(ps.stool.pair, prev.trh = NULL, taxcolor = NULL,
                         facet = FALSE, point_alpha = 0.7, showplot = T)



ps.gen.stool.20<- phyloseq_filter_prevalence(ps.stool.pair, prev.trh = 0.2, abund.trh = 0.1,
                                             threshold_condition = "AND", abund.type = "total")

ps.gen.20 = prune_samples(sample_sums(ps.gen.stool.20)>0, ps.gen.stool.20)

#### Recategorize for binary and binned####

#examine output
ps.gen.20@sam_data

ps.gen.stool.20@sam_data$Vent

#adjust values so that zeros indicate healthier

#Vent and ICU, 0-5
ps.gen.20@sam_data$Day<-as.factor(ps.gen.20@sam_data$Day)
ps.gen.20@sam_data$year<-as.factor(ps.gen.20@sam_data$year)


ps.gen.20@sam_data$Vent[ps.gen.20@sam_data$Vent<14.5]=1
ps.gen.20@sam_data$Vent[ps.gen.20@sam_data$Vent>14]=0


ps.gen.20@sam_data$ICU[ps.gen.20@sam_data$ICU<14.1]=1
ps.gen.20@sam_data$ICU[ps.gen.20@sam_data$ICU>14]=0

ps.gen.20@sam_data$Vent<-as.factor(ps.gen.20@sam_data$Vent)
ps.gen.20@sam_data$ICU<-as.factor(ps.gen.20@sam_data$ICU)

#inhosp already binned
ps.gen.20@sam_data$inhosp<-as.factor(ps.gen.20@sam_data$inhosp)


#severity 0,1,2 vs 3
ps.gen.20@sam_data$severity[ps.gen.20@sam_data$severity<1]=0
ps.gen.20@sam_data$severity[ps.gen.20@sam_data$severity==1]=2
ps.gen.20@sam_data$severity[ps.gen.20@sam_data$severity>2]=3
ps.gen.20@sam_data$severity<-as.factor(ps.gen.20@sam_data$severity)

#gender as factor
ps.gen.20@sam_data$Gender<-as.factor(ps.gen.20@sam_data$Gender)


####Full model####
#First, test entire dataset as to control for multiple comparisons across entire set
set.seed(1)
da_DAY_gen<- differentialTest(formula = ~ Gender+Age+Day*year+Day*Vent+Day*ICU+Day*inhosp+Day*severity, #the DA formula - if you have a controlling factor it must be included. 
                                    phi.formula = ~ Gender+Age+Day*year+Day*Vent+Day*ICU+Day*inhosp+Day*severity, #the DV formula - if you have a controlling factor it must be included
                                    formula_null = ~ Gender+Age, # DA controlling factors
                                    phi.formula_null = ~ Gender+Age+Day*year+Day*Vent+Day*ICU+Day*inhosp+Day*severity, #DV controlling factors
                                    test = "Wald", boot = FALSE, #wald test is "standard"
                                    data = ps.gen.20, #PS object
                                    fdr_cutoff = 0.05) #pval false discovery value
#explore significant taxa
sig<-as.data.frame(da_DAY_gen$significant_taxa)
quartz()
plot(da_DAY_gen)

#view significant models and subset only those that include a significant DA interaction with day 
# Done by pasting into excel and filtering, as "summary.bbdml" class does not  convent to dataframe
da_DAY_gen$significant_models


###Model ASVs individually####
#since "differentialTest()" function does not allow for means comparisons of interactions, 
#individually model taxa that had significant interaction terms with "Day"

#Vent: 20,29,31,55, Allistipes
#year: 196,62,65,200, Porphyromonas
#severity: 31, Lachnoclostridium, Prevotella


#duplicate ps object
ps.s<-ps.gen.20

#FOR EACH ASV, compare individually modeled taxa to null model to ensure significance

#null model
null<-bbdml(formula = ASV62 ~ 1,
         phi.formula = ~ 1,
         data = ps.s)

set.seed(1)
#model individual taxa with factors
q<-bbdml(formula = ASV62 ~ Age+Gender+Day*year+Day*Vent+Day*ICU+Day*inhosp+Day*severity,
         phi.formula = ~ Age+Gender+Day*year+Day*Vent+Day*ICU+Day*inhosp+Day*severity,
         data = ps.s)
summary(q)

#compare model with factors to null model
lrtest(mod_null = null, mod = q) 
###all tested ASVs were significant compared to NULL model 


#Obtain predicted mean values from model
p<-plot(q)
p$data$samples==rownames(ps.s@sam_data)
meta<-as.data.frame(ps.s@sam_data)
mod.RA<-p$data

#ensure factors are properly coded
mod.RA$Age<-as.numeric(meta$Age)
mod.RA$severity<-as.factor(meta$severity)
mod.RA$ICU<-as.factor(meta$ICU)
mod.RA$ID<-meta$ID
mod.RA$Day<-as.factor(meta$Day)
mod.RA$ID<-as.factor(meta$ID)
mod.RA$Gender<-as.factor(meta$Gender)
mod.RA$inhosp<-as.factor(meta$inhosp)
mod.RA$Vent<-as.factor(meta$Vent)
mod.RA$year<-as.factor(meta$year)

mod.RA2<-mod.RA
mod.RA2$m<-(mod.RA2$ymin+mod.RA2$ymax)/2


#test  final modeled values with type III sums of squares to ensure significance
aov = aov(m~Gender+Age+Day*year+Day*Vent+Day*ICU+Day*inhosp+Day*severity, data=mod.RA2)
Anova(aov, type = "III") # Type 3 SS for unbalanced designs, package car
#ASV's 196, 65, 25,29,31,55 and Allistipes genus no longer significant

#remaining significant interaction terms:
#ASV62: year, vent
#ASV20: severity
#ASV200: severity
#Lachnoclostridium (genus level): ICU, Severity
#Porphyromonas (genus level): year
#Prevotella (genus level): Severity

#for a particular ASV or Genus, test each variable that has a significant interaction with "Day"

em<-emmeans(aov, pairwise ~ Day*ICU, adjust = "tukey")
#ASV 200 not significant in means comparision between "Day"

figure<-as.data.frame(em$emmeans)

figure$interact<-paste(figure$Day,figure$Vent)

#for each significant interaction with "Day" output a table for each, label a,b,c,d,etc
#a: ASV62=Day*year
#b: ASV20=Day*severity
#c: ASV62=Day*vent

c<-figure

c$ASV<-"ASV62"
c$Genus<-"Bacteroides"
c$variable<-colnames(c)[2]

colnames(c)[2]<-"varval"


##Return to line 132 and repeat with all ASVs

####create a table of all predicted interaction means####

bind1<-rbind(a,b,c)

meanstable<-bind1

write.csv(meanstable,"ASVdayinteractionmeans10.20.20.csv")

####figure####
meanstable$dayASV<-paste(meanstable$varval,meanstable$Day)
quartz()
ggplot(meanstable, aes(x=dayASV, y=emmean, fill=varval))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(x=dayASV, ymin=emmean-SE, ymax=emmean+SE), width=0.1, size=.5)+
  facet_wrap(~variable+Genus+ASV,scales = "free")+#, nrow = 4, ncol = 4,)+
  ylab("Predicted Relative Abundance")+
  xlab("Day")+
  scale_x_discrete(labels= c("1","5","1","5","1","5","1","5",
                             "1","5","1","5","1","5","1","5","1","5","1","5"))
