#Code first drafted 10.2.20 by Suzanne Fleishman (suf44@psu.edu)
#Finalized 12/3/20 


#### Setup ####

### Clear workspace ###
rm(list=ls())

setwd("CorncobSept2020/Outputs")


### Packages ###
#Corncob Tutorial and more info: https://rdrr.io/github/bryandmartin/corncob/f/vignettes/corncob-intro.Rmd
#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")
# devtools::install_github("bryandmartin/corncob")
### Load Packages ###
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("corncob","phyloseq", "ggplot2","DESeq2","microbiomeSeq","metagMisc","tidyverse","car","forcats")
ipak(packages)

### Prep Phyloseq Objects ###
ps<-readRDS("Full-Phyloseq.rds")

#Separate by sample type
ps.stool<-subset_samples(ps, Type =="stool")
ps.trach<-subset_samples(ps, Type =="tracheal")

new.meta<-read.csv("metadata_updatesept.csv")

ps.stool@sam_data$Basespace_ID==new.meta$Basespace_ID

row.names(new.meta)<-row.names(ps.stool@otu_table)

ps.stool@sam_data<-sample_data(new.meta)

#change factor types
quartz()
plot(ps.stool@sam_data$inhosp,ps.stool@sam_data$Death)

ps.stool@sam_data$Gender<-as.factor(ps.stool@sam_data$Gender)
ps.stool@sam_data$Day<-as.factor(ps.stool@sam_data$Day)

#change inhosp so 1 is lived and 0 is died in hosp
ps.stool@sam_data$inhosp[ps.stool@sam_data$inhosp==1]="1D"
ps.stool@sam_data$inhosp[ps.stool@sam_data$inhosp==0]="2A"

ps.stool@sam_data$inhosp<-as.factor(ps.stool@sam_data$inhosp)

#change severity to 3 categories 0 = 0, 2 = 1 or 2, 2 = 3
ps.stool@sam_data$severity[ps.stool@sam_data$severity<1]=0
ps.stool@sam_data$severity[ps.stool@sam_data$severity==1]=2
ps.stool@sam_data$severity[ps.stool@sam_data$severity>2]=3
ps.stool@sam_data$severity<-as.factor(ps.stool@sam_data$severity)


#subset for taxa in 20% of samples
ps.stool.20<- phyloseq_filter_prevalence(ps.stool, prev.trh = 0.2, abund.trh = 0.1,
                                         threshold_condition = "AND", abund.type = "total")


ps.gen.stool.20<- phyloseq_filter_prevalence(ps.stool.gen, prev.trh = 0.2, abund.trh = 0.1,
                                             threshold_condition = "AND", abund.type = "total")


####Full model test####
#ASV

set.seed(1)
da_asv<- differentialTest(formula = ~ Death+ICU+Vent+inhosp+severity+Gender+Age+Day+ID, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ Death+ICU+Vent+inhosp+severity+Gender+Age+Day+ID, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ Gender+Age+Day+ID, # DA controlling factors
                          phi.formula_null = ~ Gender+Age+Day+ID, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.stool.20, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#explore significant taxa
quartz()
plot(da_asv)
print(da_asv)


####create table of output####
#since "summary.bbml" class cannot be coerced to a data frame, must examine all models and copy/paste to excel to obtain significance
da_asv$all_models
#in excel:
###create a column for each model output 
###a column for significance (1=significant) 
###a "taxa" column with taxa information, via ps.stool.20@tax_table
#Filter for each main effect of interest and save as .csv separately for plotting (e.g. vent, severity, Death)
#upload tables for each main effect separately and plot
allASV2<-as.data.frame(read.csv("ventASV.csv"))



p = ggplot(data=allASV2,
           aes(x = taxa,y = estimate, ymin = (estimate+SE), ymax = (estimate-SE) ))+
  #scale_y_continuous(breaks=c(-4,-2,0,2,4))+
  #ylim(-4.9,4.9)+
  geom_bar(stat="identity",aes(fill=Phylum))+
  geom_errorbar(aes(x=taxa, ymin=estimate-SE, ymax=estimate+SE),width=0,size=.2)+
  xlab('Genus')+ ylab("Log2fold Change")+
  
  # geom_errorbar(aes(ymin=xmin, ymax=xmax,col=sig),width=0.5,cex=.5)+ 
  facet_wrap(~variable,strip.position="top",scales = "free_x",ncol=16) +
  scale_fill_manual(values = c("#0072B2","#009E73", 
                               "#D55E00" )) + #"#E69F00",
  theme(plot.title=element_text(size=16,face="bold"),
        #axis.text.y=,
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"))+
  #strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold")) +
  coord_flip()
quartz()
p+aes(x = fct_inorder(taxa))+geom_hline(yintercept = 0, color = "black",size=.25) 

