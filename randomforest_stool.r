if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("phyloseq", "DESeq2"))

library(phyloseq)
library(tidyverse)
library(randomForest)

rf <- function(variable, data){
  a = as_tibble(sample_data(data))
  a_variable <- a %>% select(variable)
  model_rf <- as.data.frame(cbind(otu_table(data), a_variable))
  x <- model_rf[,1:ncol(model_rf)-1]
  y <- model_rf[,ncol(model_rf)] #Here, you can identify which colume will be used as classifier
  dataset_rf <- as.data.frame(cbind(x,y))
  dim(dataset_rf)
  set.seed(1)
  bestMtry <- tuneRF(x,y, mtryStart = 7,stepFactor = 1.5, improve = 1e-5, ntree = 1001)
  if(is.factor(dataset_rf$y) == TRUE){ 
  dataset_rf$y <- as.factor(dataset_rf$y)}
  else {dataset_rf$y <- as.numeric(dataset_rf$y)}
  rf <- randomForest(y ~ ., data= dataset_rf, ntree=10001, mtry=4, importance=TRUE)
  varimp <- varImpPlot(rf)
  varimp <- as.data.frame(varimp)
  varimp_plot <- as.data.frame(cbind(varimp, tax_table(data)))
  if(is.factor(dataset_rf$y) ==TRUE){
    varimp_plot <-subset(varimp_plot, varimp_plot$MeanDecreaseAccuracy >0)
    varimp_plot <- data.frame(predictors = rownames(varimp_plot), varimp_plot$Phylum,
                              varimp_plot$Family, 
                              varimp_plot$Genus, 
                              varimp_plot$MeanDecreaseAccuracy)
    varimp_plot <- cbind(varimp_plot, 1:nrow(varimp_plot))
    varimp_plot$predictors <- factor(varimp_plot$predictors, 
                            levels = varimp_plot$predictors[order(varimp_plot$varimp_plot.MeanDecreaseAccuracy, decreasing = FALSE)])
    plot<-ggplot(varimp_plot, aes(x=varimp_plot$predictors, y = varimp_plot$varimp_plot.MeanDecreaseAccuracy, fill= varimp_plot$varimp_plot.Phylum)) + 
      geom_bar(stat = "identity",width=0.7) + 
      coord_flip()+ 
      xlab("ASV") +
      ylab("Mean Decrease Accuracy") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
            text=element_text(size=14, color="black")) 
    p1 <- plot+theme(legend.position = "none")+ scale_fill_manual(values = c("#0072B2", "#009E73", "#D55E00"))
    p1
    }
  else {varimp_plot <-subset(varimp_plot, varimp_plot$`%IncMSE` >0)
  varimp_plot <- data.frame(predictors = rownames(varimp_plot), varimp_plot$Phylum,
                            varimp_plot$Family, 
                            varimp_plot$Genus, 
                            varimp_plot$`%IncMSE`)
  varimp_plot <- cbind(varimp_plot, 1:nrow(varimp_plot))
  colnames(varimp_plot)[6] <- "y"
  varimp_plot$predictors <- factor(varimp_plot$predictors, 
                          levels = varimp_plot$predictors[order(varimp_plot$varimp_plot...IncMSE., decreasing = FALSE)])
  plot<-ggplot(varimp_plot, aes(x=varimp_plot$predictors, y = varimp_plot$varimp_plot...IncMSE., fill= varimp_plot$varimp_plot.Phylum)) + 
    geom_bar(stat = "identity",width=0.7) + 
    coord_flip()+ 
    xlab("ASV") +
    ylab("%Increased Purity") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
          text=element_text(size=14, color="black")) 
  p1 <- plot+theme(legend.position = "none")+ scale_fill_manual(values = c("#0072B2", "#009E73", "#D55E00"))
  p1}
  print(p1)
  }

a <- rf("Vent", ps.perc.20)
b <- rf("ARDS", ps.perc.20)
c <- rf("inhosp", ps.perc.20)
d <- rf("severity", ps.perc.20)
e <- rf("Death", ps.perc.20)
