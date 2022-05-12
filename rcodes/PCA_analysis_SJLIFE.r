######################
## PCA analysis     ##
## Date: 05/10/2022 ##
## Achal Neupane    ##
######################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/PCA")
PCA <- read.table("final-PCAS.eigenvec", header =T, stringsAsFactors=FALSE)
dim(PCA)
thousandG.ethnicty <- read.table("integrated_call_samples_v3.20130502.ALL.panel", header = T , stringsAsFactors = FALSE)
head(thousandG.ethnicty)


PCA$Population <- "SJLIFE"
PCA$Population <- thousandG.ethnicty$super_pop[match(PCA$IID, thousandG.ethnicty$sample)]
PCA <- PCA[c("FID", "IID", c(paste0("PC", 1:10), "Population"))]
PCA$KEY <- paste(PCA$FID, PCA$IID, sep =":")
PCA$Population <- as.character(PCA$Population)
PCA$Population[is.na(PCA$Population)] <- "SJLIFE"
# write.table(PCA, "SJLIFE_PCA_with_sample_population.txt", sep ="\t", col.names = T, quote = F)

## Change the target to reorder the dots on top of the PCA plot
target <- c("AFR", "AMR", "EAS", "EUR", "SAS", "SJLIFE")
PCA$Population <- factor(PCA$Population, levels = target)
PCA <- PCA[order(-as.numeric(factor(PCA$Population))),]

library(ggplot2)
## Generate a new file that has IID, PC1,PC2, and a new column Population 
p <- ggplot(PCA, aes(x=PC1, y=PC2, color=Population)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("SJLIFE_4481") +
  scale_color_manual(values = c(AFR='violet', AMR='blue', EAS='red', EUR="black", SAS="pink", SJLIFE="brown")) +
  annotate("text", x=0.02, y=-0.01, label="AFR", size=4, color = "violet") +
  annotate("text", x=0.005, y=0.024, label="AMR", size=4, color = "blue") +
  annotate("text", x=0.011, y=0.033, label="EAS", size=4, color = "red") +
  annotate("text", x=-0.008, y=0.004, label="EUR", size=4, color = "black") +
  annotate("text", x=0.0072, y=0.013, label="SAS", size=4, color = "pink") +
  annotate("text", x=-0.01, y=0, label="SJLIFE", size=4, color = "brown") +
  theme_classic() +
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face = "plain", size=12),
        axis.text=element_text(face = "plain"),axis.title = element_text(face = "plain"),plot.title = element_text(face = "plain", hjust = 0.5,size=13))

p


# ggsave("SJLIFE-ALL-Population.jpg", plot = p, device = NULL, scale = 1, width = 12, height = 8, dpi = 600, limitsize = TRUE)



## SD FILTER
# Samples within SD cutoff in reference to CEU HAPMAP samples
NHW_SAMPLES_EUR <- PCA[PCA$Population == "EUR",]
AA_SAMPLES_AFR <- PCA[PCA$Population == "AFR",]
ASIAN_SAMPLES_EAS <- PCA[PCA$Population == "EAS",]

p.sd <- p

## (HAPMAP)
## Select NHW samples only
library(tidyverse)

## NHW
SDSelection.Table.NHW <- list()
SD.cutoff.all <- 3:5
for (i in 1:length(SD.cutoff.all)){
  SD.cutoff <- SD.cutoff.all[i]  
  PC1min <- (mean(NHW_SAMPLES_EUR$PC1) - (SD.cutoff*sd(NHW_SAMPLES_EUR$PC1)))
  PC1max <- (mean(NHW_SAMPLES_EUR$PC1) + (SD.cutoff*sd(NHW_SAMPLES_EUR$PC1)))
  PC2min <- (mean(NHW_SAMPLES_EUR$PC2) - (SD.cutoff*sd(NHW_SAMPLES_EUR$PC2)))
  PC2max <- (mean(NHW_SAMPLES_EUR$PC2) + (SD.cutoff*sd(NHW_SAMPLES_EUR$PC2)))
  
  SDSelection <- PCA[PCA$PC1 > PC1min &
                       PCA$PC1 < PC1max &
                       PCA$PC2 > PC2min &
                       PCA$PC2 < PC2max, ]
  
  SDSelection.Table.NHW[[i]] <- as.vector(SDSelection$IID)
  
  p.sd <- p.sd + annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
                          fill=NA, colour="black") +
    annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff.all[i]), size=4, color = "black")
}



## AA
SDSelection.Table.AA <- list()
SD.cutoff.all <- 3:5
for (i in 1:length(SD.cutoff.all)){
  SD.cutoff <- SD.cutoff.all[i]  
  PC1min <- (mean(AA_SAMPLES_AFR$PC1) - (SD.cutoff*sd(AA_SAMPLES_AFR$PC1)))
  PC1max <- (mean(AA_SAMPLES_AFR$PC1) + (SD.cutoff*sd(AA_SAMPLES_AFR$PC1)))
  PC2min <- (mean(AA_SAMPLES_AFR$PC2) - (SD.cutoff*sd(AA_SAMPLES_AFR$PC2)))
  PC2max <- (mean(AA_SAMPLES_AFR$PC2) + (SD.cutoff*sd(AA_SAMPLES_AFR$PC2)))
  
  SDSelection <- PCA[PCA$PC1 > PC1min & 
                       PCA$PC1 < PC1max &
                       PCA$PC2 > PC2min &
                       PCA$PC2 < PC2max,]
  
  SDSelection.Table.AA[[i]] <- as.vector(SDSelection$IID)
  
  p.sd <- p.sd + annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
                          fill=NA, colour="violet") +
    annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff.all[i]), size=4, color = "black")
}



## Asian
SDSelection.Table.Asian <- list()
SD.cutoff.all <- 3:5
for (i in 1:length(SD.cutoff.all)){
  SD.cutoff <- SD.cutoff.all[i]  
  PC1min <- (mean(ASIAN_SAMPLES_EAS$PC1) - (SD.cutoff*sd(ASIAN_SAMPLES_EAS$PC1)))
  PC1max <- (mean(ASIAN_SAMPLES_EAS$PC1) + (SD.cutoff*sd(ASIAN_SAMPLES_EAS$PC1)))
  PC2min <- (mean(ASIAN_SAMPLES_EAS$PC2) - (SD.cutoff*sd(ASIAN_SAMPLES_EAS$PC2)))
  PC2max <- (mean(ASIAN_SAMPLES_EAS$PC2) + (SD.cutoff*sd(ASIAN_SAMPLES_EAS$PC2)))
  
  SDSelection <- PCA[PCA$PC1 > PC1min & 
                       PCA$PC1 < PC1max &
                       PCA$PC2 > PC2min &
                       PCA$PC2 < PC2max,]
  
  SDSelection.Table.Asian[[i]] <- as.vector(SDSelection$IID)
  
  p.sd <- p.sd + annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
                          fill=NA, colour="green") +
    annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff.all[i]), size=4, color = "red")
}

p.sd


##################################
## Extract samples based on PCA ##
##################################
## See the regions to be selected from each Ethnicity
# p.sd.reportedAFRICAN + annotate("rect", xmin=-0.015, xmax=-0.0025, ymin=-0.008, ymax=-0.0014,
#                                 fill=NA, colour="blue") +
#   annotate("text", x=PC1max, y=PC2max, label=paste0("selected: ",SD.cutoff.all[i]), size=4, color = "black")

PCA.SJLIFE <- PCA[PCA$Population == "SJLIFE",]

SD.cutoff <- 3 
## AFRICAN
PC1min <- (mean(AA_SAMPLES_AFR$PC1) - (SD.cutoff*sd(AA_SAMPLES_AFR$PC1)))
PC1max <- (mean(AA_SAMPLES_AFR$PC1) + (SD.cutoff*sd(AA_SAMPLES_AFR$PC1)))
PC2min <- (mean(AA_SAMPLES_AFR$PC2) - (SD.cutoff*sd(AA_SAMPLES_AFR$PC2)))
PC2max <- (mean(AA_SAMPLES_AFR$PC2) + (SD.cutoff*sd(AA_SAMPLES_AFR$PC2)))

# For African American, we will keep the subset in the shaded region within the blue rectangle
SELECTED.p <- p.sd + annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
                                              colour="blue", alpha = .3) 
PCA.AFR <- PCA.SJLIFE[PCA.SJLIFE$PC1 > PC1min &
                        PCA.SJLIFE$PC1 < PC1max &
                        PCA.SJLIFE$PC2 > PC2min &
                        PCA.SJLIFE$PC2 < PC2max,]
PCA.AFR$Ethnicity <- "AFR"

ALL.PCA.AFR  <- PCA[PCA$PC1 > PC1min &
                      PCA$PC1 < PC1max &
                      PCA$PC2 > PC2min &
                      PCA$PC2 < PC2max,]



## ASIAN
PC1min <- (mean(ASIAN_SAMPLES_EAS$PC1) - (SD.cutoff*sd(ASIAN_SAMPLES_EAS$PC1)))
PC1max <- (mean(ASIAN_SAMPLES_EAS$PC1) + (SD.cutoff*sd(ASIAN_SAMPLES_EAS$PC1)))
PC2min <- (mean(ASIAN_SAMPLES_EAS$PC2) - (SD.cutoff*sd(ASIAN_SAMPLES_EAS$PC2)))
PC2max <- (mean(ASIAN_SAMPLES_EAS$PC2) + (SD.cutoff*sd(ASIAN_SAMPLES_EAS$PC2)))

# For Asian, we will keep the subset in the shaded green square because the center of SD is a little bit displaced from our center
SELECTED.p <- SELECTED.p + annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
                                    colour="blue", alpha = .3) 


PCA.EAS <- PCA.SJLIFE[PCA.SJLIFE$PC1 > PC1min &
                        PCA.SJLIFE$PC1 < PC1max &
                        PCA.SJLIFE$PC2 > PC2min &
                        PCA.SJLIFE$PC2 < PC2max,]

PCA.EAS$Ethnicity <- "EAS"

ALL.PCA.EAS  <- PCA[PCA$PC1 > PC1min &
             PCA$PC1 < PC1max &
             PCA$PC2 > PC2min &
             PCA$PC2 < PC2max,]


## NHW
PC1min <- (mean(NHW_SAMPLES_EUR$PC1) - (SD.cutoff*sd(NHW_SAMPLES_EUR$PC1)))
PC1max <- (mean(NHW_SAMPLES_EUR$PC1) + (SD.cutoff*sd(NHW_SAMPLES_EUR$PC1)))
PC2min <- (mean(NHW_SAMPLES_EUR$PC2) - (SD.cutoff*sd(NHW_SAMPLES_EUR$PC2)))
PC2max <- (mean(NHW_SAMPLES_EUR$PC2) + (SD.cutoff*sd(NHW_SAMPLES_EUR$PC2)))

# For EUR, we will keep samples within 5SD
SELECTED.p <- SELECTED.p + annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
                                    colour="blue", alpha = .3) 

PCA.EUR <- PCA.SJLIFE[PCA.SJLIFE$PC1 > PC1min &
                        PCA.SJLIFE$PC1 < PC1max &
                        PCA.SJLIFE$PC2 > PC2min &
                        PCA.SJLIFE$PC2 < PC2max,]

PCA.EUR$Ethnicity <- "EUR"

ALL.PCA.EUR  <- PCA[PCA$PC1 > PC1min &
                      PCA$PC1 < PC1max &
                      PCA$PC2 > PC2min &
                      PCA$PC2 < PC2max,]


SELECTED.p

ggsave("SJLIFE-ALL-COHORT-3-5SD-plot.jpg", plot = SELECTED.p, device = NULL, scale = 1, width = 10, height = 8, dpi = 600, limitsize = TRUE)

##############################################
## Annotate SJLIFE samples by PCA ethnicity ##
##############################################
# read admixture results file
admixture.df <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/MERGED_SJLIFE_PLINK_PER_CHR/final.5.Q_header2_SJLIFE_only", header = T )

# Recode ethnicities based on 80% or hiher cutoff; if nor column value is greater than 80%, insert MIX
admixture.df$Ethnicity.80perc <- names(admixture.df)[-1][max.col(admixture.df[-1] > 0.8)]
admixture.df$Ethnicity.80perc[!rowSums(admixture.df[2:6] > 0.8)] <- "MIX"



## PCA ethnicities based on 3SD cutoff
PCA.Ethnicities <- rbind.data.frame(PCA.EUR, PCA.AFR, PCA.EAS)
table(PCA.Ethnicities$Ethnicity)
# AFR  EAS  EUR 
# 645   16 3439

admixture.df$Ethnicity.PCA <- PCA.Ethnicities$Ethnicity[match(admixture.df$INDIVIDUAL, PCA.Ethnicities$IID)]

# Check how many are true
table(admixture.df$Ethnicity.80perc == admixture.df$Ethnicity.PCA)
# FALSE  TRUE 
# 81  4019

admixture.df$TRUE_FALSE <- admixture.df$Ethnicity.80perc == admixture.df$Ethnicity.PCA
table(admixture.df$Ethnicity.PCA, admixture.df$TRUE_FALSE)
#       FALSE TRUE
# AFR    53  592
# EAS     0   16
# EUR    28 3411

write.table(admixture.df, "SJLIFE_4481_Admixture_PCA_ethnicity.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)

## Samples to include in a second round of PCA
tt <- rbind.data.frame(ALL.PCA.AFR, ALL.PCA.EAS, ALL.PCA.EUR)
write.table(tt[1:2], "samples.to.exclude.round2.pca.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)

#########################
## Cleaned PCA round 1 ##
#########################
PCA <- read.table("final_cleaned1-PCAS.eigenvec", header =T, stringsAsFactors=FALSE)
dim(PCA)
thousandG.ethnicty <- read.table("integrated_call_samples_v3.20130502.ALL.panel", header = T , stringsAsFactors = FALSE)
head(thousandG.ethnicty)


PCA$Population <- "SJLIFE"
PCA$Population <- thousandG.ethnicty$super_pop[match(PCA$IID, thousandG.ethnicty$sample)]
PCA <- PCA[c("FID", "IID", c(paste0("PC", 1:10), "Population"))]
PCA$KEY <- paste(PCA$FID, PCA$IID, sep =":")
PCA$Population <- as.character(PCA$Population)
PCA$Population[is.na(PCA$Population)] <- "SJLIFE"
# write.table(PCA, "SJLIFE_PCA_with_sample_population.txt", sep ="\t", col.names = T, quote = F)

# ## Change the target to reorder the dots on top of the PCA plot
# target <- c("AFR", "AMR", "EAS", "EUR", "SAS", "SJLIFE")
# PCA$Population <- factor(PCA$Population, levels = target)
# PCA <- PCA[order(-as.numeric(factor(PCA$Population))),]

library(ggplot2)
## Generate a new file that has IID, PC1,PC2, and a new column Population 
p <- ggplot(PCA, aes(x=PC1, y=PC2, color=Population)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("SJLIFE_4481") +
  scale_color_manual(values = c(AFR='violet', AMR='blue', EAS='red', EUR="black", SAS="pink", SJLIFE="brown")) +
  annotate("text", x=0.02, y=-0.01, label="", size=4, color = "violet") +
  annotate("text", x=0.005, y=0.024, label="", size=4, color = "blue") +
  annotate("text", x=0.011, y=0.033, label="", size=4, color = "red") +
  annotate("text", x=-0.008, y=0.004, label="", size=4, color = "black") +
  annotate("text", x=0.0072, y=0.013, label="", size=4, color = "pink") +
  annotate("text", x=-0.01, y=0, label="", size=4, color = "brown") +
  theme_classic() +
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face = "plain", size=12),
        axis.text=element_text(face = "plain"),axis.title = element_text(face = "plain"),plot.title = element_text(face = "plain", hjust = 0.5,size=13))

p


ggsave("SJLIFE-ALL-Population_CLEANED_round1.jpg", plot = p, device = NULL, scale = 1, width = 10, height = 8, dpi = 600, limitsize = TRUE)



