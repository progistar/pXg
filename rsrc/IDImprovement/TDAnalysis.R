
library(ggplot2)
library(RColorBrewer)
library("dplyr")
library(readxl)
library(latex2exp)
library("stringr")
library(ggpmisc)
library(ComplexHeatmap)
library(ggpubr)
library(plyr)
library(heatmaply)
library(enrichplot)

## COLOR SETTING
brewer.pal(9, "Set1")
redC <- "#E41A1C" ## FOR Novel things
blueC <- "#377EB8" ## FOR Protein coding things
greenC <- "#4DAF4A" ## FOR ELSE
purpleC <- "#984EA3" ##
orangeC <- "#FF7F00"
yellowC <- "#FFFF33"
brownC <- "#A65628"
pinkC <- "#F781BF"
grayC <- "#999999"

sampleCa <- c(`Subject 1` = redC, `Subject 2`=blueC, `Subject 3`=greenC, `Subject 4`=purpleC)
sampleC <- c(redC, blueC, greenC, purpleC)

staticThemeLeftTop <- theme(text = element_text(size=25, color = "black")
                            , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                            legend.justification = c("right", "top"),
                            legend.position= c(.22, .98),
                            legend.text = element_text(size=20, color = "black"),
                            legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeRightTop <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                             legend.justification = c("right", "top"),
                             legend.text = element_text(size=20, color = "black"))

staticThemeRightBottom <- theme(text = element_text(size=25, color = "black")
                                , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                                legend.justification = c("right", "bottom"),
                                legend.position= c(.98, .05),
                                legend.text = element_text(size=20, color = "black"),
                                legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeRight <- theme(text = element_text(size=25, color = "black")
                          , axis.text.x = element_text(size=20, color="black"), 
                          axis.text.y = element_text(size=20, color="black"),
                          legend.justification = c("right"),
                          legend.text = element_text(size=20, color = "black"))


staticThemeTop <- theme(text = element_text(size=25, color = "black")
                        , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                        legend.justification = c("center"),
                        legend.position= "top",
                        legend.text = element_text(size=20, color = "black"),
                        legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeBottom <- theme(text = element_text(size=25, color = "black")
                           , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                           legend.justification = c("center"),
                           legend.position= "bottom",
                           legend.text = element_text(size=20, color = "black"),
                           legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeNone <- theme(text = element_text(size=25, color = "black")
                         , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                         legend.position = "none")

setwd("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy")

## ALC score distribution

td_score <- read_excel(path = "TargetDecoyProperties.xlsx", sheet = "DIFF_nocut")

td_score_both <- td_score[td_score$Category == "Both but decoy win" | td_score$Category == "Both but target win", ]
rtPlot <- ggplot(data = td_score_both, aes(x=MaxTargetScore, y=MaxDecoyScore)) +
  theme_bw() +
  geom_point(size=1, alpha = 0.6) +
  staticThemeRight +
  theme(plot.margin = margin(0, 0.5, 0, 0, "in")) +
  labs(y= "ALC score in decoy", x = NULL) +
  labs(x= "ALC score in target", x = NULL) +
  geom_segment(x=0, y=0, xend=100, yend=100) +
  facet_grid(cols = vars(`Sample`), rows = vars(`Category`), scales = "free")
rtPlot
ggsave("both_nocut_score_scatter.png", plot = rtPlot, width = 15, height = 7, units = "in", dpi = 1200)

td_score_decoy <- td_score[td_score$Category == "Both but decoy win" | td_score$Category == "Decoy only", ]
breakers <- seq(0, 100, 5)

tdPlot <- ggplot(data=td_score_decoy, aes(x=MaxDecoyScore, fill=Category)) +
  scale_fill_manual(values=c(blueC, redC)) +
  geom_histogram(alpha = 0.5, position = "identity", breaks = breakers) +
  theme_bw() +
  labs(y="PSM", x = "Average local confidence (ALC)") +
  facet_grid(rows = vars(`Sample`), scales = "free") +
  staticThemeRightTop

tdPlot
ggsave("both_nocut_score_hist.png", plot = tdPlot, width = 15, height = 7, units = "in", dpi = 1200)

## Delta ALC score distribution

td_score <- read_excel(path = "TargetDecoyProperties.xlsx", sheet = "DIFF_nocut")

td_score_both <- td_score[td_score$Category == "Both but target win", ]
td_score_both <- td_score_both[td_score_both$MaxTargetScore >= 80, ]
td_score_both$Category <- "80 ~ 84"
td_score_both[td_score_both$MaxTargetScore >= 85 & td_score_both$MaxTargetScore < 90,]$Category <- "85 ~ 89"
td_score_both[td_score_both$MaxTargetScore >= 90,]$Category <- "90 ~"
max(td_score_both$MaxTargetScore - td_score_both$MaxDecoyScore)
breakers <- seq(0, 50, 5)
tdPlot <- ggplot(data=td_score_both, aes(x=MaxTargetScore - MaxDecoyScore, fill=Category)) +
  scale_fill_manual(values=c(blueC, redC, greenC)) +
  geom_histogram(alpha = 0.5, position = "dodge", breaks = breakers) +
  theme_bw() +
  labs(y="PSM", x = "Delta ALC") +
  facet_grid(rows = vars(`Sample`), cols = vars(`Category`), scales = "free") +
  staticThemeRightTop
tdPlot
ggsave("both_nocut_target_score_diff_hist.png", plot = tdPlot, width = 15, height = 7, units = "in", dpi = 1200)

tdPlot <- ggplot(data=td_score_both, aes(x=MaxTargetRNA - MaxDecoyRNA, fill=Category)) +
  scale_fill_manual(values=c(blueC, redC, greenC)) +
  geom_histogram(alpha = 0.5, position = "dodge", binwidth = 1) +
  theme_bw() +
  labs(y="PSM", x = "Delta log₂(reads + 1)") +
  facet_grid(rows = vars(`Sample`), cols = vars(`Category`), scales = "free") +
  staticThemeRightTop
tdPlot
ggsave("both_nocut_target_rna_diff_hist.png", plot = tdPlot, width = 15, height = 7, units = "in", dpi = 1200)


## RNA Analysis

td_score <- read_excel(path = "TargetDecoyProperties.xlsx", sheet = "DIFF_nocut")

td_score_both <- td_score[td_score$Category == "Both but decoy win" | td_score$Category == "Both but target win", ]

rtPlot <- ggplot(data = td_score_both, aes(x=MaxDecoyScore, y=MaxDecoyRNA)) +
  theme_bw() +
  geom_point(size=1, alpha = 0.6) +
  staticThemeRight +
  theme(plot.margin = margin(0, 0.5, 0, 0, "in")) +
  labs(y= "ALC score in decoy", x = NULL) +
  labs(x= "ALC score in target", x = NULL) +
  geom_segment(x=0, y=0, xend=100, yend=100) +
  facet_grid(cols = vars(`Sample`), rows = vars(`Category`), scales = "free")
rtPlot
ggsave("both_nocut_rna_scatter.png", plot = rtPlot, width = 15, height = 7, units = "in", dpi = 1200)

td_score_decoy <- td_score[td_score$Category == "Both but decoy win" | td_score$Category == "Decoy only", ]

tdPlot <- ggplot(data=td_score_decoy, aes(x=MaxDecoyRNA, fill=Category)) +
  scale_fill_manual(values=c(blueC, redC)) +
  geom_histogram(alpha = 0.5, position = "identity", binwidth = 2) +
  theme_bw() +
  labs(y="PSM", x = "Average local confidence (ALC)") +
  facet_grid(rows = vars(`Sample`), scales = "free") +
  staticThemeRightTop

tdPlot
ggsave("both_nocut_rna_hist.png", plot = tdPlot, width = 15, height = 7, units = "in", dpi = 1200)


########## TD Distribution in Percolator #############

label_preprocess <- function(target, decoy, sample_name) {
  dataS_target <- read.csv(file = target, header = T, sep="\t", as.is = as.double())
  dataS_decoy <- read.csv(file = decoy, header = T, sep="\t", as.is = as.double())
  dataS_target <- data.frame(dataS_target)
  dataS_decoy <- data.frame(dataS_decoy)
  dataS_target$Category <- "Target"
  dataS_decoy$Category <- "Decoy"
  dataS <- rbind(dataS_target, dataS_decoy)
  dataS$Sample <- sample_name
  
  dataS
}

feat1res <- read_excel(path = "BAAnalysis.xlsx", sheet = "Feat2_selected")
feat1res$Category <- "Decoy"
feat1res[feat1res$Label == 1, ]$Category <- "Target"

feat1res$Category <- factor(x = feat1res$Category, levels = c("Target", "Decoy"))


tdPlot <- ggplot(data=feat1res[feat1res$IsCanonical == T, ], aes(x=percolator_score, fill = Category)) +
  scale_fill_manual(values=c(blueC, redC)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 100) +
  #xlim(-minimax, minimax) +
  theme_bw() +
  labs(y="PSM", x = "Score") +
  facet_grid(rows = vars(`Sample`), scales = "free") +
  staticThemeRightTop

tdPlot

ggsave("TD_graph_feat5.png", plot = tdPlot, width = 15, height = 20, units = "in", dpi = 600)

########## BA analysis ##########
#S1_p001 <- S1_p001[, -c(44,45,46,47,48,49)]
#S2_p001 <- S2_p001[, -c(44,45,46,47,48,49)]
#S3_p001 <- S3_p001[, -c(44,45,46,47,48,49)]
#S4_p001 <- S4_p001[, -c(44,45,46,47,48,49)]

feat1res <- read_excel(path = "BAAnalysis.xlsx", sheet = "Feat1")

dohh2 <- read_excel(path = "BAAnalysis_10samples.xlsx", sheet = "DOHH2.cuevas")
hbl1 <- read_excel(path = "BAAnalysis_10samples.xlsx", sheet = "HBL1.cuevas")
sudhl4 <- read_excel(path = "BAAnalysis_10samples.xlsx", sheet = "SUDHL4.cuevas")

thp1 <- read_excel(path = "BAAnalysis_10samples.xlsx", sheet = "THP1_1.scull")
thp2 <- read_excel(path = "BAAnalysis_10samples.xlsx", sheet = "THP1_2.scull")
thp3 <- read_excel(path = "BAAnalysis_10samples.xlsx", sheet = "THP1_3.scull")

feat1res$Label <- as.character(feat1res$Label)

#SA
feat1res$percolator_score
feat1res$BestDeltaRT
feat1res$SA
feat1res$mLog2BestELRank
feat1res$Log2Reads
feat1res$Log2MeanQScore
feat1res$`ALC (%)`

mapPlot <- ggplot(data=feat1res, aes(x=Label, y=Log2Reads, fill=Label)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  staticThemeRightTop +
  labs(y= "Log2MeanQScore", x = "Peptide length") +
  facet_grid(cols = vars(IsCanonical), rows = vars(Sample))

mapPlot

ggsave("tmp.png", plot = mapPlot, width = 15, height = 15, units = "in", dpi = 600)

scull_res_cal <- function (res) {
  res$BA <- "Binder"
  res[res$`MHC-I` == "NB", ]$BA <- "Non-binder"
  
  return(res)
}

cuevas_res_cal <- function (res) {
  new_data  <- res[res$isdecoy == F, ]
  new_data$BA <- "Binder"
  new_data[new_data$`MHC-I` == "NB",]$BA <- "Non-binder"
  new_data$ppm <- 1000000 * (new_data$experimentalmasstocharge-new_data$calculatedmasstocharge)/new_data$calculatedmasstocharge
  return(new_data)
}
dohh2 <- cuevas_res_cal(dohh2)
hbl1 <- cuevas_res_cal(hbl1)
sudhl4 <- cuevas_res_cal(sudhl4)
thp1 <- scull_res_cal(thp1)
thp2 <- scull_res_cal(thp2)
thp3 <- scull_res_cal(thp3)

pep_cal <- function (res) {
  res$BA <- "Non-binder"
  res[res$EL_Rank < 2, ]$BA <- "Binder"
  
  qvals <- seq(0.01, 1, 0.01)
  qData <- data.frame(matrix(nrow=0, ncol=10))
  colnames(qData) <- c("PEP", "Sample", "binder-ratioPSM","binder-ratioPeptide", "MAP", "peptide", "binderPSM","PSM","class", "MedianRNA")
  sample_list <- c("B-LCL1","B-LCL2","B-LCL3","B-LCL4",
                   "THP1-1","THP1-2","THP1-3",
                   "DOHH2","HBL1","SUDHL4")
  for(qval in qvals) {
    #subtmp <- res[res$`pep` <= qval, ]
    subtmp <- res[res$`q-value` <= qval, ]
    #subtmp <- subtmp[subtmp$Length == 8 | subtmp$Length == 9 | subtmp$Length == 10 | subtmp$Length == 11, ]
    subtmp <- subtmp[subtmp$IsCanonical == F, ]
    
    for(sample_name in sample_list) {
      sample <- subtmp[subtmp$Sample == sample_name, ]
      peptide <- sample[!duplicated(sample$InferredPeptide),]
      
      nBinder_psm <- nrow(sample[sample$BA == "Binder", ])
      nPSM <- nrow(sample)
      binder_ratio_psm <- nBinder_psm/nPSM
      
      nMAP <- nrow(peptide[peptide$BA == "Binder",])
      nPeptide <- nrow(peptide)
      binder_ratio_peptide <- nMAP/nPeptide
      
      log2Reads <- log2(median(peptide$Reads)+1)
      
      
      if(is.nan(binder_ratio_psm)) {
        binder_ratio_psm <- NA
      }
      qData[nrow(qData)+1, ] <- c(qval, sample_name, binder_ratio_psm, binder_ratio_peptide, 
                                  nMAP,nPeptide, nBinder_psm,nPSM,"noncanonical", log2Reads)
    }
    
    #subtmp <- res[res$`pep` <= qval, ]
    subtmp <- res[res$`q-value` <= qval, ]
    #subtmp <- subtmp[subtmp$Length == 8 | subtmp$Length == 9 | subtmp$Length == 10 | subtmp$Length == 11, ]
    subtmp <- subtmp[subtmp$IsCanonical == T, ]
    
    for(sample_name in sample_list) {
      sample <- subtmp[subtmp$Sample == sample_name, ]
      peptide <- sample[!duplicated(sample$InferredPeptide),]
      
      nBinder_psm <- nrow(sample[sample$BA == "Binder", ])
      nPSM <- nrow(sample)
      binder_ratio_psm <- nBinder_psm/nPSM
      
      nMAP <- nrow(peptide[peptide$BA == "Binder",])
      nPeptide <- nrow(peptide)
      binder_ratio_peptide <- nMAP/nPeptide
      
      log2Reads <- log2(median(peptide$Reads)+1)
      if(is.nan(binder_ratio_psm)) {
        binder_ratio_psm <- NA
      }
      qData[nrow(qData)+1, ] <- c(qval, sample_name, binder_ratio_psm, binder_ratio_peptide, 
                                  nMAP,nPeptide, nBinder_psm,nPSM,"canonical", log2Reads)
    }
    
  }
  
  qData$`PEP` <- as.double(qData$`PEP`)
  qData$MAP <- as.double(qData$MAP)
  qData$peptide <- as.double(qData$peptide)
  qData$`binder-ratioPeptide` <- as.double(qData$`binder-ratioPeptide`)
  qData$`binder-ratioPSM` <- as.double(qData$`binder-ratioPSM`)
  qData$binderPSM <- as.double(qData$binderPSM)
  
  return(qData)
}


feat4Data <- pep_cal(feat4res)
feat5Data <- pep_cal(feat5res)
feat8Data <- pep_cal(feat8res)
feat9Data <- pep_cal(feat9res)

feat4Data$Feature <- "Feature 4"
feat5Data$Feature <- "Feature 5"
feat8Data$Feature <- "Feature 8"
feat9Data$Feature <- "Feature 9"
feat45 <- rbind(feat4Data, feat5Data, feat8Data, feat9Data)

feat4Data[feat4Data$PEP == 0.01, ]
feat5Data[feat5Data$PEP == 0.01, ]
feat8Data[feat8Data$PEP == 0.01, ]
feat9Data[feat9Data$PEP == 0.1, ]

# HLB1 = 79
# DOHH2 = 60
# SHUDHL4 = 51
# THP1-1 = 138
# THP1-2 = 99
# THP1-3 = 194
feat4Data$`binder-ratioPeptide`
linePlot <- ggplot(data=feat45[feat45$PEP <= 0.2 & feat45$class == "noncanonical", ], aes(x=`PEP`, y=`binder-ratioPeptide`, group=Sample, color=Sample)) +
  theme_bw() +
  #scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set3") +
  geom_line() +
  geom_point() +
  theme(text = element_text(size=20)) +
  staticThemeRightTop +
  labs(y= "Binder ratio", x = "PEP") +
  facet_grid(cols = vars(Feature), scales = "free")

linePlot

ggsave("PEP_Plot_Canonical_02_Peptide.png", plot = linePlot, width = 15, height = 7, units = "in", dpi = 600)

g <- ggplot(data = feat45[feat45$PEP == 0.01 & feat45$class == "noncanonical", ], aes(x=PEP, y=`binder-ratioPeptide`, fill=Feature)) +
  theme_bw() +
  geom_bar(position= position_dodge2(),stat = "identity")+
  scale_color_brewer(palette="Set3") +
  xlab("") +
  ylab("MAP") +
  staticThemeRightTop +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in"), axis.text.x = element_blank()) +
  facet_grid(cols = vars(Sample), scales = "free")
g

ggsave("MAP_Plot_Canonical_02.png", plot = g, width = 15, height = 7, units = "in", dpi = 600)

linePlot <- ggplot(data=dohh2[dohh2$NewType == "canonical", ], aes(x=`peaks:peptidescore`, y=`binder-ratioPeptide`, group=Sample, color=Sample)) +
  theme_bw() +
  #scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set3") +
  geom_line() +
  geom_point() +
  theme(text = element_text(size=20)) +
  staticThemeRightTop +
  labs(y= "Binder ratio", x = "PEP") +
  facet_grid(cols = vars(Feature), scales = "free")

linePlot

## AI대학원 3,000 정밀 1,000 혁신허브 3,300
## PCC => PEAKS 1,700, Desktop upgrade
## 

Cuevas_DOHH2_NC_Binderatio <- nrow(dohh2[dohh2$NewType == "noncanonical" & dohh2$BA == "Binder", ]) / nrow(dohh2[dohh2$NewType == "noncanonical", ])
Cuevas_DOHH2_NC_MAP <- nrow(dohh2[dohh2$NewType == "noncanonical" & dohh2$BA == "Binder", ])

Cuevas_HBL1_NC_Binderatio <- nrow(hbl1[hbl1$NewType == "noncanonical" & hbl1$BA == "Binder", ]) / nrow(hbl1[hbl1$NewType == "noncanonical", ])
Cuevas_HBL1_NC_MAP <- nrow(hbl1[hbl1$NewType == "noncanonical" & hbl1$BA == "Binder", ])

Cuevas_SUDHL4_NC_Binderatio <- nrow(sudhl4[sudhl4$NewType == "noncanonical" & sudhl4$BA == "Binder", ]) / nrow(sudhl4[sudhl4$NewType == "noncanonical", ])
Cuevas_SUDHL4_NC_MAP <- nrow(sudhl4[sudhl4$NewType == "noncanonical" & sudhl4$BA == "Binder", ])

Cuevas_HBL1_NC_MAP
Cuevas_SUDHL4_NC_MAP
THP1_1_NC_Binderatio <- 138 / 339
THP1_1_NC_MAP <- 138

THP1_2_NC_Binderatio <- 99 / 248
THP1_2_NC_MAP <- 99

THP1_3_NC_Binderatio <- 138 / 516
THP1_3_NC_MAP <- 194


Cuevas_DOHH2_C_Binderatio <- nrow(dohh2[dohh2$NewType == "canonical" & dohh2$BA == "Binder", ]) / nrow(dohh2[dohh2$NewType == "canonical", ])
Cuevas_DOHH2_C_MAP <- nrow(dohh2[dohh2$NewType == "canonical" & dohh2$BA == "Binder", ])

Cuevas_HBL1_C_Binderatio <- nrow(hbl1[hbl1$NewType == "canonical" & hbl1$BA == "Binder", ]) / nrow(hbl1[hbl1$NewType == "canonical", ])
Cuevas_HBL1_C_MAP <- nrow(hbl1[hbl1$NewType == "canonical" & hbl1$BA == "Binder", ])

Cuevas_SUDHL4_C_Binderatio <- nrow(sudhl4[sudhl4$NewType == "canonical" & sudhl4$BA == "Binder", ]) / nrow(sudhl4[sudhl4$NewType == "canonical", ])
Cuevas_SUDHL4_C_MAP <- nrow(sudhl4[sudhl4$NewType == "canonical" & sudhl4$BA == "Binder", ])

Cuevas_DOHH2_C_Binderatio
Cuevas_HBL1_C_Binderatio
Cuevas_SUDHL4_C_Binderatio

THP1_1_C_Binderatio <- 5916 / 7239
THP1_1_C_MAP <- 5916

THP1_2_C_Binderatio <- 4009 / 4559
THP1_2_C_MAP <- 4009

THP1_3_C_Binderatio <- 8637 / 10209
THP1_3_C_MAP <- 8637


linePlot <- ggplot(data=feat5Data[feat5Data$class == "noncanonical",], aes(x=`binder-ratioPeptide`, y=`MAP`, group=Sample, color=Sample)) +
  theme_bw() +
  #scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set3") +
  geom_line() +
  geom_point() +
  scale_x_reverse() +
  theme(text = element_text(size=20)) +
  staticThemeRightTop +
  geom_point(aes(x=Cuevas_DOHH2_NC_Binderatio, y=Cuevas_DOHH2_NC_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=Cuevas_HBL1_NC_Binderatio, y=Cuevas_HBL1_NC_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=Cuevas_SUDHL4_NC_Binderatio, y=Cuevas_SUDHL4_NC_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=THP1_1_NC_Binderatio, y=THP1_1_NC_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=THP1_2_NC_Binderatio, y=THP1_2_NC_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=THP1_3_NC_Binderatio, y=THP1_3_NC_MAP), colour="blue", shape = 15) +
  labs(y= "MAP", x = "Binder ratio")

linePlot

ggsave("Compare_Noncanonical.png", plot = linePlot, width = 15, height = 7, units = "in", dpi = 600)

linePlot <- ggplot(data=feat5Data[feat5Data$class == "canonical",], aes(x=`binder-ratioPeptide`, y=`MAP`, group=Sample, color=Sample)) +
  theme_bw() +
  #scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set3") +
  geom_line() +
  geom_point() +
  scale_x_reverse() +
  theme(text = element_text(size=20)) +
  staticThemeRightTop +
  geom_point(aes(x=Cuevas_DOHH2_C_Binderatio, y=Cuevas_DOHH2_C_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=Cuevas_HBL1_C_Binderatio, y=Cuevas_HBL1_C_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=Cuevas_SUDHL4_C_Binderatio, y=Cuevas_SUDHL4_C_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=THP1_1_C_Binderatio, y=THP1_1_C_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=THP1_2_C_Binderatio, y=THP1_2_C_MAP), colour="blue", shape = 15) +
  geom_point(aes(x=THP1_3_C_Binderatio, y=THP1_3_C_MAP), colour="blue", shape = 15) +
  labs(y= "MAP", x = "Binder ratio")

linePlot

ggsave("Compare_Canonical.png", plot = linePlot, width = 15, height = 7, units = "in", dpi = 600)

## Median
Cuevas_DOHH2_NC_Binderatio_Upper <- nrow(dohh2[dohh2$NewType == "noncanonical" & dohh2$BA == "Binder" & dohh2$`peaks:peptidescore` >= median(dohh2$`peaks:peptidescore`), ]) / nrow(dohh2[dohh2$NewType == "noncanonical" & dohh2$`peaks:peptidescore` >= median(dohh2$`peaks:peptidescore`), ])
Cuevas_DOHH2_NC_Binderatio_Lower <- nrow(dohh2[dohh2$NewType == "noncanonical" & dohh2$BA == "Binder" & dohh2$`peaks:peptidescore` < median(dohh2$`peaks:peptidescore`), ]) / nrow(dohh2[dohh2$NewType == "noncanonical" & dohh2$`peaks:peptidescore` < median(dohh2$`peaks:peptidescore`), ])
Cuevas_HBL1_NC_Binderatio_Upper <- nrow(hbl1[hbl1$NewType == "noncanonical" & hbl1$BA == "Binder" & hbl1$`peaks:peptidescore` >= median(hbl1$`peaks:peptidescore`), ]) / nrow(hbl1[hbl1$NewType == "noncanonical" & hbl1$`peaks:peptidescore` >= median(hbl1$`peaks:peptidescore`), ])
Cuevas_HBL1_NC_Binderatio_Lower <- nrow(hbl1[hbl1$NewType == "noncanonical" & hbl1$BA == "Binder" & hbl1$`peaks:peptidescore` < median(hbl1$`peaks:peptidescore`), ]) / nrow(hbl1[hbl1$NewType == "noncanonical" & hbl1$`peaks:peptidescore` < median(hbl1$`peaks:peptidescore`), ])
Cuevas_SUDHL4_NC_Binderatio_Upper <- nrow(sudhl4[sudhl4$NewType == "noncanonical" & sudhl4$BA == "Binder" & sudhl4$`peaks:peptidescore` >= median(sudhl4$`peaks:peptidescore`), ]) / nrow(sudhl4[sudhl4$NewType == "noncanonical" & sudhl4$`peaks:peptidescore` >= median(sudhl4$`peaks:peptidescore`), ])
Cuevas_SUDHL4_NC_Binderatio_Lower <- nrow(sudhl4[sudhl4$NewType == "noncanonical" & sudhl4$BA == "Binder" & sudhl4$`peaks:peptidescore` < median(sudhl4$`peaks:peptidescore`), ]) / nrow(sudhl4[sudhl4$NewType == "noncanonical" & sudhl4$`peaks:peptidescore` < median(sudhl4$`peaks:peptidescore`), ])

Scull_THP1_NC_Binderatio_Upper <- nrow(thp1[thp1$NewType == "noncanonical" & thp1$BA == "Binder" & thp1$`-10lgP` >= median(thp1$`-10lgP`), ]) / nrow(thp1[thp1$NewType == "noncanonical" & thp1$`-10lgP` >= median(thp1$`-10lgP`), ])
Scull_THP1_NC_Binderatio_Lower <- nrow(thp1[thp1$NewType == "noncanonical" & thp1$BA == "Binder" & thp1$`-10lgP` < median(thp1$`-10lgP`), ]) / nrow(thp1[thp1$NewType == "noncanonical" & thp1$`-10lgP` < median(thp1$`-10lgP`), ])
Scull_THP2_NC_Binderatio_Upper <- nrow(thp2[thp2$NewType == "noncanonical" & thp2$BA == "Binder" & thp2$`-10lgP` >= median(thp2$`-10lgP`), ]) / nrow(thp2[thp2$NewType == "noncanonical" & thp2$`-10lgP` >= median(thp2$`-10lgP`), ])
Scull_THP2_NC_Binderatio_Lower <- nrow(thp2[thp2$NewType == "noncanonical" & thp2$BA == "Binder" & thp2$`-10lgP` < median(thp2$`-10lgP`), ]) / nrow(thp2[thp2$NewType == "noncanonical" & thp2$`-10lgP` < median(thp2$`-10lgP`), ])
Scull_THP3_NC_Binderatio_Upper <- nrow(thp3[thp3$NewType == "noncanonical" & thp3$BA == "Binder" & thp3$`-10lgP` >= median(thp3$`-10lgP`), ]) / nrow(thp3[thp3$NewType == "noncanonical" & thp3$`-10lgP` >= median(thp3$`-10lgP`), ])
Scull_THP3_NC_Binderatio_Lower <- nrow(thp3[thp3$NewType == "noncanonical" & thp3$BA == "Binder" & thp3$`-10lgP` < median(thp3$`-10lgP`), ]) / nrow(thp3[thp3$NewType == "noncanonical" & thp3$`-10lgP` < median(thp3$`-10lgP`), ])

Cuevas_DOHH2_C_Binderatio_Upper <- nrow(dohh2[dohh2$NewType == "canonical" & dohh2$BA == "Binder" & dohh2$`peaks:peptidescore` >= median(dohh2$`peaks:peptidescore`), ]) / nrow(dohh2[dohh2$NewType == "canonical" & dohh2$`peaks:peptidescore` >= median(dohh2$`peaks:peptidescore`), ])
Cuevas_DOHH2_C_Binderatio_Lower <- nrow(dohh2[dohh2$NewType == "canonical" & dohh2$BA == "Binder" & dohh2$`peaks:peptidescore` < median(dohh2$`peaks:peptidescore`), ]) / nrow(dohh2[dohh2$NewType == "canonical" & dohh2$`peaks:peptidescore` < median(dohh2$`peaks:peptidescore`), ])
Cuevas_HBL1_C_Binderatio_Upper <- nrow(hbl1[hbl1$NewType == "canonical" & hbl1$BA == "Binder" & hbl1$`peaks:peptidescore` >= median(hbl1$`peaks:peptidescore`), ]) / nrow(hbl1[hbl1$NewType == "canonical" & hbl1$`peaks:peptidescore` >= median(hbl1$`peaks:peptidescore`), ])
Cuevas_HBL1_C_Binderatio_Lower <- nrow(hbl1[hbl1$NewType == "canonical" & hbl1$BA == "Binder" & hbl1$`peaks:peptidescore` < median(hbl1$`peaks:peptidescore`), ]) / nrow(hbl1[hbl1$NewType == "canonical" & hbl1$`peaks:peptidescore` < median(hbl1$`peaks:peptidescore`), ])
Cuevas_SUDHL4_C_Binderatio_Upper <- nrow(sudhl4[sudhl4$NewType == "canonical" & sudhl4$BA == "Binder" & sudhl4$`peaks:peptidescore` >= median(sudhl4$`peaks:peptidescore`), ]) / nrow(sudhl4[sudhl4$NewType == "canonical" & sudhl4$`peaks:peptidescore` >= median(sudhl4$`peaks:peptidescore`), ])
Cuevas_SUDHL4_C_Binderatio_Lower <- nrow(sudhl4[sudhl4$NewType == "canonical" & sudhl4$BA == "Binder" & sudhl4$`peaks:peptidescore` < median(sudhl4$`peaks:peptidescore`), ]) / nrow(sudhl4[sudhl4$NewType == "canonical" & sudhl4$`peaks:peptidescore` < median(sudhl4$`peaks:peptidescore`), ])

Scull_THP1_C_Binderatio_Upper <- nrow(thp1[thp1$NewType == "canonical" & thp1$BA == "Binder" & thp1$`-10lgP` >= median(thp1$`-10lgP`), ]) / nrow(thp1[thp1$NewType == "canonical" & thp1$`-10lgP` >= median(thp1$`-10lgP`), ])
Scull_THP1_C_Binderatio_Lower <- nrow(thp1[thp1$NewType == "canonical" & thp1$BA == "Binder" & thp1$`-10lgP` < median(thp1$`-10lgP`), ]) / nrow(thp1[thp1$NewType == "canonical" & thp1$`-10lgP` < median(thp1$`-10lgP`), ])
Scull_THP2_C_Binderatio_Upper <- nrow(thp2[thp2$NewType == "canonical" & thp2$BA == "Binder" & thp2$`-10lgP` >= median(thp2$`-10lgP`), ]) / nrow(thp2[thp2$NewType == "canonical" & thp2$`-10lgP` >= median(thp2$`-10lgP`), ])
Scull_THP2_C_Binderatio_Lower <- nrow(thp2[thp2$NewType == "canonical" & thp2$BA == "Binder" & thp2$`-10lgP` < median(thp2$`-10lgP`), ]) / nrow(thp2[thp2$NewType == "canonical" & thp2$`-10lgP` < median(thp2$`-10lgP`), ])
Scull_THP3_C_Binderatio_Upper <- nrow(thp3[thp3$NewType == "canonical" & thp3$BA == "Binder" & thp3$`-10lgP` >= median(thp3$`-10lgP`), ]) / nrow(thp3[thp3$NewType == "canonical" & thp3$`-10lgP` >= median(thp3$`-10lgP`), ])
Scull_THP3_C_Binderatio_Lower <- nrow(thp3[thp3$NewType == "canonical" & thp3$BA == "Binder" & thp3$`-10lgP` < median(thp3$`-10lgP`), ]) / nrow(thp3[thp3$NewType == "canonical" & thp3$`-10lgP` < median(thp3$`-10lgP`), ])

feat4resPEP001 <- feat4res[feat4res$pep <= 0.01, ]
pXg_DOHH2 <- feat4resPEP001[feat4resPEP001$Sample == "DOHH2", ]
pXg_DOHH2 <- pXg_DOHH2[order(pXg_DOHH2$percolator_score, decreasing = T), ]
pXg_DOHH2 <- pXg_DOHH2[!duplicated(pXg_DOHH2$InferredPeptide), ]
pXg_HBL1 <- feat4resPEP001[feat4resPEP001$Sample == "HBL1", ]
pXg_HBL1 <- pXg_HBL1[order(pXg_HBL1$percolator_score, decreasing = T), ]
pXg_HBL1 <- pXg_HBL1[!duplicated(pXg_HBL1$InferredPeptide), ]
pXg_SUDHL4 <- feat4resPEP001[feat4resPEP001$Sample == "SUDHL4", ]
pXg_SUDHL4 <- pXg_SUDHL4[order(pXg_SUDHL4$percolator_score, decreasing = T), ]
pXg_SUDHL4 <- pXg_SUDHL4[!duplicated(pXg_SUDHL4$InferredPeptide), ]


boxplot(pXg_DOHH2[pXg_DOHH2$IsCanonical == T, ]$`ALC (%)`, pXg_DOHH2[pXg_DOHH2$IsCanonical == F, ]$`ALC (%)`)
t.test(pXg_DOHH2[pXg_DOHH2$IsCanonical == T, ]$`ALC (%)`, pXg_DOHH2[pXg_DOHH2$IsCanonical == F, ]$`ALC (%)`)

boxplot(log2(pXg_DOHH2[pXg_DOHH2$IsCanonical == T, ]$Reads+1), log2(pXg_DOHH2[pXg_DOHH2$IsCanonical == F, ]$Reads+1))
t.test(log2(pXg_DOHH2[pXg_DOHH2$IsCanonical == T, ]$Reads+1), log2(pXg_DOHH2[pXg_DOHH2$IsCanonical == F, ]$Reads+1))

boxplot(pXg_DOHH2[pXg_DOHH2$IsCanonical == T, ]$percolator_score, pXg_DOHH2[pXg_DOHH2$IsCanonical == F, ]$percolator_score)
t.test(pXg_DOHH2[pXg_DOHH2$IsCanonical == T, ]$percolator_score, pXg_DOHH2[pXg_DOHH2$IsCanonical == F, ]$percolator_score)

boxplot(pXg_HBL1[pXg_HBL1$IsCanonical == T, ]$`ALC (%)`, pXg_HBL1[pXg_HBL1$IsCanonical == F, ]$`ALC (%)`)
t.test(pXg_HBL1[pXg_HBL1$IsCanonical == T, ]$`ALC (%)`, pXg_HBL1[pXg_HBL1$IsCanonical == F, ]$`ALC (%)`)

boxplot(log2(pXg_HBL1[pXg_HBL1$IsCanonical == T, ]$Reads+1), log2(pXg_HBL1[pXg_HBL1$IsCanonical == F, ]$Reads+1))
t.test(log2(pXg_HBL1[pXg_HBL1$IsCanonical == T, ]$Reads+1), log2(pXg_HBL1[pXg_HBL1$IsCanonical == F, ]$Reads+1))

boxplot(pXg_HBL1[pXg_HBL1$IsCanonical == T, ]$`percolator_score`, pXg_HBL1[pXg_HBL1$IsCanonical == F, ]$`percolator_score`)
t.test(pXg_HBL1[pXg_HBL1$IsCanonical == T, ]$`percolator_score`, pXg_HBL1[pXg_HBL1$IsCanonical == F, ]$`percolator_score`)

boxplot(pXg_SUDHL4[pXg_SUDHL4$IsCanonical == T, ]$`ALC (%)`, pXg_SUDHL4[pXg_SUDHL4$IsCanonical == F, ]$`ALC (%)`)
t.test(pXg_SUDHL4[pXg_SUDHL4$IsCanonical == T, ]$`ALC (%)`, pXg_SUDHL4[pXg_SUDHL4$IsCanonical == F, ]$`ALC (%)`)

boxplot(log2(pXg_SUDHL4[pXg_SUDHL4$IsCanonical == T, ]$Reads+1), log2(pXg_SUDHL4[pXg_SUDHL4$IsCanonical == F, ]$Reads+1))
t.test(log2(pXg_SUDHL4[pXg_SUDHL4$IsCanonical == T, ]$Reads+1), log2(pXg_SUDHL4[pXg_SUDHL4$IsCanonical == F, ]$Reads+1))

#####
boxplot(pXg_DOHH2[pXg_DOHH2$IsCanonical == T, ]$ppm, pXg_DOHH2[pXg_DOHH2$IsCanonical == F, ]$ppm)
t.test(pXg_DOHH2[pXg_DOHH2$IsCanonical == T, ]$ppm, pXg_DOHH2[pXg_DOHH2$IsCanonical == F, ]$ppm)


boxplot(dohh2[dohh2$NewType == "canonical", ]$`peaks:peptidescore`, dohh2[dohh2$NewType == "noncanonical", ]$`peaks:peptidescore`)
t.test(dohh2[dohh2$NewType == "canonical", ]$`peaks:peptidescore`, dohh2[dohh2$NewType == "noncanonical", ]$`peaks:peptidescore`)

boxplot(dohh2[dohh2$NewType == "canonical", ]$ppm, dohh2[dohh2$NewType == "noncanonical", ]$ppm)
t.test(dohh2[dohh2$NewType == "canonical", ]$ppm, dohh2[dohh2$NewType == "noncanonical", ]$ppm)


boxplot(dohh2[dohh2$NewType == "canonical" & dohh2$BA == "Binder", ]$`peaks:peptidescore`, dohh2[dohh2$NewType == "noncanonical" & dohh2$BA == "Binder", ]$`peaks:peptidescore`)
t.test(dohh2[dohh2$NewType == "canonical" & dohh2$BA == "Binder", ]$`peaks:peptidescore`, dohh2[dohh2$NewType == "noncanonical" & dohh2$BA == "Binder", ]$`peaks:peptidescore`)


pXg_DOHH2_NC_Binderatio_Upper <- nrow(pXg_DOHH2[pXg_DOHH2$IsCanonical == F & pXg_DOHH2$EL_Rank < 2 & pXg_DOHH2$percolator_score >= median(pXg_DOHH2$percolator_score), ]) / nrow(pXg_DOHH2[pXg_DOHH2$IsCanonical == F & pXg_DOHH2$percolator_score >= median(pXg_DOHH2$percolator_score), ])
pXg_DOHH2_NC_Binderatio_Lower <- nrow(pXg_DOHH2[pXg_DOHH2$IsCanonical == F & pXg_DOHH2$EL_Rank < 2 & pXg_DOHH2$percolator_score < median(pXg_DOHH2$percolator_score), ]) /nrow(pXg_DOHH2[pXg_DOHH2$IsCanonical == F & pXg_DOHH2$percolator_score < median(pXg_DOHH2$percolator_score), ])
pXg_HBL1_NC_Binderatio_Upper <- nrow(pXg_HBL1[pXg_HBL1$IsCanonical == F & pXg_HBL1$EL_Rank < 2 & pXg_HBL1$percolator_score >= median(pXg_HBL1$percolator_score), ]) / nrow(pXg_HBL1[pXg_HBL1$IsCanonical == F & pXg_HBL1$percolator_score >= median(pXg_HBL1$percolator_score), ])
pXg_HBL1_NC_Binderatio_Lower <- nrow(pXg_HBL1[pXg_HBL1$IsCanonical == F & pXg_HBL1$EL_Rank < 2 & pXg_HBL1$percolator_score < median(pXg_HBL1$percolator_score), ]) /nrow(pXg_HBL1[pXg_HBL1$IsCanonical == F & pXg_HBL1$percolator_score < median(pXg_HBL1$percolator_score), ])


pXg_DOHH2_NC_Binderatio_Upper
pXg_DOHH2_NC_Binderatio_Lower
pXg_HBL1_NC_Binderatio_Upper
pXg_HBL1_NC_Binderatio_Lower

Cuevas_DOHH2_NC_Binderatio_Upper
Cuevas_DOHH2_NC_Binderatio_Lower
Cuevas_HBL1_NC_Binderatio_Upper
Cuevas_HBL1_NC_Binderatio_Lower
Cuevas_SUDHL4_NC_Binderatio_Upper
Cuevas_SUDHL4_NC_Binderatio_Lower
Scull_THP1_NC_Binderatio_Upper
Scull_THP1_NC_Binderatio_Lower
Scull_THP2_NC_Binderatio_Upper
Scull_THP2_NC_Binderatio_Lower
Scull_THP3_NC_Binderatio_Upper
Scull_THP3_NC_Binderatio_Lower

Cuevas_DOHH2_C_Binderatio_Upper
Cuevas_DOHH2_C_Binderatio_Lower
Cuevas_HBL1_C_Binderatio_Upper
Cuevas_HBL1_C_Binderatio_Lower
Cuevas_SUDHL4_C_Binderatio_Upper
Cuevas_SUDHL4_C_Binderatio_Lower
Scull_THP1_C_Binderatio_Upper
Scull_THP1_C_Binderatio_Lower
Scull_THP2_C_Binderatio_Upper
Scull_THP2_C_Binderatio_Lower
Scull_THP3_C_Binderatio_Upper
Scull_THP3_C_Binderatio_Lower



some_list <- S_p001[S_p001$IsCanonical == F & S_p001$BA == "Binder", ]

write.table(some_list, file = "ncMAPs.tsv", sep = "\t", row.names = F)

S_p001$BA <- "Non-binder"
S_p001[S_p001$EL_Rank < 2,]$BA <- "Binder"

S_p001$LengthCategory <- 8
S_p001[S_p001$Length == 9, ]$LengthCategory <- 9
S_p001[S_p001$Length == 10, ]$LengthCategory <- 10
S_p001[S_p001$Length == 11, ]$LengthCategory <- 11
S_p001[S_p001$Length > 11, ]$LengthCategory <- 12

S_p001$BA <- factor(x = S_p001$BA, levels = c("Binder", "Non-binder"))
S_p001$Length <- as.character(S_p001$Length)
S_p001$Length <- factor(S_p001$Length, levels=c("8", "9", "10", "11", "12", "13", "14", "15"))

subtmp <- S_p001
subtmp <- S_p001[S_p001$`q-value` < 0.01, ]
subtmp <- subtmp[subtmp$Length == 8 | subtmp$Length == 9 | subtmp$Length == 10 | subtmp$Length == 11, ]

subtmp <- subtmp[order(subtmp$percolator_score, decreasing = TRUE), ]
subtmp <- subtmp[!duplicated(subtmp[c("Sample", "InferredPeptide")]), ]


mapPlot <- ggplot(data=subtmp, aes(x=Length, y=log(BestScore+1.5, 2), fill=Length)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  scale_y_continuous(breaks=seq(from=0, to=7, by= 1), limits = c(0,7)) +
  staticThemeNone +
  labs(y= "Log₂(1.5+EL_Rank)", x = "Peptide length") +
  facet_grid(cols = vars(Sample), rows = vars(IsCanonical)) +
  geom_hline(yintercept=log2(1.5+2.0), linetype="dashed", color = "blue", size=0.5) +
  geom_hline(yintercept=log2(1.5+0.5), linetype="dashed", color = "red", size=0.5)

mapPlot

ggsave("FDR_p001_BADist.png", plot = mapPlot, width = 15, height = 7, units = "in", dpi = 600)



#### mzID converter
BiocManager::install("mzID")

library(mzID)

mzResults <- mzID("/Users/gistar/projects/pXg/PreviousStudyResource/THP1_3rdExp_BB72_Filter_TR3peptides_1_1_0.mzid")
mzResults
mzResultsDT <- flatten(mzResults)
mzResultsDT$peplen <- nchar(mzResultsDT$pepseq)
write.table(mzResultsDT[, -c(17)], file = "THP1_3_Filter3.csv", quote = T, row.names = F)

nonPassID <- mzResultsDT[mzResultsDT$passthreshold == F & mzResultsDT$peplen >= 8 & mzResultsDT$peplen <= 11 & mzResultsDT$rank == 1, ]
passID <- mzResultsDT[mzResultsDT$passthreshold == T & mzResultsDT$peplen >= 8 & mzResultsDT$peplen <= 11 & mzResultsDT$rank == 1, ]


