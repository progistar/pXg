

library(ggplot2)
library(RColorBrewer)
library("dplyr")
library(readxl)
library(latex2exp)
library("stringr")

setwd("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features_rank")
feat1res <- read.csv("_all.rank10.BA", header =T, sep="\t", as.is = as.double())

staticThemeTop <- theme(text = element_text(size=25, color = "black")
                        , axis.text.x = element_text(size=20, color="black"), 
                        axis.text.y = element_text(size=20, color="black"),
                        legend.justification = c("center"),
                        legend.position= "none",
                        legend.text = element_text(size=20, color = "black"))

calculate_fdr <- function (data, sample, BAFilter = F) {
  data <- data[data$Sample == sample, ]
  data <- data[order(data$percolator_score, decreasing = TRUE), ]
  method <- "All"
  if(BAFilter) {
    data <- data[data$EL_Rank < 2 , ]
    method <- "Sub"
  }
  fdrs <- seq(0.01, 0.1, 0.01)
  #fdrs <- c(0.01, 0.1, 0.01)
  
  qData <- data.frame(matrix(nrow=0, ncol=9))
  colnames(qData) <- c("FDR","pFDR", "Score", "PSM", "Binder PSM", "Peptide", "Binder peptide", "Sample", "Method")
  
  for(fdr in fdrs) {
    
    if (fdr != 0.05) {
      next
    }
    
    fdr_score <- 100
    cur_fdr <- 0
    target <- 0
    decoy <- 0
    
    print(fdr)
    for(idx in seq(1,nrow(data),1)) {
      record <- data[idx,]
      if(record$Label == 1) {
        target <- target + 1
        if(decoy/target < fdr) {
          fdr_score <- record$percolator_score
          cur_fdr <- decoy/target
        }
      } else {
        decoy <- decoy + 1
      }
    }
    
    psms <- data[data$percolator_score >= fdr_score & data$Label == 1, ]
    peptides <- psms[!duplicated(psms$InferredPeptide),]
    
    qData[nrow(qData)+1, ] <- c(fdr,
                                cur_fdr,
                                fdr_score,
                                nrow(psms),
                                nrow(psms[psms$EL_Rank < 2, ]),
                                nrow(peptides),
                                nrow(peptides[peptides$EL_Rank < 2, ]),
                                sample, method)
  }
  
  qData$FDR <- as.double(qData$FDR)
  qData$Score <- as.double(qData$Score)
  qData$`Binder PSM` <- as.double(qData$`Binder PSM`)
  qData$Peptide <- as.double(qData$Peptide)
  qData$`Binder peptide` <- as.double(qData$`Binder peptide`)
  
  return(qData)
}

calculate_td_ratio <- function (data, sample, BAFilter = F) {
  data <- data[data$Sample == sample, ]
  data <- data[order(data$percolator_score, decreasing = TRUE), ]
  method <- "All"
  if(BAFilter) {
    data <- data[data$EL_Rank < 2 , ]
    method <- "Sub"
  }
  
  qData <- data.frame(matrix(nrow=0, ncol=4))
  colnames(qData) <- c("Target", "Decoy", "Sample", "Method")
  
  all <- nrow(data)
  target <- nrow(data[data$Label == 1, ]) / all
  decoy <- nrow(data[data$Label == -1, ]) / all
  
  qData[nrow(qData)+1,] <- c(target, decoy, sample, method)
  return(qData)
}

#### Target / Decoy Ratio ####
sub_feat1res <- feat1res[feat1res$IsCanonical == "TRUE", ]

BA_Filter = T
C_B_LCL1_TDRatio<- calculate_fdr(sub_feat1res, "B-LCL1", BAFilter = BA_Filter)
C_B_LCL2_TDRatio<- calculate_fdr(sub_feat1res, "B-LCL2", BAFilter = BA_Filter)
C_B_LCL3_TDRatio<- calculate_fdr(sub_feat1res, "B-LCL3", BAFilter = BA_Filter)
C_B_LCL4_TDRatio<- calculate_fdr(sub_feat1res, "B-LCL4", BAFilter = BA_Filter)
C_DOHH2_TDRatio<- calculate_fdr(sub_feat1res, "DOHH2", BAFilter = BA_Filter)
C_HBL1_TDRatio<- calculate_fdr(sub_feat1res, "HBL1", BAFilter = BA_Filter)
C_SUDHL4_TDRatio<- calculate_fdr(sub_feat1res, "SUDHL4", BAFilter = BA_Filter)
C_THP1_1_TDRatio<- calculate_fdr(sub_feat1res, "THP1-1", BAFilter = BA_Filter)
C_THP1_2_TDRatio<- calculate_fdr(sub_feat1res, "THP1-2", BAFilter = BA_Filter)
C_THP1_3_TDRatio<- calculate_fdr(sub_feat1res, "THP1-3", BAFilter = BA_Filter)

C_Data_TDRatio_Sub <- rbind(C_B_LCL1_TDRatio, C_B_LCL2_TDRatio, C_B_LCL3_TDRatio, C_B_LCL4_TDRatio,
                            C_DOHH2_TDRatio, C_HBL1_TDRatio, C_SUDHL4_TDRatio,
                            C_THP1_1_TDRatio, C_THP1_2_TDRatio, C_THP1_3_TDRatio)

sub_feat1res <- feat1res[feat1res$IsCanonical == "FALSE", ]


BA_Filter = T
C_B_LCL1_TDRatio<- calculate_fdr(sub_feat1res, "B-LCL1", BAFilter = BA_Filter)
C_B_LCL2_TDRatio<- calculate_fdr(sub_feat1res, "B-LCL2", BAFilter = BA_Filter)
C_B_LCL3_TDRatio<- calculate_fdr(sub_feat1res, "B-LCL3", BAFilter = BA_Filter)
C_B_LCL4_TDRatio<- calculate_fdr(sub_feat1res, "B-LCL4", BAFilter = BA_Filter)
C_DOHH2_TDRatio<- calculate_fdr(sub_feat1res, "DOHH2", BAFilter = BA_Filter)
C_HBL1_TDRatio<- calculate_fdr(sub_feat1res, "HBL1", BAFilter = BA_Filter)
C_SUDHL4_TDRatio<- calculate_fdr(sub_feat1res, "SUDHL4", BAFilter = BA_Filter)
C_THP1_1_TDRatio<- calculate_fdr(sub_feat1res, "THP1-1", BAFilter = BA_Filter)
C_THP1_2_TDRatio<- calculate_fdr(sub_feat1res, "THP1-2", BAFilter = BA_Filter)
C_THP1_3_TDRatio<- calculate_fdr(sub_feat1res, "THP1-3", BAFilter = BA_Filter)

C_Data_TDRatio_Sub <- rbind(C_Data_TDRatio_Sub, C_B_LCL1_TDRatio, C_B_LCL2_TDRatio, C_B_LCL3_TDRatio, C_B_LCL4_TDRatio,
                            C_DOHH2_TDRatio, C_HBL1_TDRatio, C_SUDHL4_TDRatio,
                            C_THP1_1_TDRatio, C_THP1_2_TDRatio, C_THP1_3_TDRatio)
C_Data_TDRatio_Sub

fdrs <- read_excel(path = "/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/3.WhyRank10/IDRanks.xlsx", 
                   sheet = "Sheet1")

scatterFDRPlot <- ggplot(data=fdrs[fdrs$Class != "Canonical",], 
                         aes(x=Rank, y=`Binder peptide`, color=Sample)) +
  theme_bw() +
  geom_point(size=3, alpha = 0.8) + 
  staticThemeTop +
  theme(text = element_text(size=25), 
        strip.background = element_blank()) +
  facet_grid(rows = vars(`Sample`), scales = "free") +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10)) +
  labs(y= "", x = "Rank")

scatterFDRPlot

ggsave("Rank_Noncanonical.png", plot = scatterFDRPlot, width = 4, height = 15, units = "in", dpi = 600)
















