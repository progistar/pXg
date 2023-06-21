
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

setwd("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy")
feat1res <- read_excel(path = "BAAnalysis.xlsx", sheet = "Feat2_selected")

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
sub_feat1res <- feat1res[feat1res$IsCanonical == T, ]
BA_Filter = T
C_B_LCL1_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL1", BAFilter = BA_Filter)
C_B_LCL2_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL2", BAFilter = BA_Filter)
C_B_LCL3_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL3", BAFilter = BA_Filter)
C_B_LCL4_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL4", BAFilter = BA_Filter)
C_DOHH2_TDRatio<- calculate_td_ratio(sub_feat1res, "DOHH2", BAFilter = BA_Filter)
C_HBL1_TDRatio<- calculate_td_ratio(sub_feat1res, "HBL1", BAFilter = BA_Filter)
C_SUDHL4_TDRatio<- calculate_td_ratio(sub_feat1res, "SUDHL4", BAFilter = BA_Filter)
C_THP1_1_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-1", BAFilter = BA_Filter)
C_THP1_2_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-2", BAFilter = BA_Filter)
C_THP1_3_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-3", BAFilter = BA_Filter)

C_Data_TDRatio_Sub <- rbind(C_B_LCL1_TDRatio, C_B_LCL2_TDRatio, C_B_LCL3_TDRatio, C_B_LCL4_TDRatio,
                        C_DOHH2_TDRatio, C_HBL1_TDRatio, C_SUDHL4_TDRatio,
                        C_THP1_1_TDRatio, C_THP1_2_TDRatio, C_THP1_3_TDRatio)

BA_Filter = F
C_B_LCL1_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL1", BAFilter = BA_Filter)
C_B_LCL2_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL2", BAFilter = BA_Filter)
C_B_LCL3_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL3", BAFilter = BA_Filter)
C_B_LCL4_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL4", BAFilter = BA_Filter)
C_DOHH2_TDRatio<- calculate_td_ratio(sub_feat1res, "DOHH2", BAFilter = BA_Filter)
C_HBL1_TDRatio<- calculate_td_ratio(sub_feat1res, "HBL1", BAFilter = BA_Filter)
C_SUDHL4_TDRatio<- calculate_td_ratio(sub_feat1res, "SUDHL4", BAFilter = BA_Filter)
C_THP1_1_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-1", BAFilter = BA_Filter)
C_THP1_2_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-2", BAFilter = BA_Filter)
C_THP1_3_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-3", BAFilter = BA_Filter)

C_Data_TDRatio_All <- rbind(C_B_LCL1_TDRatio, C_B_LCL2_TDRatio, C_B_LCL3_TDRatio, C_B_LCL4_TDRatio,
                            C_DOHH2_TDRatio, C_HBL1_TDRatio, C_SUDHL4_TDRatio,
                            C_THP1_1_TDRatio, C_THP1_2_TDRatio, C_THP1_3_TDRatio)

sub_feat1res <- feat1res[feat1res$IsCanonical == F, ]
BA_Filter = T
C_B_LCL1_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL1", BAFilter = BA_Filter)
C_B_LCL2_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL2", BAFilter = BA_Filter)
C_B_LCL3_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL3", BAFilter = BA_Filter)
C_B_LCL4_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL4", BAFilter = BA_Filter)
C_DOHH2_TDRatio<- calculate_td_ratio(sub_feat1res, "DOHH2", BAFilter = BA_Filter)
C_HBL1_TDRatio<- calculate_td_ratio(sub_feat1res, "HBL1", BAFilter = BA_Filter)
C_SUDHL4_TDRatio<- calculate_td_ratio(sub_feat1res, "SUDHL4", BAFilter = BA_Filter)
C_THP1_1_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-1", BAFilter = BA_Filter)
C_THP1_2_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-2", BAFilter = BA_Filter)
C_THP1_3_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-3", BAFilter = BA_Filter)

NC_Data_TDRatio_Sub <- rbind(C_B_LCL1_TDRatio, C_B_LCL2_TDRatio, C_B_LCL3_TDRatio, C_B_LCL4_TDRatio,
                            C_DOHH2_TDRatio, C_HBL1_TDRatio, C_SUDHL4_TDRatio,
                            C_THP1_1_TDRatio, C_THP1_2_TDRatio, C_THP1_3_TDRatio)

BA_Filter = F
C_B_LCL1_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL1", BAFilter = BA_Filter)
C_B_LCL2_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL2", BAFilter = BA_Filter)
C_B_LCL3_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL3", BAFilter = BA_Filter)
C_B_LCL4_TDRatio<- calculate_td_ratio(sub_feat1res, "B-LCL4", BAFilter = BA_Filter)
C_DOHH2_TDRatio<- calculate_td_ratio(sub_feat1res, "DOHH2", BAFilter = BA_Filter)
C_HBL1_TDRatio<- calculate_td_ratio(sub_feat1res, "HBL1", BAFilter = BA_Filter)
C_SUDHL4_TDRatio<- calculate_td_ratio(sub_feat1res, "SUDHL4", BAFilter = BA_Filter)
C_THP1_1_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-1", BAFilter = BA_Filter)
C_THP1_2_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-2", BAFilter = BA_Filter)
C_THP1_3_TDRatio<- calculate_td_ratio(sub_feat1res, "THP1-3", BAFilter = BA_Filter)

NC_Data_TDRatio_All <- rbind(C_B_LCL1_TDRatio, C_B_LCL2_TDRatio, C_B_LCL3_TDRatio, C_B_LCL4_TDRatio,
                            C_DOHH2_TDRatio, C_HBL1_TDRatio, C_SUDHL4_TDRatio,
                            C_THP1_1_TDRatio, C_THP1_2_TDRatio, C_THP1_3_TDRatio)

C_Data_TDRatio_Sub$Type <- "Canonical"
C_Data_TDRatio_All$Type <- "Canonical"
NC_Data_TDRatio_Sub$Type <- "Noncanonical"
NC_Data_TDRatio_All$Type <- "Noncanonical"

Data_TDRatio <- rbind(C_Data_TDRatio_All, C_Data_TDRatio_Sub, NC_Data_TDRatio_All, NC_Data_TDRatio_Sub)
Data_TDRatio$Target <- as.double(Data_TDRatio$Target)
Data_TDRatio$Decoy <- as.double(Data_TDRatio$Decoy)

sd(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "All", ]$Decoy)
mean(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "All", ]$Decoy)

Data_summary <- data.frame(matrix(nrow=0, ncol=6))
colnames(Data_summary) <- c("mTarget", "mDecoy", "sdTarget", "sdDecoy", "Method", "Type")
Data_summary[1, ] <- c(mean(Data_TDRatio[Data_TDRatio$Type == "Canonical" & Data_TDRatio$Method == "All", ]$Target),
                       mean(Data_TDRatio[Data_TDRatio$Type == "Canonical" & Data_TDRatio$Method == "All", ]$Decoy),
                       sd(Data_TDRatio[Data_TDRatio$Type == "Canonical" & Data_TDRatio$Method == "All", ]$Target),
                       sd(Data_TDRatio[Data_TDRatio$Type == "Canonical" & Data_TDRatio$Method == "All", ]$Decoy),
                       "All", "Canonical")

Data_summary[2, ] <- c(mean(Data_TDRatio[Data_TDRatio$Type == "Canonical" & Data_TDRatio$Method == "Sub", ]$Target),
                       mean(Data_TDRatio[Data_TDRatio$Type == "Canonical" & Data_TDRatio$Method == "Sub", ]$Decoy),
                       sd(Data_TDRatio[Data_TDRatio$Type == "Canonical" & Data_TDRatio$Method == "Sub", ]$Target),
                       sd(Data_TDRatio[Data_TDRatio$Type == "Canonical" & Data_TDRatio$Method == "Sub", ]$Decoy),
                       "Sub", "Canonical")

Data_summary[3, ] <- c(mean(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "All", ]$Target),
                       mean(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "All", ]$Decoy),
                       sd(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "All", ]$Target),
                       sd(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "All", ]$Decoy),
                       "All", "Noncanonical")

Data_summary[4, ] <- c(mean(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "Sub", ]$Target),
                       mean(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "Sub", ]$Decoy),
                       sd(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "Sub", ]$Target),
                       sd(Data_TDRatio[Data_TDRatio$Type == "Noncanonical" & Data_TDRatio$Method == "Sub", ]$Decoy),
                       "Sub", "Noncanonical")

Data_summary$mTarget <- as.double(Data_summary$mTarget)
Data_summary$sdTarget <- as.double(Data_summary$sdTarget)
Data_summary
#### Target Decoy Ratio #####3
display.brewer.pal(n = 10, name = 'Set3')
tdRatio <- ggplot(data=Data_summary, aes(x=Method, y=mTarget)) +
  theme_bw() +
  geom_bar(stat="identity", fill = brewer.pal(n = 10, name = "Set3")[5]) +
  geom_errorbar( aes(x=Method, ymin=mTarget-sdTarget, ymax=mTarget+sdTarget), width=0.4, colour=brewer.pal(n = 10, name = "Set3")[4], alpha=0.9, size=1.5) +
  scale_y_continuous(n.breaks = 3) +
  staticThemeRightTop +
  labs(y= "", x = "") +
  facet_grid(cols = vars(Type)) +
  theme(text = element_text(size=25), strip.background = element_blank(), legend.position = "none",
        plot.margin = margin(0,0.1,0,-0.3, "in")) +
  coord_cartesian(ylim = c(0.5, 1.0)) 
tdRatio

ggsave("Figure3_TD_Ratio.png", plot = tdRatio, width = 4.5, height = 4, units = "in", dpi = 600)


#sub_feat1res <- feat1res[feat1res$percolator_score > 0 & feat1res$IsCanonical == T, ]
sub_feat1res <- feat1res[feat1res$IsCanonical == T, ]
BA_Filter = T
NC_B_LCL1_Sub<- calculate_fdr(sub_feat1res, "B-LCL1", BAFilter = BA_Filter)
NC_B_LCL2_Sub<- calculate_fdr(sub_feat1res, "B-LCL2", BAFilter = BA_Filter)
NC_B_LCL3_Sub<- calculate_fdr(sub_feat1res, "B-LCL3", BAFilter = BA_Filter)
NC_B_LCL4_Sub<- calculate_fdr(sub_feat1res, "B-LCL4", BAFilter = BA_Filter)
NC_DOHH2_Sub<- calculate_fdr(sub_feat1res, "DOHH2", BAFilter = BA_Filter)
NC_HBL1_Sub<- calculate_fdr(sub_feat1res, "HBL1", BAFilter = BA_Filter)
NC_SUDHL4_Sub<- calculate_fdr(sub_feat1res, "SUDHL4", BAFilter = BA_Filter)
NC_THP1_1_Sub<- calculate_fdr(sub_feat1res, "THP1-1", BAFilter = BA_Filter)
NC_THP1_2_Sub<- calculate_fdr(sub_feat1res, "THP1-2", BAFilter = BA_Filter)
NC_THP1_3_Sub<- calculate_fdr(sub_feat1res, "THP1-3", BAFilter = BA_Filter)

BA_Filter <- F
NC_B_LCL1_All<- calculate_fdr(sub_feat1res, "B-LCL1", BAFilter = BA_Filter)
NC_B_LCL2_All<- calculate_fdr(sub_feat1res, "B-LCL2", BAFilter = BA_Filter)
NC_B_LCL3_All<- calculate_fdr(sub_feat1res, "B-LCL3", BAFilter = BA_Filter)
NC_B_LCL4_All<- calculate_fdr(sub_feat1res, "B-LCL4", BAFilter = BA_Filter)
NC_DOHH2_All<- calculate_fdr(sub_feat1res, "DOHH2", BAFilter = BA_Filter)
NC_HBL1_All<- calculate_fdr(sub_feat1res, "HBL1", BAFilter = BA_Filter)
NC_SUDHL4_All<- calculate_fdr(sub_feat1res, "SUDHL4", BAFilter = BA_Filter)
NC_THP1_1_All<- calculate_fdr(sub_feat1res, "THP1-1", BAFilter = BA_Filter)
NC_THP1_2_All<- calculate_fdr(sub_feat1res, "THP1-2", BAFilter = BA_Filter)
NC_THP1_3_All<- calculate_fdr(sub_feat1res, "THP1-3", BAFilter = BA_Filter)


Data_Sub <- rbind(NC_B_LCL1_Sub, NC_B_LCL2_Sub, NC_B_LCL3_Sub, NC_B_LCL4_Sub, 
                  NC_DOHH2_Sub, NC_HBL1_Sub, NC_SUDHL4_Sub, 
                  NC_THP1_1_Sub, NC_THP1_2_Sub, NC_THP1_3_Sub)
Data_All <- rbind(NC_B_LCL1_All, NC_B_LCL2_All, NC_B_LCL3_All, NC_B_LCL4_All, 
                  NC_DOHH2_All, NC_HBL1_All, NC_SUDHL4_All, 
                  NC_THP1_1_All, NC_THP1_2_All, NC_THP1_3_All)

Data <- rbind(Data_Sub, Data_All)

Data_Sub[Data_Sub$FDR == 0.05 & Data_Sub$Method == "Sub", ]

#### noncanonical ncPSMs according to ncFDR #####3
scatterFDRPlot <- ggplot(data=Data, aes(x=FDR, y=`Binder PSM`, group=Method, color=Method, shape = Method)) +
  theme_bw() +
  geom_point(size=3, aes(fill=Method, shape = Method), alpha = 0.8) + 
  theme(text = element_text(size=25), strip.background = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.margin = margin(0,0,-0.2,0, "in"),
        plot.margin = margin(0,0.1,-0.25,-0.3, "in")) +
  staticThemeTop +
  facet_grid(cols = vars(`Sample`), scales = "free") +
  labs(y= "", x = "")

scatterFDRPlot
ggsave("Figure3_ncMAP_FDR.png", plot = scatterFDRPlot, width = 14, height = 4, units = "in", dpi = 600)

# HLB1 = 79
# DOHH2 = 60
# SHUDHL4 = 51
# THP1-1 = 138
# THP1-2 = 99
# THP1-3 = 194

### 5%? ###
sa_plot <- function(res, data, sample_name, score) {
  nc_score <- 0
  print(nc_score)
  c_data <- res[res$Sample == sample_name & res$EL_Rank < 101 & res$percolator_score >= score, ]
  
  c_data <- c_data[order(c_data$percolator_score, decreasing = TRUE), ]
  c_data <- c_data[!duplicated(c_data$InferredPeptide),]
  
  #nc_data <- c_data[c_data$IsCanonical == F & c_data$percolator_score > nc_score & c_data$Label == 1, ]$percolator_score
  #nc_decoy_data <- c_data[c_data$IsCanonical == F & c_data$percolator_score > nc_score & c_data$Label == -1, ]$percolator_score
  #c_decoy_data <- c_data[c_data$IsCanonical == T & c_data$percolator_score > nc_score & c_data$Label == -1, ]$percolator_score
  #c_data <- c_data[c_data$IsCanonical == T, ]$percolator_score
  
  nc_data <- c_data[c_data$IsCanonical == F & c_data$Label == 1, ]$`ALC (%)`
  nc_decoy_data <- c_data[c_data$IsCanonical == F & c_data$Label == -1, ]$`ALC (%)`
  c_decoy_data <- c_data[c_data$IsCanonical == T & c_data$Label == -1, ]$`ALC (%)`
  c_data <- c_data[c_data$IsCanonical == T, ]$`ALC (%)`
  
  print(length(nc_data)+length(c_data))
  print((length(nc_decoy_data)+length(c_decoy_data))/(length(nc_data)+length(c_data)))
  
  return(boxplot(c_data, nc_data, c_decoy_data,nc_decoy_data))
}
feat1res$`ALC (%)`
feat1res$Log2Reads
feat1res$MeanQScore

Data_Sub[Data_Sub$FDR == 0.05 & Data_Sub$Method == "Sub", ]
sa_plot(feat1res, Data, "B-LCL1", -0.08636710)
sa_plot(feat1res, Data, "B-LCL2", -0.11245800)
sa_plot(feat1res, Data, "B-LCL3", 0.23778500)
sa_plot(feat1res, Data, "B-LCL4", 0.25685900)
sa_plot(feat1res, Data, "DOHH2", 0.44799300)
sa_plot(feat1res, Data, "HBL1", 0.00791873)
sa_plot(feat1res, Data, "SUDHL4", 0.56563800)
sa_plot(feat1res, Data, "THP1-1", 0.00102963)
sa_plot(feat1res, Data, "THP1-2", -0.06280690)
sa_plot(feat1res, Data, "THP1-3", -0.03567030)



##### The number of PSMs per peptide
fdr_5_res_canonical <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "Canonical")
fdr_5_res_noncanonical <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "Noncanonical")
display.brewer.pal(n = 10, name = 'Set3')
brewer.pal(n = 10, name = "Set3")[5]

fdr_5_res <- rbind(fdr_5_res_canonical, fdr_5_res_noncanonical)
fdr_5_res$`ALC (%)`
fdr_5_res$Log2MeanQScore
generate_box_plot <- function(data) {
  data$Class <- "Canonical"
  data[data$IsCanonical == F, ]$Class <- "Nonanonical"
  mdPlot <- ggplot(data=data, aes(x=Class, y=SA, fill=Class)) +
    theme_bw() +
    scale_fill_manual(values=c(brewer.pal(n = 10, name = "Set3")[5], brewer.pal(n = 10, name = "Set3")[4])) +
    geom_violin() +
    geom_boxplot(width=.1) +
    theme(text = element_text(size=20)) +
    scale_y_continuous(n.breaks = 4) +
    ylim(0, max(data$SA)) +
    staticThemeRightTop +
    theme(axis.text.x = element_blank(), legend.title = element_blank()) +
    labs(y= "", x="")
    #stat_compare_means(method = "t.test", label = "p.format", label.x = 1.25, size = 8, hide.ns = T, label.y=5.3)
  return (mdPlot)
}
wilcox.test(fdr_5_res_canonical$SA, fdr_5_res_noncanonical$SA)
wilcox.test(fdr_5_res_canonical$`ALC (%)`, fdr_5_res_noncanonical$`ALC (%)`)

t.test(fdr_5_res_canonical$SA, fdr_5_res_noncanonical$SA)
t.test(fdr_5_res_canonical$`ALC (%)`, fdr_5_res_noncanonical$`ALC (%)`)
t.test(fdr_5_res_canonical$MeanQScore, fdr_5_res_noncanonical$MeanQScore)

alc_box <- generate_box_plot(fdr_5_res)
sa_box <- generate_box_plot(fdr_5_res)
q_box <- generate_box_plot(fdr_5_res)
sa_box
alc_box
q_box
boxes <- ggarrange(alc_box, sa_box, q_box,
                     ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
boxes
ggsave("Figure3_PSM_Qual.png", plot = boxes, width = 12, height = 5, units = "in", dpi = 600)

#### Length distribution ####
fdr_5_peptide <- fdr_5_res
fdr_5_peptide$key <- paste(fdr_5_res$InferredPeptide, fdr_5_res$Sample, sep = "_")
fdr_5_peptide <- fdr_5_peptide[!duplicated(fdr_5_peptide$key), ]

length_data <- data.frame(matrix(nrow=0, ncol=3))
colnames(length_data) <- c("Length", "Prop", "Class")

for (idx in seq(8, 15, 1)) {
  length_data[nrow(length_data)+1, ] <- c(idx, nrow(fdr_5_peptide[fdr_5_peptide$IsCanonical == F & fdr_5_peptide$Length == idx, ]) / nrow(fdr_5_peptide[fdr_5_peptide$IsCanonical == F, ]), 
                        "Noncanonical")  
}

for (idx in seq(8, 15, 1)) {
  length_data[nrow(length_data)+1, ] <- c(idx, nrow(fdr_5_peptide[fdr_5_peptide$IsCanonical == T & fdr_5_peptide$Length == idx, ]) / nrow(fdr_5_peptide[fdr_5_peptide$IsCanonical == T, ]), 
                                          "Canonical")  
}
length_data$Prop <- as.double(length_data$Prop)
length_data$Length <- as.double(length_data$Length)

f_length = 15
f_matrix <- matrix(c(nrow(fdr_5_peptide[fdr_5_peptide$IsCanonical==F & fdr_5_peptide$Length == f_length, ]),
                     nrow(fdr_5_peptide[fdr_5_peptide$IsCanonical==T & fdr_5_peptide$Length == f_length, ]),
                     nrow(fdr_5_peptide[fdr_5_peptide$IsCanonical==F & fdr_5_peptide$Length != f_length, ]),
                     nrow(fdr_5_peptide[fdr_5_peptide$IsCanonical==T & fdr_5_peptide$Length != f_length, ])), nrow=2, ncol=2, byrow = F)
fisher.test(f_matrix, alternative = "two.sided")

staticThemeRightTopInner <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), 
                             axis.text.y = element_text(size=20, color="black"),
                             legend.justification = c("right", "top"),
                             legend.position= c(.98, .99),
                             legend.text = element_text(size=20, color = "black"))

lengthPlot <- ggplot(data=length_data, aes(x=Length, y=Prop, fill = Class)) +
  scale_fill_manual(values=c(brewer.pal(n = 10, name = "Set3")[5], brewer.pal(n = 10, name = "Set3")[4])) +
  geom_bar(position = position_dodge2(), stat = "identity") +
  scale_x_continuous(n.breaks = 8) +
  ylim(0, 0.8) +
  theme_bw() +
  labs(y="", x = "") +
  theme(plot.margin = margin(0.1,0.2,-0.2,-0.2, "in"), legend.title=element_blank()) +
  staticThemeRightTopInner
lengthPlot
ggsave("Figure3_Length.png", plot = lengthPlot, width = 8, height = 6, units = "in", dpi = 600)

### get FDR record###
get_fdr_record <- function (data, sample, BAFilter = F, fdr) {
  data <- data[data$Sample == sample, ]
  data <- data[order(data$percolator_score, decreasing = TRUE), ]
  method <- "All"
  if(BAFilter) {
    data <- data[data$EL_Rank < 2 , ]
    method <- "Sub"
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
  return(as.data.frame(psms))
}

sub_feat1res <- feat1res[feat1res$IsCanonical == T, ]
BA_Filter = T
NC_B_LCL1_Sub<- get_fdr_record(sub_feat1res, "B-LCL1", BAFilter = BA_Filter, 0.05)
NC_B_LCL2_Sub<- get_fdr_record(sub_feat1res, "B-LCL2", BAFilter = BA_Filter, 0.05)
NC_B_LCL3_Sub<- get_fdr_record(sub_feat1res, "B-LCL3", BAFilter = BA_Filter, 0.05)
NC_B_LCL4_Sub<- get_fdr_record(sub_feat1res, "B-LCL4", BAFilter = BA_Filter, 0.05)
NC_DOHH2_Sub<- get_fdr_record(sub_feat1res, "DOHH2", BAFilter = BA_Filter, 0.05)
NC_HBL1_Sub<- get_fdr_record(sub_feat1res, "HBL1", BAFilter = BA_Filter, 0.05)
NC_SUDHL4_Sub<- get_fdr_record(sub_feat1res, "SUDHL4", BAFilter = BA_Filter, 0.05)
NC_THP1_1_Sub<- get_fdr_record(sub_feat1res, "THP1-1", BAFilter = BA_Filter, 0.05)
NC_THP1_2_Sub<- get_fdr_record(sub_feat1res, "THP1-2", BAFilter = BA_Filter, 0.05)
NC_THP1_3_Sub<- get_fdr_record(sub_feat1res, "THP1-3", BAFilter = BA_Filter, 0.05)

Data_Sub <- rbind(NC_B_LCL1_Sub, NC_B_LCL2_Sub, NC_B_LCL3_Sub, NC_B_LCL4_Sub, 
                  NC_DOHH2_Sub, NC_HBL1_Sub, NC_SUDHL4_Sub, 
                  NC_THP1_1_Sub, NC_THP1_2_Sub, NC_THP1_3_Sub)

write.table(Data_Sub, file = "canonical_fdr5.tsv", append = FALSE, quote = F, sep = "\t",
            na = "NA", dec = ".", row.names = F,
            col.names = TRUE)

