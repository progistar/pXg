
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
                        legend.text = element_text(size=20, color = "black"))

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

###### Figure 2A: pair count #######
pair_data <- read_excel(path = "PairCount.xlsx", sheet = "PairCount_R")
sum(pair_data[pair_data$`(spectrum, peptide, read)` <= 20, ]$Spectrum) / sum(pair_data[pair_data$`(spectrum, peptide, read)` >= 0, ]$Spectrum)

pair_bar_plot <- ggplot(data = pair_data[pair_data$`(spectrum, peptide, read)` <= 20, ], aes(x=`(spectrum, peptide, read)`, y=Spectrum)) +
  theme_bw() +
  geom_bar(stat = "identity", width = 0.8)+
  ylab("Spectrum") +
  staticThemeNone +
  xlab("Peptide-read pair count") 

pair_bar_plot

ggsave("Figure2_PairDist.png", plot = pair_bar_plot, width = 6, height = 5, units = "in", dpi = 600)


B_LCL1 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL1.predfeat.pxg", header = T, sep="\t", as.is = as.double())
B_LCL2 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL2.predfeat.pxg", header = T, sep="\t", as.is = as.double())
B_LCL3 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL3.predfeat.pxg", header = T, sep="\t", as.is = as.double())
B_LCL4 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL4.predfeat.pxg", header = T, sep="\t", as.is = as.double())
DOHH2 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/DOHH2.predfeat.pxg", header = T, sep="\t", as.is = as.double())
HBL1 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/HBL1.predfeat.pxg", header = T, sep="\t", as.is = as.double())
SUDHL4 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/SUDHL4.predfeat.pxg", header = T, sep="\t", as.is = as.double())
THP1_1 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-1.predfeat.pxg", header = T, sep="\t", as.is = as.double())
THP1_2 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-2.predfeat.pxg", header = T, sep="\t", as.is = as.double())
THP1_3 <- read.csv(file = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-3.predfeat.pxg", header = T, sep="\t", as.is = as.double())

B_LCL1$Sample <- "B-LCL1"
B_LCL2$Sample <- "B-LCL2"
B_LCL3$Sample <- "B-LCL3"
B_LCL4$Sample <- "B-LCL4"
DOHH2$Sample <- "DOHH2"
HBL1$Sample <- "HBL1"
SUDHL4$Sample <- "SUDHL4"
THP1_1$Sample <- "THP1-1"
THP1_2$Sample <- "THP1-2"
THP1_3$Sample <- "THP1-3"


generate_td_plot <- function(data) {
  data$Label <- as.character(data$Label)
  data$ALC.... <- as.double(data$ALC....)
  data <- data[data$IsCanonical == "false", ]
  data[data$Label == '1', ]$Label <- "Target"
  data[data$Label == '-1', ]$Label <- "Decoy"
  data$Label <- factor(x=data$Label, levels = c("Target", "Decoy"))
  
  tdPlot <- ggplot(data=data, aes(x=log2(ALC....), fill = Label)) +
    scale_fill_manual(values=c(blueC, redC)) +
    #geom_histogram(alpha = 0.6, position = "identity", bins = 20) +
    geom_density(alpha = 0.6, position = "identity") +
    scale_x_continuous(n.breaks = 3) +
    theme_bw() +
    labs(y="", x = "") +
    theme(plot.margin = margin(0.1,0.2,0.1,0.2, "in")) +
    staticThemeRightTop
  
  return (tdPlot)
}

B_LCL1_plot <- generate_td_plot(B_LCL1)
B_LCL2_plot <- generate_td_plot(B_LCL2)
B_LCL3_plot <- generate_td_plot(B_LCL3)
B_LCL4_plot <- generate_td_plot(B_LCL4)
DOHH2_plot <- generate_td_plot(DOHH2)
HBL1_plot <- generate_td_plot(HBL1)
SUDHL4_plot <- generate_td_plot(SUDHL4)
THP1_1_plot <- generate_td_plot(THP1_1)
THP1_2_plot <- generate_td_plot(THP1_2)
THP1_3_plot <- generate_td_plot(THP1_3)

td_plot <- ggarrange(B_LCL1_plot, B_LCL2_plot, B_LCL3_plot, B_LCL4_plot, DOHH2_plot,
                     HBL1_plot, SUDHL4_plot, THP1_1_plot, THP1_2_plot, THP1_3_plot,
                         ncol = 5, nrow = 2, common.legend = TRUE, legend = "right")
td_plot
ggsave("Figure2_TD_Log2MeanQScore_Plot.png", plot = td_plot, width = 16, height = 7, units = "in", dpi = 600)

######## ALC by Median #######3
B_LCL1$ALC....
B_LCL1$MeanQScore

display.brewer.pal(n = 8, name = 'Paired')
brewer.pal(n = 8, name = "Paired")

generate_median_plot <- function(data) {
  ## Select target only
  data <- data[data$Label == 1, ]
  data$Expression <- "High"
  median_rna <- median(data$Reads)
  
  data[data$Reads <= median_rna, ]$Expression <- "Low"
  
  data$Expression <- factor(x=data$Expression, levels = c("High", "Low"))
  
  print(median(data[data$Expression == "High", ]$ALC....))
  print(median(data[data$Expression == "Low", ]$ALC....))
  
  mdPlot <- ggplot(data=data, aes(x=Expression, y=ALC...., fill=Expression)) +
    theme_bw() +
    scale_fill_manual(values=c("#1F78B4", "#A6CEE3")) +
    geom_boxplot() +
    theme(text = element_text(size=20)) +
    scale_y_continuous(breaks = c(25,50,75,100)) +
    staticThemeRightTop +
    theme(axis.text.x = element_blank()) +
    labs(y= "", x="") +
    stat_compare_means(method = "t.test", label = "p.signif", label.x = 1.2, size = 12, hide.ns = T, label.y=90)
  return (mdPlot)
}

generate_median_qual_plot <- function(data) {
  ## Select target only
  data <- data[data$Label == -1, ]
  data$`Phred score` <- "High"
  median_phred <- median(data$MeanQScore)
  quarter <- quantile(data$MeanQScore, prob = 0.75, type = 1)
  quater <- mean(data$MeanQScore)
  
  data[data$MeanQScore <= quarter, ]$`Phred score` <- "Low"
  
  data$`Phred score` <- factor(x=data$`Phred score`, levels = c("High", "Low"))
  
  print(median(data[data$`Phred score` == "High", ]$ALC....))
  print(median(data[data$`Phred score` == "Low", ]$ALC....))
  
  mdPlot <- ggplot(data=data, aes(x=`Phred score`, y=ALC...., fill=`Phred score`)) +
    theme_bw() +
    scale_fill_manual(values=c("#1F78B4", "#A6CEE3")) +
    geom_boxplot() +
    theme(text = element_text(size=20)) +
    scale_y_continuous(breaks = c(25,50,75,100)) +
    staticThemeRightTop +
    theme(axis.text.x = element_blank()) +
    labs(y= "", x="") +
    stat_compare_means(method = "t.test", label = "p.format", label.x = 1.2, size = 10, hide.ns = T, label.y=90)
  return (mdPlot)
}

B_LCL1_md_plot
tmp_median <- median(B_LCL2[B_LCL2$Label == 1, ]$MeanQScore)
t.test(B_LCL2[B_LCL2$Label ==1 & B_LCL2$MeanQScore < tmp_median, ]$ALC....,
       B_LCL2[B_LCL2$Label ==1 & B_LCL2$MeanQScore >= tmp_median, ]$ALC....)


B_LCL1_md_plot <- generate_median_qual_plot(B_LCL1) + theme( plot.margin = margin(0,0,0,-2))
B_LCL2_md_plot <- generate_median_qual_plot(B_LCL2) + theme(axis.text.y = element_blank(), plot.margin = margin(0,0,0,-0.5))
B_LCL3_md_plot <- generate_median_qual_plot(B_LCL3) + theme(axis.text.y = element_blank(), plot.margin = margin(0,0,0,-0.5))
B_LCL4_md_plot <- generate_median_qual_plot(B_LCL4) + theme(axis.text.y = element_blank(), plot.margin = margin(0,0,0,-0.5))
DOHH2_md_plot <- generate_median_qual_plot(DOHH2) + theme(axis.text.y = element_blank(), plot.margin = margin(0,0,0,-0.5))
HBL1_md_plot <- generate_median_qual_plot(HBL1) + theme(axis.text.y = element_blank(), plot.margin = margin(0,0,0,-0.5))
SUDHL4_md_plot <- generate_median_qual_plot(SUDHL4) + theme(axis.text.y = element_blank(), plot.margin = margin(0,0,0,-0.5))
THP1_1_md_plot <- generate_median_qual_plot(THP1_1) + theme(axis.text.y = element_blank(), plot.margin = margin(0,0,0,-0.5))
THP1_2_md_plot <- generate_median_qual_plot(THP1_2) + theme(axis.text.y = element_blank(), plot.margin = margin(0,0,0,-0.5))
THP1_3_md_plot <- generate_median_qual_plot(THP1_3) + theme(axis.text.y = element_blank(), plot.margin = margin(0,0,0,-0.5))

md_plot <- ggarrange(B_LCL1_md_plot, B_LCL2_md_plot, B_LCL3_md_plot, B_LCL4_md_plot, DOHH2_md_plot,
                     HBL1_md_plot, SUDHL4_md_plot, THP1_1_md_plot, THP1_2_md_plot, THP1_3_md_plot,
                     ncol = 10, nrow = 1, common.legend = TRUE, legend = "right")
md_plot
ggsave("Figure2_MD_Read_ALC_Plot_Decoy.png", plot = md_plot, width = 22, height = 4, units = "in", dpi = 600)


#### Scatter plot of rna high/low #####
pair_data <- read_excel(path = "MedianScatterPlot.xlsx", sheet = "Phred")
pair_data$Label <- factor(pair_data$Label, levels = c("Target", "Decoy"))
scatterPhredPlot <- ggplot(data=pair_data, aes(x=Low, y=High, group=Sample, color=Sample)) +
  theme_bw() +
  scale_color_brewer(palette = "Set3") +
  geom_point(color='black', shape=21, size=3, aes(fill=Sample)) + 
  scale_fill_brewer(palette="Set3") +
  theme(text = element_text(size=20), strip.background = element_blank()) +
  scale_x_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_y_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  staticThemeRightTop +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", linewidth=0.5) +
  facet_grid(cols = vars(`Label`), scales = "free") +
  labs(y= "", x = "")

pair_data <- read_excel(path = "MedianScatterPlot.xlsx", sheet = "RNA_Expression")
pair_data$Label <- factor(pair_data$Label, levels = c("Target", "Decoy"))
scatterRNAPlot <- ggplot(data=pair_data, aes(x=Low, y=High, group=Sample, color=Sample)) +
  theme_bw() +
  scale_color_brewer(palette = "Set3") +
  geom_point(color='black', shape=21, size=3, aes(fill=Sample)) + 
  scale_fill_brewer(palette="Set3") +
  theme(text = element_text(size=20), strip.background = element_blank()) +
  scale_x_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  scale_y_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
  staticThemeRightTop +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", linewidth=0.5) +
  facet_grid(cols = vars(`Label`), scales = "free") +
  labs(y= "", x = "")

pair_plot <- ggarrange(scatterRNAPlot, scatterPhredPlot,
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")

ggsave("Figure2_MD_ALC_Scatter_Plot.png", plot = pair_plot, width = 16, height = 5, units = "in", dpi = 600)

### Target and decoy counts


### Feature plots
classCa <- c(blueC, redC)
### Supplementary Figure
Data <- rbind(B_LCL1, B_LCL2, B_LCL3, B_LCL4, DOHH2, HBL1, THP1_1, THP1_2, THP1_3)
Data[Data$Label == 1, ]$Label <- "Target"
Data[Data$Label == -1, ]$Label <- "Decoy"
Data$Label <- factor(Data$Label, levels = c("Target", "Decoy"))
SAPlot <- ggplot(data=Data, aes(x=Label, y=mLog2BestELRank, fill=Label)) +
  theme_bw() +
  scale_fill_manual(values = classCa) +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  staticThemeRightTop +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text.x = element_blank()) +
  labs(y= "SA", x="") +
  facet_grid(cols = vars(Sample))

SAPlot

ggsave("Figure2_TD_SA_Plot.png", plot = SAPlot, width = 15, height = 7, units = "in", dpi = 600)

######### Normalized weights ##########
### Figure 2B
weights_data <- read_excel(path = "PairCount.xlsx", sheet = "Weight")

weights_data[weights_data$Feature == "Log2Reads", ]$Feature <- "RNA exp."
weights_data[weights_data$Feature == "Log2MeanPhredScore", ]$Feature <- "Phred score"
weights_data[weights_data$Feature == "Absppm", ]$Feature <- "Abs(ppm)"
weights_data[weights_data$Feature == "BestDeltaRT", ]$Feature <- "Best delta RT"
weights_data[weights_data$Feature == "Charge1", ]$Feature <- "Charge 1"
weights_data[weights_data$Feature == "Charge2", ]$Feature <- "Charge 2"
weights_data[weights_data$Feature == "Charge3", ]$Feature <- "Charge 3"
weights_data[weights_data$Feature == "Charge4", ]$Feature <- "Charge 4"
weights_data[weights_data$Feature == "DeltaScore", ]$Feature <- "Delta ALC"
weights_data[weights_data$Feature == "MainScore", ]$Feature <- "ALC"

weights_data$Feature <- factor(weights_data$Feature, levels = c("ALC", "Delta ALC", "Abs(ppm)", "Charge 1", "Charge 2", "Charge 3", "Charge 4", "SA", "Best delta RT", "RNA exp.", "Phred score"))

weight_plot <- ggplot(data=weights_data, aes(x=Sample, y=Weight, fill=Sample)) +
  theme_bw() +
  scale_fill_brewer(palette="Set3") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  staticThemeRightTop +
  theme(strip.background = element_blank(), axis.text.x = element_blank()) +
  labs(y= "Normalized weight", x="") +
  facet_grid(cols = vars(Feature))
weight_plot

ggsave("Figure2_Weights.png", plot = weight_plot, width = 24, height = 5, units = "in", dpi = 600)


######### Target decoy distribution ##############
##Figure 2D, E

feat1res <- read_excel(path = "BAAnalysis.xlsx", sheet = "Feat2_selected")
display.brewer.pal(n = 8, name = 'Paired')



generate_td_percol_plot <- function(data, sample) {
  data <- data[data$IsCanonical == F, ]
  data <- data[data$percolator_score > 0, ]
  #data <- data[data$EL_Rank < 2, ]
  data <- data[data$Sample == sample, ]
  data$Label <- as.character(data$Label)
  data[data$Label == '1', ]$Label <- "Target"
  data[data$Label == '-1', ]$Label <- "Decoy"
  data$BA <- "Non-binder"
  data[data$EL_Rank < 2, ]$BA <- "Binder"

  data[data$Label == "Target" & data$BA == "Binder", ]$Label <- "Target + Binder"
  data[data$Label == "Target" & data$BA == "Non-binder", ]$Label <- "Target + Non-binder"
  data[data$Label == "Decoy" & data$BA == "Binder", ]$Label <- "Decoy + Binder"
  data[data$Label == "Decoy" & data$BA == "Non-binder", ]$Label <- "Decoy + Non-binder"
  
  data$Label <- factor(x=data$Label, levels = c( "Target + Binder", "Target + Non-binder", "Decoy + Binder", "Decoy + Non-binder"))
  
  selectedColor <- c(brewer.pal(n = 8, name = "Paired")[2],
                     brewer.pal(n = 8, name = "Paired")[1],
                     brewer.pal(n = 8, name = "Paired")[6],
                     brewer.pal(n = 8, name = "Paired")[5])
  
  tdPlot <- ggplot(data=data, aes(x=percolator_score, fill = Label)) +
    scale_fill_manual(values = selectedColor) +
    geom_histogram(alpha = 0.6, position = "stack", bins=10) +
    scale_x_continuous(n.breaks = 3) +
    theme_bw() +
    labs(y="", x = "") +
    theme(plot.margin = margin(0.1,0.2,0.1,0.2, "in"), legend.title=element_blank()) +
    staticThemeTop
  
  return (tdPlot)
}

B_LCL1_plot <- generate_td_percol_plot(feat1res, "B-LCL1")
B_LCL2_plot <- generate_td_percol_plot(feat1res, "B-LCL2")
B_LCL3_plot <- generate_td_percol_plot(feat1res, "B-LCL3")
B_LCL4_plot <- generate_td_percol_plot(feat1res, "B-LCL4")
DOHH2_plot <- generate_td_percol_plot(feat1res, "DOHH2")
HBL1_plot <- generate_td_percol_plot(feat1res, "HBL1")
SUDHL4_plot <- generate_td_percol_plot(feat1res, "SUDHL4")
THP1_1_plot <- generate_td_percol_plot(feat1res, "THP1-1")
THP1_2_plot <- generate_td_percol_plot(feat1res, "THP1-2")
THP1_3_plot <- generate_td_percol_plot(feat1res, "THP1-3")


td_plot <- ggarrange(B_LCL1_plot, B_LCL2_plot, B_LCL3_plot, B_LCL4_plot, DOHH2_plot,
                     HBL1_plot, SUDHL4_plot, THP1_1_plot, THP1_2_plot, THP1_3_plot,
                     ncol = 5, nrow = 2, common.legend = TRUE, legend = "top")
td_plot
ggsave("Figure3_TD_Percolator_NC_0_Plot.png", plot = td_plot, width = 20, height = 6, units = "in", dpi = 600)




