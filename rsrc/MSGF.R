## Histogram for score distribution
library(ggplot2)
library(RColorBrewer)
library("dplyr")
library(readxl)
library(latex2exp)
library("stringr")       

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

setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/MSGF/unmodified")
## BA distribution
dataS1 <- read_excel(path = "MSGF_Summary.xlsx", sheet = "B-LCL1_FDR5")
dataS2 <- read_excel(path = "MSGF_Summary.xlsx", sheet = "B-LCL2_FDR5")
dataS3 <- read_excel(path = "MSGF_Summary.xlsx", sheet = "B-LCL3_FDR5")
dataS4 <- read_excel(path = "MSGF_Summary.xlsx", sheet = "B-LCL4_FDR5")

dataS1 <- dataS1[, -c(20:25)]
dataS2 <- dataS2[, -c(20:25)]
dataS3 <- dataS3[, -c(20:25)]
dataS4 <- dataS4[, -c(20:25)]

dataS1$Sample <- "B-LCL subject 1"
dataS2$Sample <- "B-LCL subject 2"
dataS3$Sample <- "B-LCL subject 3"
dataS4$Sample <- "B-LCL subject 4"

dataS1 <- dataS1[str_detect(dataS1$Peptide, "\\+", negate = T), ]
dataS1 <- dataS1[!duplicated(dataS1[,c('Peptide')]), ]
dataS2 <- dataS2[str_detect(dataS2$Peptide, "\\+", negate = T), ]
dataS2 <- dataS2[!duplicated(dataS2[,c('Peptide')]), ]
dataS3 <- dataS3[str_detect(dataS3$Peptide, "\\+", negate = T), ]
dataS3 <- dataS3[!duplicated(dataS3[,c('Peptide')]), ]
dataS4 <- dataS4[str_detect(dataS4$Peptide, "\\+", negate = T), ]
dataS4 <- dataS4[!duplicated(dataS4[,c('Peptide')]), ]

allData <- rbind(dataS1, dataS2, dataS3, dataS4)

subData <- allData

## select unmodi
subData$BestScore <- log2(subData$BestScore+1.5)

subData$Length <- as.character(subData$Length)
subData$Length <- factor(subData$Length, levels=c("8", "9", "10", "11", "12", "13", "14", "15"))

nrow(subData)

mapPlot <- ggplot(data=subData, aes(x=Length, y=BestScore, fill=Length)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  ggtitle("HLA binding affinity for canonical peptides (MS-GF+)") +
  scale_y_continuous(breaks=seq(from=0, to=7, by= 1)) +
  labs(y= TeX("$Log_{2}$(1.5+%Rank)"), x = "Peptide length") +
  annotate("text", label = nrow(subData[subData$Length == 8,]), x =1, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 9,]), x =2, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 10,]), x =3, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 11,]), x =4, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 12,]), x =5, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 13,]), x =6, y = 7, size = 8, family="serif", colour = "black") +
  geom_hline(yintercept=log2(1.5+2.0), linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=log2(1.5+0.5), linetype="dashed", color = "black", size=1) +
  annotate("text", label = nrow(subData[subData$Length == 14,]), x =7, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 15,]), x =8, y = 7, size = 8, family="serif", colour = "black")

mapPlot
ggsave("MSGF.BA.Canonical.png", plot = mapPlot, width = 10, height = 8, units = "in", dpi = 300)

nrow(subData[subData$BestScore < 1.5+2.0, ])
nrow(subData)

## MS-GF+ vs pXg length distribution
dataS1 <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg_Summary.xlsx", sheet = "B-LCL1_FDR5")
dataS2 <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg_Summary.xlsx", sheet = "B-LCL2_FDR5")
dataS3 <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg_Summary.xlsx", sheet = "B-LCL3_FDR5")
dataS4 <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg_Summary.xlsx", sheet = "B-LCL4_FDR5")

dataS1 <- dataS1[, -c(38:43)]
dataS2 <- dataS2[, -c(38:43)]
dataS3 <- dataS3[, -c(38:43)]
dataS4 <- dataS4[, -c(38:43)]

dataS1$Sample <- "B-LCL subject 1"
dataS2$Sample <- "B-LCL subject 2"
dataS3$Sample <- "B-LCL subject 3"
dataS4$Sample <- "B-LCL subject 4"

dataS1 <- dataS1[str_detect(dataS1$Peptide, "\\+", negate = T), ]
dataS1 <- dataS1[!duplicated(dataS1[,c('InferredPeptide')]), ]
dataS2 <- dataS2[str_detect(dataS2$Peptide, "\\+", negate = T), ]
dataS2 <- dataS2[!duplicated(dataS2[,c('InferredPeptide')]), ]
dataS3 <- dataS3[str_detect(dataS3$Peptide, "\\+", negate = T), ]
dataS3 <- dataS3[!duplicated(dataS3[,c('InferredPeptide')]), ]
dataS4 <- dataS4[str_detect(dataS4$Peptide, "\\+", negate = T), ]
dataS4 <- dataS4[!duplicated(dataS4[,c('InferredPeptide')]), ]

allData <- rbind(dataS1, dataS2, dataS3, dataS4)

subDataPXG <- allData
subDataPXG$BestScore <- log2(subDataPXG$BestScore+1.5)

subData <- subData[subData$BestScore < log2(1.5+2.0), ]
subDataPXG <- subDataPXG[subDataPXG$BestScore < log2(1.5+2.0), ]

subData$Class <- "cMAP (MS-GF+)"
subDataPXG$Class <- "cMAP (pXg)"
subDataPXG[subDataPXG$IsCanonical == "FALSE", ]$Class <- "ncMAP (pXg)"


subData <- rbind(subData[, c('Length', 'Class')], subDataPXG[,c('Length', 'Class')])
subData$Length <- as.numeric(subData$Length)

topRankedPlot <- ggplot(data=subData, aes(x=Length, fill=Class)) +
  scale_fill_manual(values=c(greenC, blueC, redC)) +
  theme_bw() +
  #geom_bar(aes(y= (..prop..), x=Length), position=position_dodge(), width = 0.5) +
  geom_bar(aes(y= (..count..), x=Length), position=position_dodge(), width = 0.5) +
  scale_x_continuous(breaks = seq(from=8, to=15, by=1)) +
  #scale_y_continuous(breaks=seq(from=0, to=1, by = 0.1), 
  #                   labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE),
  #                    limits=c(0,0.8)) +
  theme(text = element_text(size=20)) +
  labs(y="MAPs", x = "Length of MAP") +
  #annotate("text", label = "**", x =8, y = 0.22, size = 8, family="serif", colour = "black") +
  #annotate("text", label = "**", x =9, y = 0.58, size = 8, family="serif", colour = "black") +
  #annotate("text", label = "**", x =10, y = 0.25, size = 8, family="serif", colour = "black") +
  #annotate("text", label = "**", x =11, y = 0.13, size = 8, family="serif", colour = "black") +
  #annotate("text", label = "**", x =12, y = 0.05, size = 8, family="serif", colour = "black") +
  #annotate("text", label = "**", x =13, y = 0.03, size = 8, family="serif", colour = "black") +
  #geom_segment(aes(x = 7.75, xend =8.25, y = 0.21, yend = 0.21), color = "Black", size = 1) +
  #geom_segment(aes(x = 8.75, xend =9.25, y = 0.57, yend = 0.57), color = "Black", size = 1) +
  #geom_segment(aes(x = 9.75, xend =10.25, y = 0.24, yend = 0.24), color = "Black", size = 1) +
  #geom_segment(aes(x = 10.75, xend =11.25, y = 0.12, yend = 0.12), color = "Black", size = 1) +
  #geom_segment(aes(x = 11.75, xend =12.25, y = 0.04, yend = 0.04), color = "Black", size = 1) +
  #geom_segment(aes(x = 12.75, xend =13.25, y = 0.02, yend = 0.02), color = "Black", size = 1) +
  ggtitle("Length distributions of MAPs")

topRankedPlot
ggsave("MSGF_pXg.Length.png", plot = topRankedPlot, width = 10, height = 8, units = "in", dpi = 300)

len = 15
nrow(subData[subData$Class == "cMAP (MS-GF+)" & subData$Length == len, ])
nrow(subData[subData$Class == "cMAP (pXg)" & subData$Length == len, ])
nrow(subData[subData$Class == "ncMAP (pXg)" & subData$Length == len, ])

len = 8
matrix <- matrix(c(nrow(subData[subData$Class == "MS-GF+", ]) - nrow(subData[subData$Length == len & subData$Class == "MS-GF+", ]), 
                   nrow(subData[subData$Length == len & subData$Class == "MS-GF+", ]),
                   nrow(subData[subData$Class == "pXg", ]) - nrow(subData[subData$Length == len & subData$Class == "pXg", ]), 
                   nrow(subData[subData$Length == len & subData$Class == "pXg", ])), 
                 nrow = 2, ncol = 2)
matrix <- t(matrix)
fisher.test(matrix, alternative = "two.sided")

nrow(subData[subData$Length == 15 & subData$Class == "MS-GF+", ])

2^3.5-1.5
log2(11.5)
