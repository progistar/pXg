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

## Load data
setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/")
dataS1 <- read_excel(path = "pXg_Summary.xlsx", sheet = "B-LCL1_FDR5")
dataS2 <- read_excel(path = "pXg_Summary.xlsx", sheet = "B-LCL2_FDR5")
dataS3 <- read_excel(path = "pXg_Summary.xlsx", sheet = "B-LCL3_FDR5")
dataS4 <- read_excel(path = "pXg_Summary.xlsx", sheet = "B-LCL4_FDR5")

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

subData <- allData[allData$BestScore < 2, ]
subData <- subData[subData$IsCanonical == "FALSE", ]

## Draw
topRankedPlot <- ggplot(data=subData, aes(x=BestType, fill=Sample)) +
  scale_fill_brewer(palette="Set1", name ="Sample") +
  theme_bw() +
  geom_bar(aes(y= (..count..), x=BestType), position=position_dodge2(preserve = "single", padding = 0), width = 0.5) +
  theme(text = element_text(size=25, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"), 
        axis.text.x = element_text(size = 15, colour = "black", angle = 90))+
  labs(y="ncMAPs", x = "\nHLA type") +
  ggtitle("ncMAPs according to HLA type")

topRankedPlot
ggsave("BA.HLA.ncMAPs.png", plot = topRankedPlot, width = 10, height = 8, units = "in", dpi = 600)

nrow(subData[subData$BestType == "HLA-C05:01" & subData$Sample == "B-LCL subject 4", ])


