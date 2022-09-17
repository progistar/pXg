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

sampleC <- c(`Subject 1` = redC, `Subject 2`=blueC, `Subject 3`=greenC, `Subject 4`=purpleC)

staticThemeLeftTop <- theme(text = element_text(size=25, color = "black")
                            , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                            legend.justification = c("right", "top"),
                            legend.position= c(.22, .98),
                            legend.text = element_text(size=20, color = "black"),
                            legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeRightTop <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                             legend.justification = c("right", "top"),
                             legend.position= c(.98, .98),
                             legend.text = element_text(size=20, color = "black"),
                             legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeRightBottom <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                             legend.justification = c("right", "bottom"),
                             legend.position= c(.98, .05),
                             legend.text = element_text(size=20, color = "black"),
                             legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeRight <- theme(text = element_text(size=25, color = "black")
                                , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                                legend.justification = c("right", "bottom"),
                                legend.position= c(.98, .25),
                                legend.text = element_text(size=20, color = "black"),
                                legend.box.background = element_rect(linetype = 1, size = 1))


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


setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/")


## RNA Read Distribution
dataS1 <- read.csv(file = "pXg/S1.RAW.PEAKS.pxg.pval.dist", header = T, sep="\t", as.is = as.double())
dataS2 <- read.csv(file = "pXg/S2.RAW.PEAKS.pxg.pval.dist", header = T, sep="\t", as.is = as.double())
dataS3 <- read.csv(file = "pXg/S3.RAW.PEAKS.pxg.pval.dist", header = T, sep="\t", as.is = as.double())
dataS4 <- read.csv(file = "pXg/S4.RAW.PEAKS.pxg.pval.dist", header = T, sep="\t", as.is = as.double())

dataS1$Subject <- "Subject 1"
dataS2$Subject <- "Subject 2"
dataS3$Subject <- "Subject 3"
dataS4$Subject <- "Subject 4"

data <- rbind(dataS1, dataS2, dataS3, dataS4)
data$Subject <- factor(data$Subject, levels = c("Subject 1", "Subject 2", "Subject 3", "Subject 4"))

data$"Peptide length" <- as.character(data$PeptideLength)
data$Mock <- log2(data$Mock+1)
data$Experiment <- log2(data$Experiment+1)
#data$`Peptide length` <- factor(data$PeptideLength, levels = c('8','9','10','11','12','13', '14', '15'))
data <- data[data$ReadCount <= 30, ]

data_length89 <- data[data$PeptideLength <= 9, ]
data_length89[data_length89$PeptideLength == 8, ]$`Peptide length` <- "Length 8"
data_length89[data_length89$PeptideLength == 9, ]$`Peptide length` <- "Length 9"

g <- ggplot(data = data_length89, aes(x=ReadCount, y=Mock, fill=`Subject`)) +
  theme_bw() +
  scale_color_manual(values=sampleC) +
  geom_line(size = 1, aes(color = Subject)) +
  xlab("Read") +
  ylab(TeX("$Log_{2}$(number of peptides + 1)")) +
  scale_x_continuous(breaks=seq(from=0, to=30, by = 1)) +
  staticThemeRightTop +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in")) +
  geom_point() +
  facet_grid(rows = vars(`Peptide length`))
g
ggsave("read.png", plot = g, width = 14, height = 9, units = "in", dpi = 300)

## Decoy ratio
data <- read_excel(path = "FDR_Analysis.xlsx", sheet = "DecoyRatio")
data$TargetRatio <- data$TargetCandidate/(data$TargetCandidate+data$DecoyCandidate)
data$Rank <- as.character(data$Rank)
data$Rank <- factor(data$Rank, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
g <- ggplot(data = data[data$IsCanonical == FALSE, ], aes(x=Rank, y=TargetRatio, fill=`Subject`)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_bar(stat = "identity", position = position_dodge())+
  xlab("Rank") +
  ylab("Ratio of noncanonical target and decoy PSMs") +
  scale_y_continuous(breaks = seq(from=0.0, to=1.0, by = 0.1)) +
  #ylab(TeX(r'($\frac{PSM_{decoy}}{PSM_{target}+PSM_{decoy}})')) +
  staticThemeTop +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in")) +
  facet_grid(rows = vars(`Cutoff`))
g

data[data$IsCanonical == FALSE & data$Rank == 1 & data$Cutoff == "p < 0.01", ]$TargetRatio

ggsave("targetratio.nc.png", plot = g, width = 14, height = 9, units = "in", dpi = 300)

## TD distribution
data <- read_excel(path = "FDR_Analysis.xlsx", sheet = "TD")
data_score <- data.frame()
data$cCount
data$Class
data$Cutoff
data$Subject


for(idx in c(1:nrow(data))) {
  d <- data[idx, ]
  while(d$cCount != 0) {
    data_score <- rbind(data_score, d)
    d$cCount <- d$cCount - 1
  }
}

median_sd1 <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

mean_sd1 <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

data_score$Class <- factor(data_score$Class, levels = c("Target", "Decoy"))
data_score$Score <- as.double(data_score$Score)
g <- ggplot(data = data_score, aes(x=Subject, y=Score, fill=`Class`)) +
  theme_bw() +
  geom_violin(trim = T) +
  scale_fill_manual(values = c(blueC, redC)) +
  scale_y_continuous(breaks = seq(from = 0, to =100, by = 10)) +
  stat_summary(fun.data = mean_sd1, geom="pointrange", color = "white", position = position_dodge2(width = 0.9), ) +
  xlab("Subject") +
  ylab("ALC score") +
  staticThemeTop +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in")) +
  facet_grid(rows = vars(`Cutoff`))
g

ggsave("TD.c.png", plot = g, width = 14, height = 9, units = "in", dpi = 300)

mean(data_score[data_score$Cutoff == "All" & data_score$Subject == "Subject 1" & data_score$Class == "Target", ]$Score)
mean(data_score[data_score$Cutoff == "p < 0.01" & data_score$Subject == "Subject 1" & data_score$Class == "Target", ]$Score)
mean(data_score[data_score$Cutoff == "All" & data_score$Subject == "Subject 1" & data_score$Class == "Decoy", ]$Score)
mean(data_score[data_score$Cutoff == "p < 0.01" & data_score$Subject == "Subject 1" & data_score$Class == "Decoy", ]$Score)

mean(data_score[data_score$Cutoff == "All" & data_score$Subject == "Subject 2" & data_score$Class == "Target", ]$Score)
mean(data_score[data_score$Cutoff == "p < 0.01" & data_score$Subject == "Subject 2" & data_score$Class == "Target", ]$Score)
mean(data_score[data_score$Cutoff == "All" & data_score$Subject == "Subject 2" & data_score$Class == "Decoy", ]$Score)
mean(data_score[data_score$Cutoff == "p < 0.01" & data_score$Subject == "Subject 2" & data_score$Class == "Decoy", ]$Score)

mean(data_score[data_score$Cutoff == "All" & data_score$Subject == "Subject 3" & data_score$Class == "Target", ]$Score)
mean(data_score[data_score$Cutoff == "p < 0.01" & data_score$Subject == "Subject 3" & data_score$Class == "Target", ]$Score)
mean(data_score[data_score$Cutoff == "All" & data_score$Subject == "Subject 3" & data_score$Class == "Decoy", ]$Score)
mean(data_score[data_score$Cutoff == "p < 0.01" & data_score$Subject == "Subject 3" & data_score$Class == "Decoy", ]$Score)

mean(data_score[data_score$Cutoff == "All" & data_score$Subject == "Subject 4" & data_score$Class == "Target", ]$Score)
mean(data_score[data_score$Cutoff == "p < 0.01" & data_score$Subject == "Subject 4" & data_score$Class == "Target", ]$Score)
mean(data_score[data_score$Cutoff == "All" & data_score$Subject == "Subject 4" & data_score$Class == "Decoy", ]$Score)
mean(data_score[data_score$Cutoff == "p < 0.01" & data_score$Subject == "Subject 4" & data_score$Class == "Decoy", ]$Score)

## Binding affinity
dataS1 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL1_FDR10")
dataS2 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL2_FDR10")
dataS3 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL3_FDR10")
dataS4 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL4_FDR10")
dataS1 <- dataS1[, -c(38:43)]
dataS2 <- dataS2[, -c(38:43)]
dataS3 <- dataS3[, -c(38:43)]
dataS4 <- dataS4[, -c(38:43)]

dataS1$Sample <- "Subject 1"
dataS2$Sample <- "Subject 2"
dataS3$Sample <- "Subject 3"
dataS4$Sample <- "Subject 4"

dataS1 <- dataS1[str_detect(dataS1$Peptide, "\\+", negate = T), ]
dataS1 <- dataS1[!duplicated(dataS1[,c('InferredPeptide')]), ]
dataS2 <- dataS2[str_detect(dataS2$Peptide, "\\+", negate = T), ]
dataS2 <- dataS2[!duplicated(dataS2[,c('InferredPeptide')]), ]
dataS3 <- dataS3[str_detect(dataS3$Peptide, "\\+", negate = T), ]
dataS3 <- dataS3[!duplicated(dataS3[,c('InferredPeptide')]), ]
dataS4 <- dataS4[str_detect(dataS4$Peptide, "\\+", negate = T), ]
dataS4 <- dataS4[!duplicated(dataS4[,c('InferredPeptide')]), ]

data <- rbind(dataS1, dataS2, dataS3, dataS4)
data$Length <- as.character(data$Length)
data$Length <- factor(data$Length, levels=c("8", "9", "10", "11", "12", "13", "14", "15"))

data$BestScore <- log2(data$BestScore+1.5)
data$Class <- "Canonical"
data[data$IsCanonical == "TRUE", ]$Class <- "Canonical"
data[data$IsCanonical == "FALSE", ]$Class <- "Noncanonical"

mapPlot <- ggplot(data=data, aes(x=Length, y=BestScore, fill=Length)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  scale_y_continuous(breaks=seq(from=0, to=7, by= 1), limits = c(0,7)) +
  staticThemeNone +
  labs(y= TeX("$Log_{2}$(1.5+%Rank)"), x = "Peptide length") +
  facet_grid(cols = vars(Sample), rows = vars(Class)) +
  geom_hline(yintercept=log2(1.5+2.0), linetype="dashed", color = "blue", size=0.5) +
  geom_hline(yintercept=log2(1.5+0.5), linetype="dashed", color = "red", size=0.5)

mapPlot
ggsave("BA.png", plot = mapPlot, width = 14, height = 8, units = "in", dpi = 300)


## Ided PSMs according to rank
dataS1 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL1_FDR10")
dataS2 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL2_FDR10")
dataS3 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL3_FDR10")
dataS4 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL4_FDR10")
dataS1 <- dataS1[, -c(38:43)]
dataS2 <- dataS2[, -c(38:43)]
dataS3 <- dataS3[, -c(38:43)]
dataS4 <- dataS4[, -c(38:43)]

dataS1$Sample <- "Subject 1"
dataS2$Sample <- "Subject 2"
dataS3$Sample <- "Subject 3"
dataS4$Sample <- "Subject 4"

data <- rbind(dataS1, dataS2, dataS3, dataS4)

data$Class <- "Canonical"
data[data$IsCanonical == "TRUE", ]$Class <- "Canonical"
data[data$IsCanonical == "FALSE", ]$Class <- "Noncanonical"
data$Type <- "NB"
data[data$BestScore < 2, ]$Type <- "WB"
data[data$BestScore < 0.5, ]$Type <- "SB"
data_rank <- data.frame(Rank=numeric(), Type=character(), Class=character(), Subject=character(), PSM=numeric())

## Canonical ratio
nrow(data[data$Sample == "Subject 1" & data$Type != "NB" & data$Class == "Canonical", ])/nrow(data[data$Sample == "Subject 1" & data$Class == "Canonical", ])
nrow(data[data$Sample == "Subject 2" & data$Type != "NB" & data$Class == "Canonical", ])/nrow(data[data$Sample == "Subject 2" & data$Class == "Canonical", ])
nrow(data[data$Sample == "Subject 3" & data$Type != "NB" & data$Class == "Canonical", ])/nrow(data[data$Sample == "Subject 3" & data$Class == "Canonical", ])
nrow(data[data$Sample == "Subject 4" & data$Type != "NB" & data$Class == "Canonical", ])/nrow(data[data$Sample == "Subject 4" & data$Class == "Canonical", ])

## Noncanonical ratio
nrow(data[data$Sample == "Subject 1" & data$Type != "NB" & data$Class == "Noncanonical", ])/nrow(data[data$Sample == "Subject 1" & data$Class == "Noncanonical", ])
nrow(data[data$Sample == "Subject 2" & data$Type != "NB" & data$Class == "Noncanonical", ])/nrow(data[data$Sample == "Subject 2" & data$Class == "Noncanonical", ])
nrow(data[data$Sample == "Subject 3" & data$Type != "NB" & data$Class == "Noncanonical", ])/nrow(data[data$Sample == "Subject 3" & data$Class == "Noncanonical", ])
nrow(data[data$Sample == "Subject 4" & data$Type != "NB" & data$Class == "Noncanonical", ])/nrow(data[data$Sample == "Subject 4" & data$Class == "Noncanonical", ])


for(idx in c(1: nrow(data))) {
  sub_data <- data[idx, ]
  rank <- sub_data$Rank
  type <- sub_data$Type
  class_ <- sub_data$Class
  subject <- sub_data$Sample
  
  if(nrow(data_rank[data_rank$Rank == rank & data_rank$Type == type & data_rank$Class == class_ & data_rank$Subject == subject, ]) == 0){
    data_rank[nrow(data_rank)+1, ] = c(rank, type, class_, subject, 1)
    
  } else {
    data_rank[data_rank$Rank == rank & data_rank$Type == type & data_rank$Class == class_ & data_rank$Subject == subject, ]$PSM <- as.numeric(data_rank[data_rank$Rank == rank & data_rank$Type == type & data_rank$Class == class_ & data_rank$Subject == subject, ]$PSM)+1
  }
}
data_rank$Rank <- factor(data_rank$Rank, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
data_rank$PSM <- as.numeric(data_rank$PSM)
data_rank$Type <- factor(data_rank$Type, levels = c("SB", "WB", "NB"))

g <- ggplot(data = data_rank, aes(x=Subject, y=PSM, fill=Type)) +
  theme_bw() +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c(blueC, greenC, redC)) +
  ylab("PSM") +
  staticThemeRight +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in")) +
  facet_grid(rows = vars(`Class`), scales = "free")
g

ggsave("PSMs.png", plot = g, width = 8, height = 8, units = "in", dpi = 300)

g <- ggplot(data = data_rank, aes(x=Rank, y=PSM, fill=Type)) +
  theme_bw() +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c(blueC, greenC, redC)) +
  #scale_y_continuous(breaks = seq(from=0, to=1, by=0.2)) +
  xlab("Rank") +
  ylab("PSM") +
  staticThemeRightBottom +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in")) +
  facet_grid(cols = vars(`Subject`), rows = vars(Class), scales = "free")
g

ggsave("Rank.png", plot = g, width = 14, height = 8, units = "in", dpi = 300)

g <- ggplot(data = data_rank, aes(x=Class, y=PSM, fill=Rank)) +
  theme_bw() +
  geom_bar(position="fill", stat = "identity")+
  scale_fill_brewer(palette = "Set3")+
  scale_y_continuous(breaks = seq(from=0, to=1, by=0.2)) +
  ylab("Proportion") +
  staticThemeTop +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in")) +
  facet_grid(cols = vars(`Subject`))
g

ggsave("Rank2.png", plot = g, width = 16, height = 8, units = "in", dpi = 300)




## Sequence Length Distribution
dataS1 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL1_FDR10")
dataS2 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL2_FDR10")
dataS3 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL3_FDR10")
dataS4 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL4_FDR10")
dataS1 <- dataS1[, -c(38:43)]
dataS2 <- dataS2[, -c(38:43)]
dataS3 <- dataS3[, -c(38:43)]
dataS4 <- dataS4[, -c(38:43)]

dataS1$Subject <- "Subject 1"
dataS2$Subject <- "Subject 2"
dataS3$Subject <- "Subject 3"
dataS4$Subject <- "Subject 4"

dataS1 <- dataS1[str_detect(dataS1$Peptide, "\\+", negate = T), ]
dataS1 <- dataS1[!duplicated(dataS1[,c('InferredPeptide')]), ]
dataS2 <- dataS2[str_detect(dataS2$Peptide, "\\+", negate = T), ]
dataS2 <- dataS2[!duplicated(dataS2[,c('InferredPeptide')]), ]
dataS3 <- dataS3[str_detect(dataS3$Peptide, "\\+", negate = T), ]
dataS3 <- dataS3[!duplicated(dataS3[,c('InferredPeptide')]), ]
dataS4 <- dataS4[str_detect(dataS4$Peptide, "\\+", negate = T), ]
dataS4 <- dataS4[!duplicated(dataS4[,c('InferredPeptide')]), ]

data <- rbind(dataS1, dataS2, dataS3, dataS4)

data$Class <- "Canonical"
data[data$IsCanonical == "TRUE", ]$Class <- "Canonical"
data[data$IsCanonical == "FALSE", ]$Class <- "Noncanonical"
data$Type <- "NB"
data[data$BestScore < 2, ]$Type <- "WB"
data[data$BestScore < 0.5, ]$Type <- "SB"

dataS12 <- data[data$Subject == "Subject 1" | data$Subject == "Subject 2", ]
dataS34 <- data[data$Subject == "Subject 3" | data$Subject == "Subject 4", ]
dataS1 <- data[data$Subject == "Subject 1", ]
dataS2 <- data[data$Subject == "Subject 2", ]
dataS3 <- data[data$Subject == "Subject 3", ]
dataS4 <- data[data$Subject == "Subject 4", ]

g <- ggplot(data = dataS4[dataS4$BestScore < 2 & dataS4$Length < 14, ], aes(x=Length, fill=BestType)) +
  theme_bw() +
  geom_histogram(position = position_dodge2(preserve = "single", padding = 0))+
  scale_fill_brewer(palette = "Set1")+
  #scale_y_continuous(breaks = seq(from=0, to=1, by=0.2)) +
  scale_x_continuous(breaks = seq(from = 8, to = 13, by = 1)) +
  ylab("MAP") +
  staticThemeTop +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in")) +
  facet_grid(rows=vars(`Class`) ,scales = "free")
g
ggsave("LengthDistribution_S4.png", plot = g, width = 14, height = 8, units = "in", dpi = 300)



  
