install.packages("ggpmisc")
install.packages("heatmaply")

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


setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/")
## DeepLC
rtUnmatched <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXgResults_p001_FDR.xlsx", sheet = "Unmatched_RT_85")
rtMatched <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXgResults_p001_FDR.xlsx", sheet = "pFDR_RT")

plotData <- rtMatched[rtMatched$IsCanonical == TRUE & rtMatched$`MHC-I` != "NB", c("RT", "deeplcRT")]
plotData$Class <- "Canonical MAP"

plotData1 <- rtMatched[rtMatched$IsCanonical == FALSE & rtMatched$`MHC-I` != "NB", c("RT", "deeplcRT")]
plotData1$Class <- "Noncanonical MAP"

plotData2 <- rtUnmatched[rtUnmatched$`ALC (%)` >= 85 & rtUnmatched$`MHC-I` != "NB", c("RT", "deeplcRT")]
plotData2$Class <- "Unmatched MAP (ALC≥85)"

plotData <- rbind(plotData, plotData1, plotData2)
plotData$Class <- factor(plotData$Class, levels = c("Unmatched MAP (ALC≥85)", "Canonical MAP", "Noncanonical MAP"))

nrow(plotData[plotData$Class == "Unmatched MAP (ALC≥85)", ])
nrow(plotData[plotData$Class == "Noncanonical MAP", ])
nrow(plotData[plotData$Class == "Canonical MAP", ])
?lm
summary(lm(deeplcRT~RT, plotData[plotData$Class == "Unmatched MAP (ALC≥85)",]))
summary(lm(deeplcRT~RT, plotData[plotData$Class == "Noncanonical MAP",]))

rtPlot <- ggplot(data = plotData, aes(x=RT, y=`deeplcRT`)) +
  theme_bw() +
  geom_point(size=1, alpha = 0.6) +
  stat_smooth(method = "lm", se = T, size = 0.8) +
  staticThemeRight +
  theme(plot.margin = margin(0, 0.5, 0, 0, "in")) +
  labs(y= "Predicted retention time", x = NULL) +
  labs(x= "Observed retention time", x = NULL) +
  facet_grid(cols = vars(`Class`), scales = "free")
rtPlot
ggsave("RT_Unmatched_85.png", plot = rtPlot, width = 15, height = 7, units = "in", dpi = 1200)



## Peptide length
lenData1 <- data.frame(matrix(ncol = 3, nrow = 0))
lenData2 <- data.frame(matrix(ncol = 3, nrow = 0))
lenData3 <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(lenData1) <- c("Length", "Count","Class")
colnames(lenData2) <- c("Length", "Count","Class")
colnames(lenData3) <- c("Length", "Count","Class")

plotData1 <- rtMatched[rtMatched$IsCanonical == TRUE & rtMatched$`MHC-I` != "NA", c("Length")]
plotData2 <- rtMatched[rtMatched$IsCanonical == FALSE & rtMatched$`MHC-I` != "NA", c("Length")]
plotData3 <- rtUnmatched[rtUnmatched$`ALC (%)` >= 85 & rtUnmatched$`MHC-I` != "NA", c("Length")]

plotSum1 <- 0
plotSum2 <- 0
plotSum3 <- 0
for(idx in c(8,9,10,11,12,13,14,15)) {
  lenData1[idx-7, ]$Length <- idx
  lenData1[idx-7, ]$Count <- nrow(plotData1[plotData1$Length == idx, ])
  lenData1[idx-7, ]$Class <- "Canonical"
  plotSum1 <- lenData1[idx-7, ]$Count + plotSum1
  
  lenData2[idx-7, ]$Length <- idx
  lenData2[idx-7, ]$Count <- nrow(plotData2[plotData2$Length == idx, ])
  lenData2[idx-7, ]$Class <- "Noncanonical"
  plotSum2 <- lenData2[idx-7, ]$Count + plotSum2
  
  lenData3[idx-7, ]$Length <- idx
  lenData3[idx-7, ]$Count <- nrow(plotData3[plotData3$Length == idx, ])
  lenData3[idx-7, ]$Class <- "Unmatched (ALC≥85)"
  plotSum3 <- lenData3[idx-7, ]$Count + plotSum3
}
for(idx in c(8,9,10,11,12,13,14,15)) {
  lenData1[idx-7, ]$Count <- lenData1[idx-7, ]$Count / plotSum1
  lenData2[idx-7, ]$Count <- lenData2[idx-7, ]$Count / plotSum2
  lenData3[idx-7, ]$Count <- lenData3[idx-7, ]$Count / plotSum3
}

lenData <- rbind(lenData1, lenData2)
lenData$Class <- factor(lenData$Class, levels = c("Canonical", "Noncanonical"))

#lenData$Length <- as.character(lenData$Length)
classCa <- c(`Canonical` = blueC, `Noncanonical` = redC)

lengthPlot <- ggplot(data = lenData, aes(x=Length, y=Count, fill = Class)) +
  geom_bar(stat="identity", position = position_dodge2()) +
  scale_x_continuous(breaks = c(8,9,10,11,12,13,14,15)) +
  theme_bw() +
  staticThemeTop +
  labs(y= "MAP proportion", x = NULL) +
  labs(x= "Length", x = NULL) +
  scale_fill_manual(values = classCa)
lengthPlot

ggsave("length.png", plot = lengthPlot, width = 6, height = 7, units = "in", dpi = 300)


## PPM error
ppmUnmatched <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXgResults_p001_FDR.xlsx", sheet = "Unmatched_ppm")
ppmMatched <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXgResults_p001_FDR.xlsx", sheet = "pFDR_ppm")

ppmMatched$Class <- "Canonical MAP"
ppmMatched[ppmMatched$IsCanonical == FALSE, ]$Class <- "Noncanonical MAP"
ppmUnmatched$Class <- "Unmatched MAP (ALC≥85)"

plotData <- ppmMatched[ppmMatched$`MHC-I` != "NB", c("ppm", "Class")]
plotData <- rbind(plotData, ppmUnmatched[ppmUnmatched$`ALC (%)` >=85 & 
                                           ppmUnmatched$`MHC-I` != "NB", 
                                         c("ppm", "Class")])

summary(plotData)

plotData$Class <- factor(plotData$Class, labels = c("Canonical MAP", "Noncanonical MAP", "Unmatched MAP (ALC≥85)"))

nrow(plotData[plotData$Class == "Unmatched MAP (ALC≥85)",])
nrow(plotData[plotData$Class == "Noncanonical MAP",])
nrow(plotData[plotData$Class == "Canonical MAP",])

classCa <- c(`Canonical MAP` = blueC, `Noncanonical MAP` = redC, `Unmatched MAP (ALC≥85)` = grayC)

ppmPlot <- ggplot(data = plotData, aes(x=ppm, color = Class)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(values = classCa) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  staticThemeTop +
  labs(y= "PSM density", x = NULL) +
  labs(x= "ppm error", x = NULL)
  
ppmPlot
ggsave("ppm_Unmatched_85.png", plot = ppmPlot, width = 9, height = 8, units = "in", dpi = 1200)

##
tmp <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXgResults_nocut.xlsx", sheet = "Subject 4")
tmp <- tmp[!duplicated(tmp$InferredPeptide), ]
nrow(tmp[tmp$IsCanonical == "TRUE" & tmp$`MHC-I` != "NB", ])
nrow(tmp[tmp$IsCanonical == "TRUE", ])
nrow(tmp[tmp$IsCanonical == "FALSE" & tmp$`MHC-I` != "NB" & tmp$`Denovo score` >= 94, ])
nrow(tmp[tmp$IsCanonical == "FALSE" & tmp$`Denovo score` >= 94, ])

binderData <- read_excel(path = "OPF_Analysis.xlsx", sheet = "Binder")
binderData$Type <- factor(binderData$Type, levels = c("Unmatched", "At least 1", "p < 0.01", "At least 1 + FDR", "p < 0.01 + FDR"))

binderPlot <- ggplot(data = binderData, aes(x=Type, y=`Binder ratio`, group = Subject)) +
  theme_bw() +
  geom_line(aes(color = Subject)) +
  geom_point(aes(color = Subject)) +
  staticThemeRight +
  scale_y_continuous(breaks = seq(from=0, to=1, by=0.2), limits = c(0,1)) +
  labs(y= "Binder ratio", x = NULL) +
  facet_grid(rows = vars(`Class`)) +
  scale_color_manual(values = sampleC) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
binderPlot

ggsave("BinderRatio.png", plot = binderPlot, width = 7, height = 6, units = "in", dpi = 300)

IDData <- read_excel(path = "OPF_Analysis.xlsx", sheet = "ID")

idPlot <- ggplot(data = IDData, aes(x=Type, y=MAP, group = Subject)) +
  theme_bw() +
  geom_line(aes(color = Subject)) +
  geom_point(aes(color = Subject)) +
  staticThemeNone +
  labs(y= "Identified MAP", x = NULL) +
  scale_color_manual(values = sampleC) +
  facet_grid(rows = vars(`Class`), scales = "free") +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))

idPlot
ggsave("IDCount.png", plot = idPlot, width = 4, height = 6, units = "in", dpi = 300)

dataS1 <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "Subject 1")
dataS2 <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "Subject 2")
dataS3 <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "Subject 3")
dataS4 <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "Subject 4")
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
data_rank$isTop <- "Rank N"
data_rank[data_rank$Rank == "1", ]$isTop <- "Rank 1"

inc <- sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Canonical"  & data_rank$Subject == "Subject 1", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 1" & data_rank$Class == "Canonical", ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Canonical"  & data_rank$Subject == "Subject 2", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 2" & data_rank$Class == "Canonical", ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Canonical"  & data_rank$Subject == "Subject 3", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 3" & data_rank$Class == "Canonical", ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Canonical"  & data_rank$Subject == "Subject 4", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 4" & data_rank$Class == "Canonical", ]$PSM)
inc/4

inc <- sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Noncanonical"  & data_rank$Subject == "Subject 1", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 1" & data_rank$Class == "Noncanonical" , ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Noncanonical"  & data_rank$Subject == "Subject 2", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 2" & data_rank$Class == "Noncanonical" , ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Noncanonical"  & data_rank$Subject == "Subject 3", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 3" & data_rank$Class == "Noncanonical" , ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Noncanonical"  & data_rank$Subject == "Subject 4", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 4" & data_rank$Class == "Noncanonical" , ]$PSM)
inc/4

g <- ggplot(data = data_rank[data_rank$Type != "NB", ], aes(x=Subject, y=PSM, fill=isTop)) +
  theme_bw() +
  geom_bar(stat = "identity", width = 0.8)+
  scale_fill_manual(values = c(blueC, greenC)) +
  ylab("PSM") +
  xlab(NULL) +
  staticThemeRight +
  facet_grid(rows = vars(`Class`), scales = "free") +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +
  labs(fill = "Rank")
g

ggsave("PSMs.png", plot = g, width = 7, height = 5.5, units = "in", dpi = 300)


#figure2 <- ggarrange(binderPlot, idPlot, ncol = 2, nrow = 1, common.legend = T)
#ggsave("Figure2.png", plot = figure2, width = 12, height = 6, units = "in", dpi = 300)

## Event catalog
dataS1 <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "Subject 1")
dataS2 <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "Subject 2")
dataS3 <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "Subject 3")
dataS4 <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "Subject 4")

dataS1 <- dataS1[, -c(38:43)]
dataS2 <- dataS2[, -c(38:43)]
dataS3 <- dataS3[, -c(38:43)]
dataS4 <- dataS4[, -c(38:43)]

dataS1$Subject <- "Subject 1"
dataS2$Subject <- "Subject 2"
dataS3$Subject <- "Subject 3"
dataS4$Subject <- "Subject 4"

allData <- rbind(dataS1, dataS2, dataS3, dataS4)

subData <- allData[allData$BestScore < 2, ]
subData$Type <- "WB"
subData[subData$BestScore < 2, ]$Type <- "WB"
subData[subData$BestScore < 0.5, ]$Type <- "SB"
subData$Class <- "Canonical"
subData[subData$IsCanonical == "FALSE", ]$Class <- "Noncanonical"
subData$Subject <- factor(subData$Subject, levels = c("Subject 1", "Subject 2", "Subject 3", "Subject 4"))

### Unique selection
subDataUnique <- subData[subData$Class == "Noncanonical", ]

subDataUnique <- subDataUnique[subDataUnique$EventCount == 1, ]
subDataUnique <- subDataUnique[subDataUnique$GeneIDCount <= 1, ]
subDataUnique <- subDataUnique[subDataUnique$GenomicLociCount <= 1, ]

tmp <- subDataUnique[subDataUnique$GeneIDCount > 1 | subDataUnique$EventCount > 1 | subDataUnique$GenomicLociCount > 1, ]
nrow(tmp[!duplicated(tmp$InferredPeptide), ])

tmp <- tmp[!duplicated(tmp$InferredPeptide), ]

# meta characteristics
subDataUnique$AS <- FALSE
subDataUnique$asRNA <- FALSE

subDataUnique[subDataUnique$Events == "PC;sense;AS", ]$AS <- TRUE
subDataUnique[subDataUnique$Events == "ncRNA;asRNA;AS", ]$AS <- TRUE
subDataUnique[subDataUnique$Events == "ncRNA;sense;AS", ]$AS <- TRUE
subDataUnique[subDataUnique$Events == "IR;sense;AS", ]$AS <- TRUE
subDataUnique[subDataUnique$Events == "IR;asRNA;AS", ]$AS <- TRUE

subDataUnique[subDataUnique$Events == "PC;asRNA", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "5`-UTR;asRNA", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "3`-UTR;asRNA", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "IR;asRNA", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "IR;asRNA;AS", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "ncRNA;asRNA;AS", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "ncRNA;asRNA", ]$asRNA <- TRUE


subDataUnique[subDataUnique$Events == "5`-UTR;sense", ]$Events <- "5`-UTR"
subDataUnique[subDataUnique$Events == "3`-UTR;sense", ]$Events <- "3`-UTR"
subDataUnique[subDataUnique$Events == "ncRNA;sense", ]$Events <- "ncRNA"
subDataUnique[subDataUnique$Events == "FS;sense", ]$Events <- "FS"
subDataUnique[subDataUnique$Events == "IR;sense", ]$Events <- "IR"
subDataUnique[subDataUnique$Events == "IR;sense;AS", ]$Events <- "IR"
subDataUnique[subDataUnique$Events == "IR;asRNA;AS", ]$Events <- "IR"
subDataUnique[subDataUnique$Events == "IGR;sense", ]$Events <- "IGR"
subDataUnique[subDataUnique$Events == "PC;asRNA", ]$Events <- "Coding"
subDataUnique[subDataUnique$Events == "PC;sense;AS", ]$Events <- "Coding"
subDataUnique[subDataUnique$Events == "ncRNA;sense;AS", ]$Events <- "ncRNA"
subDataUnique[subDataUnique$Events == "ncRNA;asRNA;AS", ]$Events <- "ncRNA"
subDataUnique[subDataUnique$Events == "5`-UTR;asRNA", ]$Events <- "5`-UTR"
subDataUnique[subDataUnique$Events == "3`-UTR;asRNA", ]$Events <- "3`-UTR"
subDataUnique[subDataUnique$Events == "IR;asRNA", ]$Events <- "IR"
subDataUnique[subDataUnique$Events == "ncRNA;asRNA", ]$Events <- "ncRNA"
subDataUnique[subDataUnique$Events == "unknown", ]$Events <- "Unknown"
subDataUnique[subDataUnique$Events == "PC;sense", ]$Events <- "Coding"

subDataUnique$Mutation <- TRUE
subDataUnique[subDataUnique$Mutations == "-", ]$Mutation <- FALSE

subDataUnique$Events <- factor(subDataUnique$Events, 
                                             levels = c("Coding", "5`-UTR", "ncRNA", "FS", "3`-UTR", "IR", "IGR", 
                                                        "Unknown"))

subDataUnique$InferredPeptide
subDataUnique$Reads

peptideLevel <- ddply(subDataUnique, .(InferredPeptide, Subject, asRNA, AS, Mutation, GenomicLoci, Events, Reads), nrow)
peptideLevel$`Subject 1` <- 0
peptideLevel$`Subject 2` <- 0
peptideLevel$`Subject 3` <- 0
peptideLevel$`Subject 4` <- 0

peptideLevel[peptideLevel$Subject == "Subject 1",]$`Subject 1` <- peptideLevel[peptideLevel$Subject == "Subject 1",]$V1
peptideLevel[peptideLevel$Subject == "Subject 2",]$`Subject 2` <- peptideLevel[peptideLevel$Subject == "Subject 2",]$V1
peptideLevel[peptideLevel$Subject == "Subject 3",]$`Subject 3` <- peptideLevel[peptideLevel$Subject == "Subject 3",]$V1
peptideLevel[peptideLevel$Subject == "Subject 4",]$`Subject 4` <- peptideLevel[peptideLevel$Subject == "Subject 4",]$V1

#write.table(peptideLevel, file = "ncMAPs.tsv", sep = "\t")

peptideLevel <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "ncMAP")
peptideLevel$Events <- factor(peptideLevel$Events, 
                              levels = c("Coding", "5`-UTR", "FS", "IR", "ncRNA", "3`-UTR", "IGR", "Unknown"))
peptideLevel <- peptideLevel[order(peptideLevel$Events, -peptideLevel$PSM1, -peptideLevel$PSM2, -peptideLevel$PSM3, -peptideLevel$PSM4),]

peptideLevel$Event <- 8
peptideLevel[peptideLevel$Events == "5`-UTR", ]$Event <- 7
peptideLevel[peptideLevel$Events == "FS", ]$Event <- 6
peptideLevel[peptideLevel$Events == "IR", ]$Event <- 5
peptideLevel[peptideLevel$Events == "ncRNA", ]$Event <- 4
peptideLevel[peptideLevel$Events == "3`-UTR", ]$Event <- 3
peptideLevel[peptideLevel$Events == "IGR", ]$Event <- 2
peptideLevel[peptideLevel$Events == "Unknown", ]$Event <- 1

peptideLevel[peptideLevel$PSM1 > 2, ]$PSM1 <- 3
peptideLevel[peptideLevel$PSM2 > 2, ]$PSM2 <- 3
peptideLevel[peptideLevel$PSM3 > 2, ]$PSM3 <- 3
peptideLevel[peptideLevel$PSM4 > 2, ]$PSM4 <- 3

peptideLevel$Read1 <- log2(peptideLevel$Read1+1)
peptideLevel$Read2 <- log2(peptideLevel$Read2+1)
peptideLevel$Read3 <- log2(peptideLevel$Read3+1)
peptideLevel$Read4 <- log2(peptideLevel$Read4+1)

nrow(peptideLevel[peptideLevel$Event == 8, ])
nrow(peptideLevel[peptideLevel$Event == 7, ])
nrow(peptideLevel[peptideLevel$Event == 6, ])
nrow(peptideLevel[peptideLevel$Event == 5, ])
nrow(peptideLevel[peptideLevel$Event == 4, ])
nrow(peptideLevel[peptideLevel$Event == 3, ])
nrow(peptideLevel[peptideLevel$Event == 2, ])
nrow(peptideLevel[peptideLevel$Event == 1, ])
nrow(peptideLevel[peptideLevel$asRNA1 == 1 | peptideLevel$asRNA2 == 1 | peptideLevel$asRNA3 == 1 | peptideLevel$asRNA4 == 1, ])

nrow(peptideLevel[peptideLevel$Event == 8 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 7 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 6 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 5 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 4 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 3 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 2 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 1 & peptideLevel$`Observed subjects` > 1, ])

pdf(file = "Events.pdf", width = 12, height = 4.5)
Heatmap(t(peptideLevel[, c("Event")]), width = 0.5, height = unit(0.3, "cm"), col = brewer.pal(n=8, name="Set1"),
        show_row_names = F, 
        cluster_columns = F, show_column_dend = F, heatmap_legend_param = list(
          title = "Event", labels = c("Coding", "5`-UTR", "FS", "IR", "ncRNA", "3`-UTR", "IGR", 
                                      "Unknown"),
          color_bar = "discrete", direction = "horizontal", ncol = 2, nrow = 4
        )) %v%
  Heatmap(t(peptideLevel[, c("PSM1", "PSM2", "PSM3", "PSM4")]), width = 0.5, height = unit(2.2, "cm"), 
          row_title = "PSM", row_title_rot = 90, show_row_names = F, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Blues", n=4), 
          cluster_columns = F, show_column_dend = F, cluster_rows = F, heatmap_legend_param = list(
            title = "Number of PSMs", at = c(0, 1, 2, 3), 
            labels = c("0", "1", "2", ">2"), legend_height = unit(1, "cm"),
            color_bar = "discrete", nrow = 1
          )) %v%
  Heatmap(t(peptideLevel[, c("Read1","Read2","Read3","Read4")]), width = 0.5, height = unit(2.2, "cm"),
          show_row_names = F,
          row_title = "Read", row_title_rot = 90, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Reds", n=8),
          cluster_columns = F, show_column_dend = F, cluster_rows = F, heatmap_legend_param = list(
            title = "log2(number of reads+1)",legend_height = unit(1, "cm"), nrow = 1, direction = "horizontal"
          )) %v%
  Heatmap(t(peptideLevel[, c("Mutation1","Mutation2","Mutation3","Mutation4")]), width = 0.5, height = unit(2.2, "cm"),
          show_row_names = F,
          cluster_columns = F, show_column_dend = F, cluster_rows = F, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          row_title = "Mutation", row_title_rot = 90,
          col = brewer.pal(name="Blues", n=5),
          heatmap_legend_param = list(
            title = "Number of mutations", at = c(0, 1, 2, 3, 4), 
            labels = c("0", "1", "2", "3", "4"), legend_height = unit(1, "cm"), 
            color_bar = "discrete", nrow = 1
          ))  %v%
  Heatmap(t(peptideLevel[, c("AS1","AS2","AS3","AS4")]), width = 0.5, height = unit(1.3, "cm"),
          show_row_names = F,
          row_title = "AS", row_title_rot = 90, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Reds", n=3),
          cluster_columns = F, show_column_dend = F, cluster_rows = F, heatmap_legend_param = list(
            title = "Is AS", at = c(0, 1), 
            labels = c("N", "Y"), legend_height = unit(1, "cm"),
            color_bar = "discrete", nrow = 1
          )) %v%
  Heatmap(t(peptideLevel[, c("asRNA1","asRNA2","asRNA3","asRNA4")]), width = 0.5, height = unit(1.3, "cm"),
          show_row_names = F,
          row_title = "asRNA", row_title_rot = 90, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Blues", n=3),
          cluster_columns = F, show_column_dend = F, cluster_rows = F,
          heatmap_legend_param = list(
            title = "Is asRNA", at = c(0, 1), 
            labels = c("N", "Y"), legend_height = unit(1, "cm"),
            color_bar = "discrete", nrow = 1
          )) 

dev.off()

event_category <- c("ncRNA", "5`-UTR", "3`-UTR", "Coding", "FS", "IGR", "IR", "Unknown")
nPeptideLevel <- peptideLevel
nPeptideLevel$E <- ""
nPeptideLevel[nPeptideLevel$Events == "Unknown", ]$E <- "Unknown"
nPeptideLevel[nPeptideLevel$Events == "IGR", ]$E <- "IGR"
nPeptideLevel[nPeptideLevel$Events == "IR", ]$E <- "IR"
nPeptideLevel[nPeptideLevel$Events == "ncRNA", ]$E <- "ncRNA"
nPeptideLevel[nPeptideLevel$Events == "5`-UTR", ]$E <- "5`-UTR"
nPeptideLevel[nPeptideLevel$Events == "3`-UTR", ]$E <- "3`-UTR"
nPeptideLevel[nPeptideLevel$Events == "FS", ]$E <- "FS"
nPeptideLevel[nPeptideLevel$Events == "Coding", ]$E <- "Coding"

eventFeature <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(eventFeature) <- c("Event", "Category", "Count", "p value")
for(event in event_category) {
  num_ncMAPs <- nrow(nPeptideLevel)
  num_events <- nrow(nPeptideLevel[nPeptideLevel$E == event, ])
  num_asRNAs <- nrow(nPeptideLevel[nPeptideLevel$asRNA1 + nPeptideLevel$asRNA2 + nPeptideLevel$asRNA3 + nPeptideLevel$asRNA4 > 0, ])
  num_asRNAs_in_events <- nrow(nPeptideLevel[nPeptideLevel$asRNA1 + nPeptideLevel$asRNA2 + nPeptideLevel$asRNA3 + nPeptideLevel$asRNA4 > 0
                                            & nPeptideLevel$E == event, ])
  pval <- phyper(num_asRNAs_in_events-1, num_asRNAs, num_ncMAPs - num_events, num_events, lower.tail = F, log.p = FALSE) 
  eventFeature[nrow(eventFeature) + 1, ] <- c(event, "asRNA",num_asRNAs_in_events, pval)
  
  num_ASs <- nrow(nPeptideLevel[nPeptideLevel$AS1 + nPeptideLevel$AS2 + nPeptideLevel$AS3 + nPeptideLevel$AS4 > 0, ])
  num_ASs_in_events <- nrow(nPeptideLevel[nPeptideLevel$AS1 + nPeptideLevel$AS2 + nPeptideLevel$AS3 + nPeptideLevel$AS4 > 0
                                             & nPeptideLevel$E == event, ])
  pval <- phyper(num_ASs_in_events-1, num_ASs, num_ncMAPs - num_events, num_events, lower.tail = F, log.p = FALSE) 
  eventFeature[nrow(eventFeature) + 1, ] <- c(event, "AS",num_ASs_in_events, pval)
}
eventFeature$`p value` <- as.numeric(eventFeature$`p value`)
eventFeature$Count <- as.numeric(eventFeature$Count)
p.adjust(eventFeature[eventFeature$Category == "AS", ]$`p value`, method = "bonferroni")
p.adjust(eventFeature[eventFeature$Category == "asRNA", ]$`p value`, method = "bonferroni")

eventFeatureComb <- eventFeature
eventFeatureComb$Event <- paste(eventFeatureComb$Event, eventFeatureComb$Category, sep = ";")

g <- ggplot(eventFeature[eventFeature$Count != 0, ], aes(x=Count, y=Event, size = as.numeric(-log10(`p value`)), color = Category)) + 
  geom_point() +
  theme_bw() +
  ylab("") +
  #staticThemeRightBottom
  staticThemeNone
g
ggsave("ComplexEvents.png", plot = g, width = 5, height = 8, units = "in", dpi = 300)


