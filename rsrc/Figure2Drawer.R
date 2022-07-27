install.packages("latex2exp")
install.packages("stringr")                        # Install stringr package
install.packages("svglite")
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

staticThemeNone <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                             legend.position = "none")


setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified")

scores <- read.csv(file = "pXg/S4.UNMOD.PEAKS.pxg.fdr.dist", header = T, sep = "\t", as.is = as.double())
allResult <- read.csv(file = "PEAKS/S4.UNMOD.PEAKS.csv", header = T, sep = ",", as.is = as.double())

uniqueAllScores <- allResult[!duplicated(paste(allResult$Fraction, allResult$Source.File, allResult$Scan, sep = "_")), ]

allPSMs <- uniqueAllScores$Denovo.score
allPSMs <- data.frame(allPSMs)
colnames(allPSMs) <- "Score"

allScores <- allPSMs %>% dplyr::count(Score)
colnames(allScores) <- c("Score", "Count")

cPSM_targets <- scores[, c(1,2)]
cPSM_decoys <- scores[, c(1,3)]
ncPSM_targets <- scores[, c(1,4)]
ncPSM_decoys <- scores[, c(1,5)]

cPSM_targets$Class <- "Target"
cPSM_decoys$Class <- "Decoy"
ncPSM_targets$Class <- "Target"
ncPSM_decoys$Class <- "Decoy"

fieldNames <- c("Score", "Count", "Class")

colnames(cPSM_targets) <- fieldNames
colnames(cPSM_decoys) <- fieldNames
colnames(ncPSM_targets) <- fieldNames
colnames(ncPSM_decoys) <- fieldNames

c_scoreTable <- rbind(cPSM_targets, cPSM_decoys)
nc_scoreTable <- rbind(ncPSM_targets, ncPSM_decoys)
c_scoreTable$Class <- factor(c_scoreTable$Class, c("Target", "Decoy"))
nc_scoreTable$Class <- factor(nc_scoreTable$Class, c("Target", "Decoy"))

## score distribution size = 1000x800
g <- ggplot(data = allScores, aes(x=Score, y=Count)) +
  scale_fill_manual(values=c(grayC)) +
  theme_bw() +
  geom_bar(stat="identity") +
  staticThemeLeftTop +
  labs(y="PSMs", x = "ALC score")
g

ggsave("score.all.S4.png", plot = g, width = 8, height = 6, units = "in", dpi = 300)

g <- ggplot(data = c_scoreTable, aes(x=Score, y=Count, fill=Class)) +
  scale_fill_manual(values=c(blueC, redC)) +
  theme_bw() +
  geom_bar(stat="identity", position=position_dodge()) + 
  staticThemeLeftTop +
  labs(y="PSMs", x = "ALC score")

g
ggsave("score.canonical.S4.png", plot = g, width = 8, height = 6, units = "in", dpi = 300)

## Noncanonical
## FDR threshold
fdr <- 0.10
fdrScore <- 0
targetCount <- 0
decoyCount <- 0
fdrCutoff <- 0

for(alc in c(100:1)) {
  ## target count
  if(length(ncPSM_targets[ncPSM_targets$Score == alc, "Count"]) != 0){
    targetCount = targetCount + ncPSM_targets[ncPSM_targets$Score == alc, "Count"]
  }
  ## decoy count
  if(length(ncPSM_decoys[ncPSM_decoys$Score == alc, "Count"]) != 0){
    decoyCount = decoyCount + ncPSM_decoys[ncPSM_decoys$Score == alc, "Count"]
  }
  
  if(targetCount > 0) {
    cutFDR <- decoyCount/targetCount
    if(cutFDR <= fdr) {
      fdrScore <- alc
      fdrCutoff <- cutFDR
    }
  }
}
fdrScore
fdrCutoff

g <- ggplot(data = nc_scoreTable, aes(x=Score, y=Count, fill=Class)) +
  scale_fill_manual(values=c(blueC,redC)) +
  theme_bw() +
  geom_bar(stat="identity", position=position_dodge()) +
  staticThemeLeftTop +
  labs(y="PSMs", x = "ALC score")

g
ggsave("score.noncanonical.S4.png", plot = g, width = 8, height = 6, units = "in", dpi = 300)

## RNA-read distrubition
data <- read.csv(file = "pXg/S4.UNMOD.PEAKS.pxg.pval.dist", header = T, sep="\t", as.is = as.double())

data$"Peptide length" <- as.character(data$PeptideLength)
data$Mock <- log2(data$Mock+1)
data$Experiment <- log2(data$Experiment+1)
data$`Peptide length` <- factor(data$PeptideLength, levels = c('8','9','10','11','12','13', '14', '15'))
data <- data[data$ReadCount <= 30, ]


g <- ggplot(data = data, aes(x=ReadCount, y=Mock, fill=`Peptide length`)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_line(size = 1, aes(color=`Peptide length`)) +
  xlab("Read") +
  ylab(TeX("$Log_{2}$(number of peptides + 1)")) +
  scale_x_continuous(breaks=seq(from=0, to=30, by = 1)) +
  staticThemeRightTop +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in")) +
  geom_point()
g
ggsave("read.S4.png", plot = g, width = 12, height = 6, units = "in", dpi = 300)

## PPM-tolerance
ppmData <-  read_excel(path = "pXgResults.xlsx", sheet = "B-LCL4")
unIdedData <- read.csv(file = "PEAKS/S4.UNMOD.PEAKS.csv.unided", header = T, sep = ",", as.is = as.double())

ppmData <- ppmData[!duplicated(paste(ppmData$Fraction, ppmData$`Source File`, ppmData$Scan, sep = "_")), ]
uniqueUnIdedData <- unIdedData[!duplicated(paste(unIdedData$Fraction, unIdedData$Source.File, unIdedData$Scan, sep = "_")), ]


subjectNames <- c("B-LCL subject1", "B-LCL subject2", "B-LCL subject3", "B-LCL subject4")
subjectFDRScore <- c(75, 78, 90, 84)

unIdedData <- uniqueUnIdedData$ppm
unIdedData <- data.frame(unIdedData)

cPPM <- ppmData[ppmData$IsCanonical == T, ]
cPPM <- cPPM$ppm
ncPPM <- ppmData[ppmData$IsCanonical == F, ]
ncFDRPPM <- ncPPM[ncPPM$`Denovo score` >= 84, ]$ppm

cPPM <- data.frame(cPPM)
ncFDRPPM <- data.frame(ncFDRPPM)

colnames(unIdedData) <- "PPM tolerance"
colnames(cPPM) <- "PPM tolerance"
colnames(ncFDRPPM) <- "PPM tolerance"

unIdedData$Class <- "Unidentified PSM"
cPPM$Class <- "cPSM"
ncFDRPPM$Class <- "ncPSM"

data <- rbind(unIdedData, cPPM)
data <- rbind(data, ncFDRPPM)

data$Class <- factor(data$Class, c("Unidentified PSM", "cPSM", "ncPSM"))

g <- ggplot(data=data) +
  geom_density(aes(x=`PPM tolerance`, colour = Class), size = 1) +
  scale_color_manual(values=c(grayC, blueC, redC)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(from=0, to=1, by = 0.1), 
                     labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE),
                     limits=c(0,0.5)) +
  staticThemeRightTop +
  labs(y="Density", x = "PPM error")

g
ggsave("ppm.S4.png", plot = g, width = 8, height = 6, units = "in", dpi = 300)

ks.test(data[data$Class == "Unidentified PSM", ]$`PPM tolerance`, data[data$Class == "cPSM", ]$`PPM tolerance`)
ks.test(data[data$Class == "Unidentified PSM", ]$`PPM tolerance`, data[data$Class == "ncPSM", ]$`PPM tolerance`)
ks.test(data[data$Class == "cPSM", ]$`PPM tolerance`, data[data$Class == "ncPSM", ]$`PPM tolerance`)

## NetMHCpan BA distribution
setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/")
dataS1 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL1_FDR10")
dataS2 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL2_FDR10")
dataS3 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL3_FDR10")
dataS4 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL4_FDR10")

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
## select unmodi
allData$Length <- as.character(allData$Length)
allData$Length <- factor(allData$Length, levels=c("8", "9", "10", "11", "12", "13", "14", "15"))

allData$BestScore <- log2(allData$BestScore+1.5)
allData$Class <- "Canonical"
allData[allData$IsCanonical == "TRUE", ]$Class <- "Canonical"
allData[allData$IsCanonical == "FALSE", ]$Class <- "Noncanonical"

subData <- allData[allData$Sample == "B-LCL subject 4", ]
nrow(subData[subData$Length == 8 & subData$IsCanonical == "FALSE", ])
mapPlot <- ggplot(data=subData, aes(x=Length, y=BestScore, fill=Length)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  scale_y_continuous(breaks=seq(from=0, to=7, by= 1), limits = c(0,7)) +
  staticThemeNone +
  labs(y= TeX("$Log_{2}$(1.5+%Rank)"), x = "Peptide length") +
  facet_grid(. ~ Class) +
  geom_hline(yintercept=log2(1.5+2.0), linetype="dashed", color = "blue", size=1) +
  geom_hline(yintercept=log2(1.5+0.5), linetype="dashed", color = "red", size=1)

mapPlot
ggsave("BA.S4.png", plot = mapPlot, width = 14, height = 8, units = "in", dpi = 300)

nrow(subData[subData$BestScore < log2(3.5) & subData$Length == 8, ])
nrow(subData)


## Peptide length distribution
nrow(subData[subData$Class == "Noncanonical", ])
nrow(subData[subData$Class == "Canonical", ])

## select length
## discard over length 13 because of small
## Of course, it slightly twists the ratio
subData <- subData[subData$Length < 14, ]

len = 9
matrix <- matrix(c(nrow(subData[subData$IsCanonical == "TRUE", ]) - nrow(subData[subData$Length == len & subData$IsCanonical == "TRUE", ]), 
                   nrow(subData[subData$Length == len & subData$IsCanonical == "TRUE", ]),
                   nrow(subData[subData$IsCanonical == "FALSE", ]) - nrow(subData[subData$Length == len & subData$IsCanonical == "FALSE", ]), 
                   nrow(subData[subData$Length == len & subData$IsCanonical == "FALSE", ])), 
                 nrow = 2, ncol = 2)
matrix <- t(matrix)
fisher.test(matrix)

topRankedPlot <- ggplot(data=subData, aes(x=Length, fill=Class)) +
  scale_fill_manual(values=c(blueC, redC)) +
  theme_bw() +
  geom_bar(aes(y= (..prop..), x=Length), position=position_dodge(), width = 0.5) +
  scale_x_continuous(breaks = seq(from=8, to=13, by=1)) +
  scale_y_continuous(breaks=seq(from=0, to=1, by = 0.1), 
                     labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE),
                     limits=c(0,0.7)) +
  theme(text = element_text(size=25), 
        plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 30, color = "darkblue")) +
  labs(y="Proportion", x = "Peptide length") +
 annotate("text", label = "*", x =8, y = 0.22, size = 8, family="serif", colour = "black") +
 annotate("text", label = "*", x =9, y = 0.70, size = 8, family="serif", colour = "black") +
#annotate("text", label = "*", x =10, y = 0.2, size = 8, family="serif", colour = "black") +
#annotate("text", label = "*", x =11, y = 0.1, size = 8, family="serif", colour = "black") +
#annotate("text", label = "*", x =12, y = 0.05, size = 8, family="serif", colour = "black") +
#annotate("text", label = "*", x =13, y = 0.04, size = 8, family="serif", colour = "black") +
#annotate("text", label = "*", x =14, y = 0.03, size = 8, family="serif", colour = "black") +
 geom_segment(aes(x = 7.75, xend =8.25, y = 0.21, yend = 0.21), color = "Black", size = 1) +
 geom_segment(aes(x = 8.75, xend =9.25, y = 0.69, yend = 0.69), color = "Black", size = 1) 
#geom_segment(aes(x = 9.75, xend =10.25, y = 0.19, yend = 0.19), color = "Black", size = 1) +
#geom_segment(aes(x = 10.75, xend =11.25, y = 0.09, yend = 0.09), color = "Black", size = 1) +
#geom_segment(aes(x = 11.75, xend =12.25, y = 0.04, yend = 0.04), color = "Black", size = 1) +
#geom_segment(aes(x = 12.75, xend =13.25, y = 0.03, yend = 0.03), color = "Black", size = 1) +
#geom_segment(aes(x = 13.75, xend =14.25, y = 0.02, yend = 0.02), color = "Black", size = 1)

topRankedPlot
ggsave("BA.Length.png", plot = topRankedPlot, width = 10, height = 8, units = "in", dpi = 300)



## HLA length distribution
subData <- allData[allData$BestScore < 2, ]
subData$Class <- "cMAP"
subData[subData$IsCanonical == "FALSE", ]$Class <- "ncMAP"
#subData$Length <- as.character(subData$Length)
nrow(subData[subData$Class == "ncMAP", ])
nrow(subData[subData$Class == "cMAP", ])
## select length
## discard over length 13 because of small
## Of course, it slightly twists the ratio
subData <- subData[subData$Length < 14, ]
subData <- subData[str_detect(subData$BestType, "HLA-A", negate = F), ]

topRankedPlot <- ggplot(data=subData, aes(x=Length, fill=BestType)) +
  scale_fill_brewer(palette="Set1", name ="HLA type") +
  theme_bw() +
  geom_bar(aes(y= (..prop..), x=Length), position=position_dodge2(preserve = "single", padding = 0), width = 0.5) +
  scale_x_continuous(breaks = seq(from=8, to=13, by=1)) +
  scale_y_continuous(breaks=seq(from=0, to=1, by = 0.1), 
                     labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE),
                     limits=c(0,0.7)) +
  theme(text = element_text(size=25), 
        plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 30, color = "darkblue"))+
  labs(y="Proportion", x = "Peptide length") 

topRankedPlot
ggsave("BA.HLA-A.Length.png", plot = topRankedPlot, width = 10, height = 8, units = "in", dpi = 300)

## HLA type
## unique selection
subData <- allData[allData$BestScore < 2 & allData$EventCount == 1 & allData$GeneIDCount <= 1 & allData$GenomicLociCount <= 1, ]


subData$Class <- "cMAP"
subData[subData$IsCanonical == "FALSE", ]$Class <- "ncMAP"
#subData$Length <- as.character(subData$Length)
nrow(subData[subData$Class == "ncMAP", ])
nrow(subData[subData$Class == "cMAP", ])
subData <- subData[subData$Class == "ncMAP",]

subData$Sample <- with(subData, factor(subData$Sample, levels=c("B-LCL subject 1", "B-LCL subject 2", "B-LCL subject 3", "B-LCL subject 4")))

nrow(subData)

topRankedPlot <- ggplot(data=subData, aes(x=BestType, fill=Sample)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  geom_bar(aes(y= (..count..), x=BestType), position=position_dodge2(preserve = "single", padding = 0, reverse = F), width = 1) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 90, size = 15)) +
  labs(y="ncMAPs", x = "\nHLA type") + 
  ggtitle("Unique ncMAPs according to HLA type")

topRankedPlot
ggsave("BA.HLA.uncMAPs.png", plot = topRankedPlot, width = 10, height = 8, units = "in", dpi = 300)

## Event enumerate
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
subData$Class <- "Canonical"
subData[subData$IsCanonical == "FALSE", ]$Class <- "Noncanonical"
#subData$Length <- as.character(subData$Length)
nrow(subData[subData$Class == "Noncanonical", ])
nrow(subData[subData$Class == "Canonical", ])

nrow(subData[subData$Sample == "B-LCL subject 1" & subData$Class == "Canonical", ])
nrow(subData[subData$Sample == "B-LCL subject 2" & subData$Class == "Canonical", ])
nrow(subData[subData$Sample == "B-LCL subject 3" & subData$Class == "Canonical", ])
nrow(subData[subData$Sample == "B-LCL subject 4" & subData$Class == "Canonical", ])

nrow(subData[subData$Sample == "B-LCL subject 1" & subData$Class == "Noncanonical", ])
nrow(subData[subData$Sample == "B-LCL subject 2" & subData$Class == "Noncanonical", ])
nrow(subData[subData$Sample == "B-LCL subject 3" & subData$Class == "Noncanonical", ])
nrow(subData[subData$Sample == "B-LCL subject 4" & subData$Class == "Noncanonical", ])


subData <- subData[subData$Class == "Noncanonical", ]
subData$Sample <- factor(subData$Sample, levels = c("B-LCL subject 1", "B-LCL subject 2", "B-LCL subject 3", "B-LCL subject 4"))
### Unique selection
subDataUnique <- subData[subData$EventCount == 1, ]
subDataUnique <- subDataUnique[subDataUnique$GeneIDCount <= 1, ]
subDataUnique <- subDataUnique[subDataUnique$GenomicLociCount <= 1, ]
subDataUnique <- subDataUnique[!duplicated(subDataUnique[,c('InferredPeptide')]), ]
nrow(subDataUnique)
### shared
subDataShared <- subData[subData$EventCount > 1 | subData$GeneIDCount > 1 | subData$GenomicLociCount > 1, ]
subDataShared <- subDataShared[!duplicated(subDataShared[,c('InferredPeptide')]), ]
nrow(subDataShared)

byEventCount <- subDataShared[subDataShared$EventCount > 1, ]
byGeneIDCount <- subDataShared[subDataShared$GeneIDCount > 1, ]
byGenomicLociCount <- subDataShared[subDataShared$GenomicLociCount > 1, ]

nrow(subData)
nrow(subDataUnique)
nrow(subDataShared)

nrow(byEventCount)
nrow(byGeneIDCount)
nrow(byGenomicLociCount)

## shared venn count
inter3 <- subDataShared[subDataShared$EventCount > 1 & subDataShared$GeneIDCount > 1 & subDataShared$GenomicLociCount > 1,]
interEventLoci <- subDataShared[subDataShared$EventCount > 1 & subDataShared$GenomicLociCount > 1,]
interEventGene <- subDataShared[subDataShared$EventCount > 1 & subDataShared$GeneIDCount > 1,]
interLociGene <- subDataShared[subDataShared$GeneIDCount > 1 & subDataShared$GenomicLociCount > 1,]
eventOnly <- subDataShared[subDataShared$EventCount > 1 & subDataShared$GeneIDCount <= 1 & subDataShared$GenomicLociCount <= 1,]

nrow(eventOnly)
nrow(inter3)
nrow(interEventLoci)
nrow(interEventGene)
nrow(interLociGene)

subDataUnique[subDataUnique$Mutations != "-", ]$Events <- paste(subDataUnique[subDataUnique$Mutations != "-", ]$Events, "mutations", sep=";")
subDataUnique$Events <- factor(subDataUnique$Events, 
                               levels = c("5'-UTR;sense", "proteincoding;sense;mutations","noncoding;sense", "frameshift;sense",
                                          "intronic;sense", "3'-UTR;sense","intergenic;sense", "5'-UTR;sense;mutations",
                                          "proteincoding;sense;alternativesplicing", "3'-UTR;antisense", "intronic;antisense",
                                          "proteincoding;antisense;mutations", "noncoding;sense;mutations", "5'-UTR;antisense",
                                          "intergenic;sense;mutations","noncoding;antisense","unknown"))
nrow(subDataUnique[subDataUnique$Events == "proteincoding;sense;mutations" ,])
nrow(subDataUnique[subDataUnique$Events == "5UTR;sense" ,])
nrow(subDataUnique[subDataUnique$Events == "5UTR;sense;mutations" ,])
nrow(subDataUnique[subDataUnique$Events == "noncoding;sense" ,])
nrow(subDataUnique[subDataUnique$Events == "noncoding;sense;mutations" ,])
nrow(subDataUnique[subDataUnique$Events == "frameshift;sense" ,])
nrow(subDataUnique[subDataUnique$Events == "frameshift;sense;mutations" ,])

nrow(subDataUnique[subDataUnique$Sample == "B-LCL subject 1", ])
nrow(subDataUnique[subDataUnique$Sample == "B-LCL subject 2", ])
nrow(subDataUnique[subDataUnique$Sample == "B-LCL subject 3", ])
nrow(subDataUnique[subDataUnique$Sample == "B-LCL subject 4", ])
## mutated
topRankedPlot <- ggplot(data=subDataUnique, aes(x=Events, fill=Sample)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  geom_bar(aes(x= (..count..), y=reorder(Events, desc(Events))), position=position_dodge2(preserve = "single", padding = 0, reverse = T), width = 1) +
  theme(text = element_text(size=20), 
        axis.text.y = element_text(size = 20, colour = "black")) +
  ggtitle("Category of unique ncMAPs") +
  labs(y="Event", x = "Unique ncMAPs") 

topRankedPlot
ggsave("BA.ncMAPs.Events.png", plot = topRankedPlot, width = 10, height = 16, units = "in", dpi = 600)
dev.off()

topRankedPlot <- ggplot(data=subDataUnique, aes(x=BestType, fill=Events)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  geom_bar(aes(x= (..count..), y=BestType), position=position_dodge2(preserve = "single", padding = 0, reverse = T), width = 1) +
  theme(text = element_text(size=20), 
        axis.text.y = element_text(size = 15)) +
  ggtitle("Category of unique ncMAPs") +
  labs(y="Event", x = "ncMAPs") 

topRankedPlot

nrow(subDataUnique[subDataUnique$BestType == "HLA-A03:01", ])
nrow(subDataUnique[subDataUnique$BestType == "HLA-A03:01" & subDataUnique$Events == "5UTR;sense", ])
nrow(subDataUnique[subDataUnique$BestType == "HLA-A03:01" & subDataUnique$Events == "frameshift;sense", ])
nrow(subDataUnique[subDataUnique$BestType == "HLA-A03:01" & subDataUnique$Events == "proteincoding;sense;mutations", ])


nrow(subDataUnique[subDataUnique$BestType == "HLA-B44:03", ])
nrow(subDataUnique[subDataUnique$BestType == "HLA-B44:03" & subDataUnique$Events == "5UTR;sense", ])
nrow(subDataUnique[subDataUnique$BestType == "HLA-B44:03" & subDataUnique$Events == "frameshift;sense", ])
nrow(subDataUnique[subDataUnique$BestType == "HLA-B44:03" & subDataUnique$Events == "proteincoding;sense;mutations", ])

nrow(subDataUnique[subDataUnique$BestType == "HLA-C16:01", ])
nrow(subDataUnique[subDataUnique$BestType == "HLA-C16:01" & subDataUnique$Events == "5UTR;sense", ])
nrow(subDataUnique[subDataUnique$BestType == "HLA-C16:01" & subDataUnique$Events == "frameshift;sense", ])
nrow(subDataUnique[subDataUnique$BestType == "HLA-C16:01" & subDataUnique$Events == "proteincoding;sense;mutations", ])


## Rank
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

allData <- rbind(dataS1, dataS2, dataS3, dataS4)

for(rank in c(3:10)) {
  allData[allData$Rank == rank, ]$Rank <- 3
}

for(rank in c(1:3)) {
  rankRatio <- nrow(allData[allData$Rank == rank & allData$`MHC-I` != "NB", ]) / nrow(allData[allData$Rank == rank, ])  
  print(paste(rank, rankRatio, sep="_"))
}

allData$Rank <- as.character(allData$Rank)
allData[allData$Rank == 3, ]$Rank <- "Others"
allData$`MHC-I` <- factor(allData$`MHC-I`, levels = c("SB", "WB", "NB"))

nrow(allData[allData$Rank == 1, ])
nrow(allData[allData$Rank == 2, ])
nrow(allData[allData$Rank == "Others", ])
nrow(allData[allData$Rank == 1, ]) / nrow(allData)  
nrow(allData[allData$Rank == 2, ]) / nrow(allData)  
nrow(allData[allData$Rank == "Others", ]) / nrow(allData)  

topRankedPlot <- ggplot(data=allData, aes(x=Rank, fill=`MHC-I`)) +
  scale_fill_manual("Class",values = c("SB" = blueC, "WB" = greenC, "NB" = redC))+
  theme_bw() +
  geom_bar(aes(y= (..count..), x=Rank)) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 15)) +
  labs(y="PSMs", x = "Candidate rank") + 
  ggtitle("Binding prediction according to candidate rank")

topRankedPlot
ggsave("BA.Ranks.png", plot = topRankedPlot, width = 10, height = 8, units = "in", dpi = 300)

## IEDB
IEDB <- read.csv(file = "/Users/gistar/projects/pXg/IEDB/HLA-A0101.csv", header = T, sep = ",", as.is = as.double())
IEDB$Epitope.2
IEDB$Epitope.10


## chromosomes
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

allData <- rbind(dataS1, dataS2, dataS3, dataS4)

subData <- allData[allData$BestScore < 2, ]
subData$Class <- "Canonical"
subData[subData$IsCanonical == "FALSE", ]$Class <- "Noncanonical"
#subData$Length <- as.character(subData$Length)
nrow(subData[subData$Class == "Noncanonical", ])
nrow(subData[subData$Class == "Canonical", ])

subData <- subData[subData$Class == "Noncanonical", ]
nrow(subData[subData$Sample == "B-LCL subject 1", ])
nrow(subData[subData$Sample == "B-LCL subject 2", ])
nrow(subData[subData$Sample == "B-LCL subject 3", ])
nrow(subData[subData$Sample == "B-LCL subject 4", ])

nrow(subData)

subDataUnique <- subData[subData$EventCount == 1, ]
subDataUnique <- subDataUnique[subDataUnique$GeneIDCount <= 1, ]
subDataUnique <- subDataUnique[subDataUnique$GenomicLociCount <= 1, ]
subDataUnique <- subDataUnique[!duplicated(subDataUnique[,c('InferredPeptide')]), ]
nrow(subDataUnique)
subDataUnique$Chromosome <- str_split(subDataUnique$GenomicLoci, pattern = ":", simplify = T)[,1]
subDataUnique$Region <- "unknown"
subDataUnique[subDataUnique$Events == "3UTR;antisense" |subDataUnique$Events == "3UTR;sense", ]$Region <- "3'-UTR"
subDataUnique[subDataUnique$Events == "5UTR;antisense" |subDataUnique$Events == "5UTR;sense", ]$Region <- "5'-UTR"
subDataUnique[subDataUnique$Events == "frameshift;antisense" |subDataUnique$Events == "frameshift;sense", ]$Region <- "frameshift"
subDataUnique[subDataUnique$Events == "intron;antisense" |subDataUnique$Events == "intron;sense", ]$Region <- "intron"
subDataUnique[subDataUnique$Events == "intergenic;antisense" |subDataUnique$Events == "intergenic;sense", ]$Region <- "intergenic"
subDataUnique[subDataUnique$Events == "proteincoding;antisense" |subDataUnique$Events == "proteincoding;sense"|subDataUnique$Events == "proteincoding;sense;alternativesplicing", ]$Region <- "proteincoding"
subDataUnique[subDataUnique$Events == "noncoding;antisense" |subDataUnique$Events == "noncoding;sense", ]$Region <- "noncoding"

subDataUnique$Region <- factor(subDataUnique$Region, levels=c("5'-UTR", "proteincoding", "noncoding", "frameshift",
                                                              "intron", "3'-UTR", "intergenic", "unknown"))

subDataUnique$Chromosome <- factor(subDataUnique$Chromosome, levels = c("chr2", "chr1", "chr6", "chr8",
                                                                        "chr3", "chr11", "chr5", "chr19",
                                                                        "chr7", "chr12", "chr13", "chr14",
                                                                        "chr4", "chr16", "chr17", "chr18", "chr20", "chr22",
                                                                        "chrx", "chr10", "chr9", "chr15",
                                                                        "chr21", "-"))

topRankedPlot <- ggplot(data=subDataUnique, aes(x=Chromosome, fill = Region)) +
  scale_fill_brewer(palette = "Set1")+
  theme_bw() +
  geom_bar(aes(y= (..count..), x=Chromosome)) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 15)) +
  labs(y="ncMAPs", x = "Chromosome") + 
  ggtitle("")

topRankedPlot
ggsave("unique.ncMAP.chromosome.svg", plot = topRankedPlot, width = 16, height = 8, units = "in", dpi = 300)

theme_get()
