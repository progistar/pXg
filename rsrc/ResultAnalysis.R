library(ggplot2)
library(RColorBrewer)

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


setwd("C:/Users/progi/Desktop/Projects/pXg/Laumont_NatCommun2016/Results")

pResult <- read.csv(file = "subjectM.5ppm.002.peaks", header = T, sep = "\t", as.is = as.double())
result <- read.csv(file = "subjectM.5ppm.002.rep1.netMHCpan.laumont.rt.pXg", header = T, sep = "\t", as.is = as.double())

pResult <- data.frame(pResult)
result <- data.frame(result)

pUniRes <- pResult[!duplicated(pResult[,c('Fraction','Source.File','Scan')]),]
pUniRes <- pUniRes[pUniRes$Length > 7 & pUniRes$Length < 16, ]
pUniRes$Class <- "All PSMs"



alUniRes <- result[!duplicated(result[,c('Fraction','Source.File','Scan')]), ]
alUniRes <- alUniRes[alUniRes$Length > 7 & alUniRes$Length < 16, ]
alUniRes$Class <- "All PSMs"

pcUniRes <- alUniRes[alUniRes$Events == "proteincoding;sense", ]
pcUniRes$Class <- "Protein coding PSMs"

nvUniRes <- alUniRes[alUniRes$Events != "proteincoding;sense", ]
nvUniRes <- nvUniRes[nvUniRes$FastaIDCount == 0, ]
nvUniRes$Class <- "Novel PSMs"

## PPM TOL
subData <- pResult[, c("ppm", "Length" )]
subData$Class <- "All"

tmpData <- pcUniRes[pcUniRes$BestScore < 2 & pcUniRes$ALC.... >= 50, c("ppm", "Length")]
tmpData$Class <- "Protein coding"

subData <- rbind(subData, tmpData)

tmpData <- nvUniRes[nvUniRes$BestScore < 2 & nvUniRes$ALC.... >= 50, c("ppm", "Length")]
tmpData$Class <- "Novel"

subData <- rbind(subData, tmpData)

subData$ppm <- as.numeric(as.character(subData$ppm))

t.test(subData[subData$Class=="Protein coding",]$ppm, subData[subData$Class=="All",]$ppm)
t.test(subData[subData$Class=="Novel",]$ppm, subData[subData$Class=="All",]$ppm)
t.test(subData[subData$Class=="Protein coding",]$ppm, subData[subData$Class=="Novel",]$ppm)

ppmPlot <- ggplot(data=subData) +
  geom_density(aes(x=ppm, colour = Class), size = 1) +
  scale_color_manual(values=c(greenC, redC, blueC)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(from=0, to=1, by = 0.1), 
                     labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE),
                     limits=c(0,0.3)) +
  theme(text = element_text(size=25), 
        plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 30, color = "darkblue")) +
  labs(y="Density", x = "PPM error")

ppmPlot

## Optionally
subData <- rbind(pUniRes[, c("ALC....", "Class")], pcUniRes[, c("ALC....", "Class")], nvUniRes[, c("ALC....", "Class")])
subData$Class <- factor(subData$Class, c("All PSMs", "Protein coding PSMs", "Novel PSMs"))

# Score violin plot
## ALC Scores in Subject1-1
## for RNA-seq replicate 1 
nrow(pUniRes)
nrow(nvUniRes)
nrow(pcUniRes)

alcScoreDist <- 
  ggplot(data=subData, aes(x=Class, y=ALC...., fill = as.character(Class))) +
  theme_bw() +
  scale_fill_manual(values=c(greenC, redC, blueC)) +
  geom_violin(position = "dodge") +
  theme(text = element_text(size=25), legend.position = "none") +
  labs(x="", y = "ALC score") +
  scale_y_continuous(breaks=seq(from=0, to=100, by= 10)) +
  annotate("text", label = nrow(pUniRes), x =1, y = -3, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(pcUniRes), x =2, y = -3, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(nvUniRes), x =3, y = -3, size = 8, family="serif", colour = "black")
  

alcScoreDist

# PSMs according to Class
pcUniRes$EventCount

allPC <- pcUniRes[!duplicated(pcUniRes[,c('InferredPeptide')]), 
                  c("Length", "Class", "BestType", "BestScore", "GenomicLociCount", "ALC....", "EventCount", "Laumont" )]
allPC$Event <- "Protein coding"
allPC$Filter <- "pXg"

mapPC <- allPC[allPC$BestScore < 2, ]
mapPC$Filter <- "MAPs"

mapHighPC <- mapPC[mapPC$ALC.... >= 50, ]
mapHighPC$Filter <- "ALC >= 50"

allNV <- nvUniRes[!duplicated(nvUniRes[,c('InferredPeptide')]), 
                  c("Length", "Class", "BestType", "BestScore", "GenomicLociCount", "ALC....", "EventCount", "Laumont" )]
allNV$Event <- "Novel"
allNV$Filter <- "pXg"

mapNV <- allNV[allNV$BestScore < 2, ]
mapNV$Filter <- "MAPs"

mapHighNV <- mapNV[mapNV$ALC.... >= 50, ]
mapHighNV$Filter <- "ALC >= 50"

subData <- rbind(allPC, mapPC, mapHighPC, allNV, mapNV, mapHighNV)

subData$Event <- factor(subData$Event, levels=c("Protein coding", "Novel"))
subData$Filter <- factor(subData$Filter, levels=c("pXg", "MAPs", "ALC >= 50"))

topRankedPlot <- ggplot(data=subData, aes(x=Filter, fill=Event)) +
  scale_fill_manual(values=c(blueC, redC)) +
  theme_bw() +
  geom_bar(aes(y= (..count..), x=Filter), position=position_dodge(), width = 0.5) +
  #geom_text(stat='count', aes(label = scales::comma(..count.., accuracy = 1)), vjust=-0.5, hjust=c(1.0,-0.5,1.0,-0.5,1.0,-0.5,1.0,-0.5), family = "serif", size = 6) +
  geom_text(stat='count', aes(label = ..count..), 
            vjust=-0.5, hjust=c(1.0,0,1.0,0,1.0,0), family = "serif", size = 6) +
  theme(text = element_text(size=25)) +
  labs(y="Peptides", x = "")

topRankedPlot



## MAP score plot
subData <- allNV
subData$Length <- as.character(subData$Length)
subData$Length <- factor(subData$Length, levels=c("8", "9", "10", "11", "12", "13", "14", "15"))

subData$BestScore <- log2(subData$BestScore+1.5)
mapPlot <- ggplot(data=subData, aes(x=Length, y=BestScore, fill=Length)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  scale_y_continuous(breaks=seq(from=0, to=7, by= 1)) +
  labs(y="NetMHCpan\nlog2-transform (1.5+%Rank)", x = "Peptide length") +
  annotate("text", label = nrow(subData[subData$Length == 8,]), x =1, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 9,]), x =2, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 10,]), x =3, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 11,]), x =4, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 12,]), x =5, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 13,]), x =6, y = 7, size = 8, family="serif", colour = "black") +
  geom_hline(yintercept=log2(3.5), linetype="dashed", color = "blue", size=1) +
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=1)
  #annotate("text", label = nrow(subData[subData$Length == 14,]), x =7, y = 7, size = 8, family="serif", colour = "black") +
  #annotate("text", label = nrow(subData[subData$Length == 15,]), x =8, y = 7, size = 8, family="serif", colour = "black")

mapPlot

# Score violin plot
## ALC Scores in Subject1-1 in binding only
## for RNA-seq replicate 1 

mapHighUniPC <- pcUniRes[pcUniRes$BestScore < 2 & pcUniRes$ALC.... >= 50, ]
mapHighUniPC$Class <- "Protein coding"

mapHighUniNV <- nvUniRes[nvUniRes$BestScore < 2 & nvUniRes$ALC.... >= 50, ]
mapHighUniNV$Class <- "Novel"


subData <- rbind(mapHighUniNV, mapHighUniPC)
subData$Class <- factor(subData$Class, levels=c("Protein coding", "Novel"))

tTest <- t.test(subData[subData$Class == "Novel", ]$ALC...., subData[subData$Class == "Protein coding", ]$ALC....)

alcScoreDist <- 
  ggplot(data=subData, aes(x=Class, y=ALC...., fill = as.character(Class))) +
  theme_bw() +
  scale_fill_manual(values=c(redC, blueC)) +
  geom_violin(position = "dodge") +
  theme(text = element_text(size=25), plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 30, color = "darkblue"),
        legend.position = "none") +
  labs(x="", y = "ALC score") +
  scale_y_continuous(breaks=seq(from=50, to=100, by= 10)) +
  annotate("text", label = nrow(subData[subData$Class == "Protein coding", ]), x =1, y = 45, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Class == "Novel", ]), x =2, y = 45, size = 8, family="serif", colour = "black") +
  annotate("text", label = paste("p-value =",format(tTest$p.value, decimal.mark=".", digits = 3),sep = " "), x =1.5, y = 105, size = 8, family="serif", colour = "black")


alcScoreDist



## RT Plot
subData <- mapHighUniNV
subData <- mapHighUniPC
nrow(subData)

subData$deeplcRT <- as.double(as.character(subData$deeplcRT))

corr <- cor.test(subData$RT, subData$deeplcRT)

rtPlot <- ggplot(data=subData, aes(x=RT, y=deeplcRT)) +
  geom_point(colour=blueC,size=1, alpha =0.3) +
  theme_bw() +
  theme(text = element_text(size=25)) +
  labs(y="Predicted retention time", x = "Observed retention time") +
  geom_smooth(method=lm, se=FALSE, formula = y ~ x, colour=redC, alpha = 0.5) +
  xlim(c(10,100)) +
  ylim(c(10,100)) +
  annotate("text", label = paste("r =",format(corr$estimate[[1]], decimal.mark=".", digits = 3),sep = " "), 
           x =75, y = 90, size = 8, family="serif", colour = "black", fontface = 'italic')

rtPlot


## Binding affinity length distribution
mapHighUniPC <- pcUniRes[pcUniRes$BestScore < 2 & pcUniRes$ALC.... >= 50, ]
mapHighUniPC <- mapHighUniPC[!duplicated(mapHighUniPC[,c('InferredPeptide')]), ]
mapHighUniPC$Event <- "Protein coding"
mapHighUniPC <- mapHighUniPC[, c("Event", "Length")]

mapHighUniNV <- nvUniRes[nvUniRes$BestScore < 2 & nvUniRes$ALC.... >= 50, ]
mapHighUniNV <- mapHighUniNV[!duplicated(mapHighUniNV[,c('InferredPeptide')]), ]
mapHighUniNV$Event <- "Novel"
mapHighUniNV <- mapHighUniNV[, c("Event", "Length")]

subData <- rbind(mapHighUniPC, mapHighUniNV)
subData$Event <- factor(subData$Event, levels=c("Protein coding", "Novel"))


## Fisher exact test
len = 13
matrix <- matrix(c(nrow(subData[subData$Event == "Protein coding", ]), nrow(subData[subData$Length == len & subData$Event == "Protein coding", ]),
                   nrow(subData[subData$Event == "Novel", ]), nrow(subData[subData$Length == len & subData$Event == "Novel", ])), 
                 nrow = 2, ncol = 2)
matrix <- t(matrix)
fisher.test(matrix)

topRankedPlot <- ggplot(data=subData, aes(x=Length, fill=Event)) +
  scale_fill_manual(values=c(blueC, redC)) +
  theme_bw() +
  geom_bar(aes(y= (..prop..), x=Length), position=position_dodge(), width = 0.5) +
  scale_x_continuous(breaks = seq(from=8, to=14, by=1)) +
  scale_y_continuous(breaks=seq(from=0, to=1, by = 0.1), 
                     labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE),
                     limits=c(0,0.7)) +
  theme(text = element_text(size=25), 
        plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 30, color = "darkblue")) +
  labs(y="Proportion", x = "Peptide length") +
  annotate("text", label = "*", x =8, y = 0.25, size = 8, family="serif", colour = "black") +
  annotate("text", label = "*", x =9, y = 0.62, size = 8, family="serif", colour = "black") +
  annotate("text", label = "*", x =10, y = 0.2, size = 8, family="serif", colour = "black") +
  annotate("text", label = "*", x =11, y = 0.1, size = 8, family="serif", colour = "black") +
  annotate("text", label = "*", x =12, y = 0.05, size = 8, family="serif", colour = "black") +
  annotate("text", label = "*", x =13, y = 0.04, size = 8, family="serif", colour = "black") +
  annotate("text", label = "*", x =14, y = 0.03, size = 8, family="serif", colour = "black") +
  geom_segment(aes(x = 7.75, xend =8.25, y = 0.24, yend = 0.24), color = "Black", size = 1) +
  geom_segment(aes(x = 8.75, xend =9.25, y = 0.61, yend = 0.61), color = "Black", size = 1) +
  geom_segment(aes(x = 9.75, xend =10.25, y = 0.19, yend = 0.19), color = "Black", size = 1) +
  geom_segment(aes(x = 10.75, xend =11.25, y = 0.09, yend = 0.09), color = "Black", size = 1) +
  geom_segment(aes(x = 11.75, xend =12.25, y = 0.04, yend = 0.04), color = "Black", size = 1) +
  geom_segment(aes(x = 12.75, xend =13.25, y = 0.03, yend = 0.03), color = "Black", size = 1) +
  geom_segment(aes(x = 13.75, xend =14.25, y = 0.02, yend = 0.02), color = "Black", size = 1)

topRankedPlot


## Event

###
library(plyr)
subData <- pcUniRes[pcUniRes$BestScore < 2 & pcUniRes$ALC.... >= 50, ]
subData <- subData[!duplicated(subData[,c('InferredPeptide')]),]
subData <- subData[subData$GeneIDCount == 1, ]

uniquePCGenes <- ddply(subData, .(GeneIDs, GeneNames), nrow)

uniquePCGenes <- uniquePCGenes[order(uniquePCGenes$V1, decreasing = TRUE), ]


subData <- nvUniRes[nvUniRes$BestScore < 2 & nvUniRes$ALC.... >= 50, ]
subData <- subData[!duplicated(subData[,c('InferredPeptide')]),]
subData <- subData[subData$GeneIDCount == 1 & subData$EventCount == 1, ]

uniqueNVGenes <- ddply(subData, .(GeneIDs, GeneNames), nrow)

uniqueNVGenes <- uniqueNVGenes[order(uniqueNVGenes$V1, decreasing = TRUE), ]

