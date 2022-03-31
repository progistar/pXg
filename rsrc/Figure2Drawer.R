## Histogram for score distribution
library(ggplot2)
library(RColorBrewer)
library("dplyr")
library(readxl)

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


setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/5.withCalibrationAddScanNumWithoutDeami")

scores <- read.csv(file = "pXg/PeptideAnnotationS1_5ppm_002_recal.scanNum.rep1.rank10.pXg.fdr.dist", header = T, sep = "\t", as.is = as.double())
allResult <- read.csv(file = "PEAKS/PeptideAnnotationS1_5ppm_002_MSFrecal_addScanNum.tsv", header = T, sep = "\t", as.is = as.double())

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
  theme(text = element_text(size=25)) +
  ggtitle("All PSMs") +
  labs(y="PSMs", x = "ALC Score")

ggsave("score.all.S1.png", plot = g, width = 10, height = 8, units = "in", dpi = 300)

g <- ggplot(data = c_scoreTable, aes(x=Score, y=Count, fill=Class)) +
  scale_fill_manual(values=c(blueC, redC)) +
  theme_bw() +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme(text = element_text(size=25)) +
  ggtitle("Canonical PSMs") +
  labs(y="PSMs", x = "ALC Score")

ggsave("score.canonical.S1.png", plot = g, width = 10, height = 8, units = "in", dpi = 300)

## Noncanonical
## FDR threshold
fdr <- 0.05
fdrScore <- 0
targetCount <- 0
decoyCount <- 0

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
    if(cutFDR < 0.05) {
      fdrScore <- alc
    }
  }
}

g <- ggplot(data = nc_scoreTable, aes(x=Score, y=Count, fill=Class)) +
  scale_fill_manual(values=c(blueC,redC)) +
  theme_bw() +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(text = element_text(size=25)) +
  ggtitle("Noncanonical PSMs") +
  labs(y="PSMs", x = "ALC Score")
ggsave("score.noncanonical.S1.png", plot = g, width = 10, height = 8, units = "in", dpi = 300)

print(fdrScore)

## RNA-read distrubition
data <- read.csv(file = "pXg/PeptideAnnotationS4_5ppm_002_recal.scanNum.rep1.rank10.pXg.pval.dist", header = T, sep="\t", as.is = as.double())

data$"Peptide length" <- as.character(data$PeptideLength)
data$Mock <- log2(data$Mock+1)
data$Experiment <- log2(data$Experiment+1)
data$`Peptide length` <- factor(data$PeptideLength, levels = c('8','9','10','11','12','13', '14', '15'))
data <- data[data$ReadCount <= 20, ]


g <- ggplot(data = data, aes(x=ReadCount, y=Mock, fill=`Peptide length`)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_line(size = 1, aes(color=`Peptide length`)) +
  xlab("Mock read") +
  ylab("Log2(number of peptides + 1)") +
  scale_x_continuous(breaks=seq(from=0, to=20, by = 1)) +
  theme(text = element_text(size=25)) +
  ggtitle("Null distribution of peptide-reads mapping") +
  geom_point()

ggsave("read.S4.png", plot = g, width = 10, height = 8, units = "in", dpi = 300)

## PPM-tolerance
ppmData <-  read_excel(path = "Laumont_Results_Without_Deami.xlsx", sheet = "B-LCL4")
allData <- read.csv(file = "PEAKS/PeptideAnnotationS4_5ppm_002_MSFrecal_addScanNum.tsv", header = T, sep = "\t", as.is = as.double())

uniqueAllData <- allData[!duplicated(paste(allData$Fraction, allData$Source.File, allData$Scan, sep = "_")), ]


subjectNames <- c("B-LCL subject1", "B-LCL subject2", "B-LCL subject3", "B-LCL subject4")
subjectFDRScore <- c(73, 77, 79, 90)

allPPM <- uniqueAllData$ppm
allPPM <- data.frame(allPPM)

cPPM <- ppmData[ppmData$IsCanonical == T, ]
cPPM <- cPPM$ppm
ncPPM <- ppmData[ppmData$IsCanonical == F, ]$ppm
ncFDRPPM <- ppmData[ppmData$IsCanonical == F, ]
ncFDRPPM <- ncFDRPPM[ncFDRPPM$`Denovo score` >= 90, ]$ppm

cPPM <- data.frame(cPPM)
ncPPM <- data.frame(ncPPM)
ncFDRPPM <- data.frame(ncFDRPPM)

colnames(allPPM) <- "PPM tolerance"
colnames(cPPM) <- "PPM tolerance"
colnames(ncPPM) <- "PPM tolerance"
colnames(ncFDRPPM) <- "PPM tolerance"

allPPM$Class <- "All PSMs"
cPPM$Class <- "cPSMs (FDR <0.05)"
ncPPM$Class <- "ncPSMs"
ncFDRPPM$Class <- "ncPSMs (FDR < 0.05)"

data <- rbind(allPPM, cPPM)
data <- rbind(data, ncPPM)
data <- rbind(data, ncFDRPPM)

data$Class <- factor(data$Class, c("All PSMs", "cPSMs (FDR <0.05)", "ncPSMs", "ncPSMs (FDR < 0.05)"))

g <- ggplot(data=data) +
  geom_density(aes(x=`PPM tolerance`, colour = Class), size = 1) +
  scale_color_manual(values=c(grayC, blueC, redC, greenC)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(from=0, to=1, by = 0.1), 
                     labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE),
                     limits=c(0,0.5)) +
  theme(text = element_text(size=25), 
        plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 30, color = "darkblue")) +
  labs(y="Density", x = "PPM error")

ggsave("ppm.S4.png", plot = g, width = 10, height = 8, units = "in", dpi = 300)


## t-test is not proper to compare two distributions... in terms of shape!
wilcox.test(data[data$Class == "All PSMs", ]$`PPM tolerance`, 
       data[data$Class == "ncPSMs", ]$`PPM tolerance`, paired = F, alternative = "two.sided")

t.test(data[data$Class == "All PSMs", ]$`PPM tolerance`, 
       data[data$Class == "ncPSMs", ]$`PPM tolerance`, paired = F, alternative = "two.sided")

t.test(data[data$Class == "All PSMs", ]$`PPM tolerance`, 
       data[data$Class == "cPSMs (FDR <0.05)", ]$`PPM tolerance`, paired = F, alternative = "two.sided")

wilcox.test(data[data$Class == "All PSMs", ]$`PPM tolerance`, 
            data[data$Class == "ncPSMs (FDR < 0.05)", ]$`PPM tolerance`, paired = F, alternative = "two.sided")

t.test(data[data$Class == "All PSMs", ]$`PPM tolerance`, 
       data[data$Class == "ncPSMs (FDR < 0.05)", ]$`PPM tolerance`, paired = F, alternative = "two.sided")

t.test(data[data$Class == "ncPSMs (FDR < 0.05)", ]$`PPM tolerance`, 
       data[data$Class == "cPSMs (FDR <0.05)", ]$`PPM tolerance`, paired = F, alternative = "two.sided")

t.test(data[data$Class == "ncPSMs", ]$`PPM tolerance`, 
       data[data$Class == "cPSMs (FDR <0.05)", ]$`PPM tolerance`, paired = F, alternative = "two.sided")
