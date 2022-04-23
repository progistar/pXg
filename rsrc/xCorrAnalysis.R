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
dataS1 <- read_excel(path = "MSGF_Summary.xlsx", sheet = "B-LCL1_FDR5")
dataS2 <- read_excel(path = "MSGF_Summary.xlsx", sheet = "B-LCL2_FDR5")
dataS3 <- read_excel(path = "MSGF_Summary.xlsx", sheet = "B-LCL3_FDR5")
dataS4 <- read_excel(path = "MSGF_Summary.xlsx", sheet = "B-LCL4_FDR5")

dataS1 <- dataS1[, -c(20:25)]
dataS2 <- dataS2[, -c(20:25)]
dataS3 <- dataS3[, -c(20:25)]
dataS4 <- dataS4[, -c(20:25)]

MSGFData <- rbind(dataS1, dataS2, dataS3, dataS4)
MSGFData$Class <- "cPSM (MS-GF+)"
xCorrMSGF <- MSGFData[, c("xCorr", "Class")]

setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/")
dataS1 <- read_excel(path = "pXg_Summary.xlsx", sheet = "B-LCL1_FDR5")
dataS2 <- read_excel(path = "pXg_Summary.xlsx", sheet = "B-LCL2_FDR5")
dataS3 <- read_excel(path = "pXg_Summary.xlsx", sheet = "B-LCL3_FDR5")
dataS4 <- read_excel(path = "pXg_Summary.xlsx", sheet = "B-LCL4_FDR5")

dataS1 <- dataS1[, -c(38:43)]
dataS2 <- dataS2[, -c(38:43)]
dataS3 <- dataS3[, -c(38:43)]
dataS4 <- dataS4[, -c(38:43)]

pXgData <- rbind(dataS1, dataS2, dataS3, dataS4)
nrow(pXgData[pXgData$IsCanonical == "TRUE", ])
nrow(pXgData[pXgData$IsCanonical == "FALSE", ])

mean(pXgData$xCorr)
mean(pXgData[pXgData$IsCanonical == "TRUE", ]$xCorr)

pXgData$Class <- "cPSM (pXg)"
pXgData[pXgData$IsCanonical == "FALSE", ]$Class <- "ncPSM (pXg)"

xCorrpXg <- pXgData[, c("xCorr","Class")]

## Combine
xCorr <- rbind(xCorrMSGF, xCorrpXg)
xCorr$Class <- factor(xCorr$Class, levels = c("cPSM (MS-GF+)", "cPSM (pXg)", "ncPSM (pXg)"))

## median for xCorr among three cases
medCanonicalpXg <- median(xCorr[xCorr$Class == "cPSM (pXg)", ]$xCorr)
medCanonicalMSGF <- median(xCorr[xCorr$Class == "cPSM (MS-GF+)", ]$xCorr)
medNoncanonicalpXg <- median(xCorr[xCorr$Class == "ncPSM (pXg)", ]$xCorr)

## Histogram 
g <- ggplot(data=xCorr, aes(x=`xCorr`, colour = Class)) +
  theme_bw() +
  theme(text = element_text(size=25)) +
  geom_histogram(alpha=0, position="identity", bins = 50, fill="transparent", lwd = 1) +
  scale_colour_manual(values = c(greenC, blueC, redC)) +
  labs(y="PSMs", x = "xCorr score") +
  geom_segment(aes(x=medCanonicalpXg, xend=medCanonicalpXg, y=0, yend=9000), color = blueC, size = 1, linetype = "dotted") +
  geom_segment(aes(x=medCanonicalMSGF, xend=medCanonicalMSGF, y=0, yend=8500), color = greenC, size = 1, linetype = "dotted") +
  geom_segment(aes(x=medNoncanonicalpXg, xend=medNoncanonicalpXg, y=0, yend=8000), color = redC, size = 1, linetype = "dotted") +
  annotate("text", label = 1.88, x = medCanonicalpXg+0.4, y = 9000, size = 7, family="serif", colour = blueC) +
  annotate("text", label = 1.96, x = medCanonicalMSGF+0.4, y = 8500, size = 7, family="serif", colour = greenC) +
  annotate("text", label = 2.01, x = medNoncanonicalpXg+0.4, y = 8000, size = 7, family="serif", colour = redC) +
  ggtitle("xCorr score distributions")
g


ggsave("pXg.MSGF.xCorr.density.png", plot = g, width = 10, height = 8, units = "in", dpi = 300)

g <- ggplot(data=xCorr, aes(x=`xCorr`, colour = Class)) +
  theme_bw() +
  theme(text = element_text(size=25)) +
  geom_density(alpha=0.3, position="identity", lwd = 1) +
  scale_colour_manual(values = c(greenC, blueC, redC)) +
  labs(y="Density", x = "xCorr score") +
  geom_segment(aes(x=medCanonicalpXg, xend=medCanonicalpXg, y=0, yend=1), color = blueC, size = 1, linetype = "dotted") +
  geom_segment(aes(x=medCanonicalMSGF, xend=medCanonicalMSGF, y=0, yend=0.95), color = greenC, size = 1, linetype = "dotted") +
  geom_segment(aes(x=medNoncanonicalpXg, xend=medNoncanonicalpXg, y=0, yend=0.9), color = redC, size = 1, linetype = "dotted") +
  annotate("text", label = 1.88, x = medCanonicalpXg+0.4, y = 1, size = 7, family="serif", colour = blueC) +
  annotate("text", label = 1.96, x = medCanonicalMSGF+0.4, y = 0.95, size = 7, family="serif", colour = greenC) +
  annotate("text", label = 2.01, x = medNoncanonicalpXg+0.4, y = 0.9, size = 7, family="serif", colour = redC) +
  ggtitle("xCorr score distributions")
g




t.test(xCorr[xCorr$Class == "Canonical (MS-GF+)", ]$xCorr, xCorr[xCorr$Class == "Noncanonical (pXg)", ]$xCorr)
