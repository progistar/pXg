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


setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXg/OriginData/")


## RNA Read Distribution
dataS1 <- read.csv(file = "S1.RAW.PEAKS.pxg.pval.1000.dist", header = T, sep="\t", as.is = as.double())
dataS2 <- read.csv(file = "S2.RAW.PEAKS.pxg.pval.1000.dist", header = T, sep="\t", as.is = as.double())
dataS3 <- read.csv(file = "S3.RAW.PEAKS.pxg.pval.1000.dist", header = T, sep="\t", as.is = as.double())
dataS4 <- read.csv(file = "S4.RAW.PEAKS.pxg.pval.1000.dist", header = T, sep="\t", as.is = as.double())

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

data[data$PeptideLength == 8 & data$ReadCount == 0, ]

tLevelData <- data.frame(Count = data$Experiment, Length=data$PeptideLength, ReadCount = data$ReadCount, Label = 'Target', Subject = data$Subject)
dLevelData <- data.frame(Count = data$Mock, Length=data$PeptideLength , ReadCount = data$ReadCount, Label = 'Decoy', Subject = data$Subject)
tdLevelData <- rbind(tLevelData, dLevelData)

sampleC <- c(`Subject 1` = redC, `Subject 2`=blueC, `Subject 3`=greenC, `Subject 4`=purpleC)
LabelC <- c(`Target`=blueC, `Decoy` = redC)

g <- ggplot(data = tdLevelData[(tdLevelData$Length <= 11  ) & tdLevelData$ReadCount <= 1000, ], aes(x=ReadCount, y=Count, fill=`Label`)) +
  theme_bw() +
  scale_color_manual(values=LabelC) +
  geom_line(size = 1, aes(color = Label), alpha = 0.8) +
  #geom_density(size = 1, aes(color = Label, y=Count), alpha = 0.8) +
  xlab("Read") +
  ylab(TeX("$Log_{2}$(number of peptides + 1)")) +
  #scale_x_continuous(breaks=seq(from=0, to=40, by = 1)) +
  staticThemeRightTop +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in")) +
  #geom_point() +
  facet_grid(rows = vars(`Subject`), cols = vars(`Length`))
g
ggsave("read.length8_11.1000.png", plot = g, width = 28, height = 9, units = "in", dpi = 300)


  
