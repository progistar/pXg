install.packages("ggpmisc")

library(ggplot2)
library(RColorBrewer)
library("dplyr")
library(readxl)
library(latex2exp)
library("stringr")
library(ggpmisc)
library(ComplexHeatmap)
library(ggpubr)

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
                          legend.justification = c("right"),
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


setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/")


scaData <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/ProteomeToolsHLA/UniqueNoncanonicalPSMs.xlsx", 
                      sheet = "SA")
scaData$Class <- "B-LCL"
scaData$SA <- scaData$`SA B-LCL`

scaDataPT <- scaData
scaDataPT$Class <- "ProteomeTools"
scaDataPT$SA <- scaDataPT$`SA ProteomeTools`
scaData <- rbind(scaData, scaDataPT)
scaData$mStatus <- 1
scaData[scaData$Mutations == "-", ]$mStatus <- 0
scaData[scaData$Mutations == "chr22:22792800T>C|chr22:22792803A>C", ]$mStatus <- 2
## Coding 14
## 5`-UTR 20
## FS 13
## ncRNA 4
## 3`-UTR 4
## IR 2
## IGR 1
## IR;AS 1
scaData$Events <- factor(scaData$Events, levels = c("5`-UTR", "Coding", "FS", "ncRNA", "3`-UTR", "IR", "IR;AS", "IGR"))
scaDataBLCL <- scaData[scaData$Class == "B-LCL", ]
scaDataBLCL$Event <- scaDataBLCL$Events
scaDataBLCL$Precursor <- factor(scaDataBLCL$Precursor,
                                levels = scaDataBLCL$Precursor[order(scaDataBLCL$Events, decreasing = FALSE)])


boxPlot_ <- ggplot(data = scaData, aes(x=Class, y=SA, fill=`Class`)) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_manual(values = c(blueC, redC)) +
  scale_y_continuous(breaks = seq(from = 0, to =1, by = 0.1), limits = c(0, 1)) +
  xlab(NULL) +
  ylab("SA") +
  staticThemeNone +
  theme(legend.key.size = unit(0.2, "in"), legend.key.width = unit(0.6, "in"))
boxPlot_

ggsave("SCA_box.png", plot = boxPlot_, width = 5, height = 8, units = "in", dpi = 300)

mutPlot <- ggplot(data=scaDataBLCL, aes(x=Precursor, y=mStatus, fill = "Count")) +
  scale_fill_grey() +
  theme_bw() +
  scale_y_continuous(breaks = seq(from = 0, to =2, by = 1), limits = c(0, 2)) +
  geom_bar(stat = "identity", position = position_dodge())+
  staticThemeRight +
  theme(text = element_text(size=20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(title = "Mutation")) +
  labs(y= NULL, x = NULL)

mutPlot

ggsave("SCA_mut.png", plot = mutPlot, width = 19, height = 1, units = "in", dpi = 300)

scaPlot <- ggplot(data=scaDataBLCL, aes(x=Precursor, y=SA, fill=Event)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_bar(stat = "identity", position = position_dodge())+
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  staticThemeRight +
  labs(y= "SA", x = "Precursor")

scaPlot
ggsave("SCA_Hist.png", plot = scaPlot, width = 20, height = 6, units = "in", dpi = 300)

scatterPlot <- ggplot(data = scaData, aes(x=`SA B-LCL`, y=`SA ProteomeTools`, fill = Events)) +
  geom_point()

scatterPlot



