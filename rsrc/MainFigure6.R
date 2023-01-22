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


IEDBData <- read_excel(path = "pXgResults_p001_FDR.xlsx", sheet = "ncMAP")
IEDBData$Event <- IEDBData$Events

ncRNA <- paste("ncRNA\n(",nrow(IEDBData[IEDBData$Events == "ncRNA",] ), ")", sep = "")
n5UTR <- paste("5`-UTR\n(",nrow(IEDBData[IEDBData$Events == "5`-UTR",] ), ")", sep = "")
FS <- paste("FS\n(",nrow(IEDBData[IEDBData$Events == "FS",] ), ")", sep = "")
IR <- paste("IR\n(",nrow(IEDBData[IEDBData$Events == "IR",] ), ")", sep = "")
Coding <- paste("Coding\n(mutated)\n(",nrow(IEDBData[IEDBData$Events == "Coding",] ), ")", sep = "")
n3UTR <- paste("3`-UTR\n(",nrow(IEDBData[IEDBData$Events == "3`-UTR",] ), ")", sep = "")
IGR <- paste("IGR\n(",nrow(IEDBData[IEDBData$Events == "IGR",] ), ")", sep = "")
Unknown <- paste("Unknown\n(",nrow(IEDBData[IEDBData$Events == "Unknown",] ), ")", sep = "")


IEDBData[IEDBData$Event == "ncRNA",]$Event <- ncRNA
IEDBData[IEDBData$Event == "5`-UTR",]$Event <- n5UTR
IEDBData[IEDBData$Event == "FS",]$Event <- FS
IEDBData[IEDBData$Event == "IR",]$Event <- IR
IEDBData[IEDBData$Event == "Coding",]$Event <- Coding
IEDBData[IEDBData$Event == "3`-UTR",]$Event <- n3UTR
IEDBData[IEDBData$Event == "IGR",]$Event <- IGR
IEDBData[IEDBData$Event == "Unknown",]$Event <- Unknown
  
  
  

IEDBData$Event <- factor(IEDBData$Event, levels = c(ncRNA, n5UTR, FS, IR, Coding, n3UTR, IGR, Unknown))
IEDBData$Category <- factor(IEDBData$Category, levels = c("Fully annotated", "Partially annotated", "NA"))
g <- ggplot(data = IEDBData, aes(x=Event, fill = Category)) +
  theme_bw() +
  staticThemeTop +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "Set1") +
  ylab("MAP") +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
g

ggsave("IEDB.png", plot = g, width = 12, height = 6, units = "in", dpi = 300)



