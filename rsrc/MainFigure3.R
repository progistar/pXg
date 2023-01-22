
install.packages("networkD3")
webshot::install_phantomjs()
install.packages("ggExtra")

library(ggExtra)
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
library(networkD3)
library(webshot)


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

staticThemeRightTop <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                             legend.justification = c("right", "top"),
                             legend.position= c(.98, .98),
                             legend.text = element_text(size=20, color = "black"),
                             legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeLeftTop <- theme(text = element_text(size=25, color = "black")
                            , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                            legend.justification = c("right", "top"),
                            legend.position= c(.22, .98),
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
                          legend.text = element_text(size=25, color = "black"))


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

tmp <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/Laumont_Comparison/Comparison_v2.xlsx", sheet = "ncMAP", skip = 1)

## Select pXg Result
pXgRes <- tmp[, c("Events", "FastaIDs", "EventCount", "GeneNameCount", "GenomicLociCount", "MHC-I")]
pXgRes <- pXgRes[pXgRes$FastaIDs == "-" & pXgRes$GeneNameCount < 2 & pXgRes$GenomicLociCount < 2 & pXgRes$EventCount < 2, ]
pXgRes <- pXgRes[pXgRes$Events != "-" & pXgRes$`MHC-I` != "NB", ]

## Select Luamont Result
lauRes <- tmp[, c("Events_L")]
lauRes <- lauRes[lauRes$Events_L != "-", ]

res <- lauRes
colnames(res) <- c("Event")
res$Category <- "Laumont"

pXgRes <- pXgRes[, c("Events")]
colnames(pXgRes) <- c("Event")
pXgRes$Category <- "pXg"

res <- rbind(res, pXgRes)

#IEDBData$Event <- factor(IEDBData$Event, levels = c(ncRNA, n5UTR, FS, IR, Coding, n3UTR, IGR, Unknown))
#IEDBData$Category <- factor(IEDBData$Category, levels = c("Fully annotated", "Partially annotated", "NA"))

res[res$Event == "3`-UTR;asRNA" | res$Event == "5`-UTR;asRNA" | 
      res$Event == "PC;asRNA", ]$Event <- "asRNA"

res[res$Event == "3`-UTR;sense", ]$Event <- "3`-UTR"
res[res$Event == "5`-UTR;sense", ]$Event <- "5`-UTR"
res[res$Event == "FS;sense", ]$Event <- "FS"
res[res$Event == "IGR;sense", ]$Event <- "IGR"
res[res$Event == "IR;sense", ]$Event <- "IR"
res[res$Event == "IR;sense;AS", ]$Event <- "IR;AS"
res[res$Event == "ncRNA;sense", ]$Event <- "ncRNA"
res[res$Event == "PC;sense" | res$Event == "PC;sense;Mutant", ]$Event <- "Coding (mutated)"
res[res$Event == "PC;sense;AS", ]$Event <- "Coding;AS"
res[res$Event == "unknown", ]$Event <- "Unknown"

g <- ggplot(data = res, aes(x=Event, fill = Category)) +
  theme_bw() +
  staticThemeRightTop +
  geom_bar(position = position_dodge2(preserve = "single")) +
  scale_fill_brewer(palette = "Set1") +
  ylab("ncMAP") +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
g

ggsave("EventComparison.png", plot = g, width = 10, height = 6, units = "in", dpi = 300)

targetEvent <- "Coding (mutated)"
allpXg <- nrow(res[res$Category == "pXg", ])
allLaumont <- nrow(res[res$Category == "Laumont", ])
targetEventNpXg <- nrow(res[res$Event == targetEvent & res$Category == "pXg", ])
targetEventNpLaumont <- nrow(res[res$Event == targetEvent & res$Category == "Laumont", ])

targetEventNpLaumont
targetEventNpXg

## Fisher exact test
x <- matrix(data = c(targetEventNpXg, targetEventNpLaumont, 
                     allpXg-targetEventNpXg, allLaumont-targetEventNpLaumont), 
            nrow = 2, ncol = 2, byrow = FALSE)

fisher.test(x)




