
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

tmp <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/Laumont_Comparison/Comparison_v2.xlsx", sheet = "Mascot-reanalysis-summary")

tmp$ALCScore <- 0
tmp[tmp$Category == "Overlap", ]$ALCScore <- tmp[tmp$Category == "Overlap", ]$`pXg_ALC (%)`
tmp[tmp$Category == "Laumont", ]$ALCScore <- tmp[tmp$Category == "Laumont", ]$`PEAKS_ALC (%)`
tmp[tmp$Category == "pXg", ]$ALCScore <- tmp[tmp$Category == "pXg", ]$`pXg_ALC (%)`

staticThemeRightTop <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), 
                             axis.text.y = element_text(size=20, color="black"),
                             legend.position= c(0.85, 0.1),
                             legend.text = element_text(size=15, color = "black"),
                             legend.title = element_text(size=0, color = "black"))

tmp$Category <- factor(tmp$Category, levels = c("pXg", "Overlap", "Laumont"))
tmp$Comparison <- "Overlap"
tmp[tmp$Category == "pXg", ]$Comparison <- "pXg only"
tmp[tmp$Category == "Laumont", ]$Comparison <- "Laumont only"
tmp$Comparison <- factor(tmp$Comparison, levels = c("pXg only", "Overlap", "Laumont only"))

nrow(tmp[tmp$Laumont_Event == "asRNA", ])
nrow(tmp[tmp$pXg_Events == "asRNA", ])
g <- ggplot(data = tmp, aes(x=MascotScore, y=ALCScore, colour=Comparison)) +
  theme_bw() +
  staticThemeRightTop +
  #staticThemeNone +
  geom_point() +
  scale_fill_brewer(palette = "Set1") +
  ylab("ALC score") +
  xlab("Mascot score")

g <- ggExtra::ggMarginal(g, type = "boxplot", groupFill = TRUE, groupColour = T)
g

mean(tmp[tmp$Category == "pXg", ]$`pXg_ALC (%)`)
mean(tmp[tmp$Category == "pXg", ]$MascotScore)
mean(tmp[tmp$Category == "Laumont", ]$ALCScore)
mean(tmp[tmp$Category == "Laumont", ]$MascotScore)

mean(tmp[tmp$Category == "Overlap", ]$`pXg_ALC (%)`)
mean(tmp[tmp$Category == "Overlap", ]$MascotScore)



ggsave("Comparison_score.png", plot = g, width = 8, height = 8, units = "in", dpi = 300)


tmp$Class <- "NA"

## Laumont => pXg
tmpLaumontTopXg <- tmp[tmp$Category == "Laumont" & tmp$pXg_Fraction != 0, ]
tmpLaumontTopXg$Class <- "Unique"
tmpLaumontTopXg[tmpLaumontTopXg$`pXg_MHC-I` == "NB", ]$Class <- "NB"
tmpLaumontTopXg[tmpLaumontTopXg$`pXg_MHC-I` != "NB" & 
                  (tmpLaumontTopXg$pXg_GeneNameCount > 1 |
                     tmpLaumontTopXg$pXg_EventCount > 1 |
                     tmpLaumontTopXg$pXg_GenomicLociCount > 1
                     ), ]$Class <- "Shared"
tmpLaumontTopXg[tmpLaumontTopXg$pXg_Match == "Equivalent AA composition", ]$pXg_Match <- "Equivalent"
tmpLaumontTopXg[tmpLaumontTopXg$pXg_Match == "Different AA composition", ]$pXg_Match <- "Different"

tmppXgToLaumont <- tmp[tmp$Category == "pXg" & tmp$Laumont_Peptide != 0, ]
tmppXgToLaumont$Class <- "Unique"
tmppXgToLaumont[tmppXgToLaumont$pXg_Match == "Equivalent AA composition", ]$pXg_Match <- "Equivalent"
tmppXgToLaumont[tmppXgToLaumont$pXg_Match == "Different AA composition", ]$pXg_Match <- "Different"

tmpChange <- rbind(tmppXgToLaumont, tmpLaumontTopXg)

tmpChange$pXg_Match <- factor(tmpChange$pXg_Match, c("Different", "Equivalent"))
tmpChange$Comparison <- "pXg only"
tmpChange[tmpChange$Category == "Laumont", ]$Comparison <- "Laumont only"

a <- tmpChange[tmpChange$Class == "Unique",]
write.table(a, file = "73unique_changes.tsv", sep = "\t")


nrow(tmpChange[tmpChange$Class == "Unique" & tmpChange$pXg_Events == "Coding", ])
nrow(tmpChange[tmpChange$Class == "Unique" & tmpChange$Laumont_Event == "Coding", ])

g <- ggplot(data = tmpChange, aes(x=Class, fill = Class)) +
  theme_bw() +
  #staticThemeRight +
  staticThemeNone +
  geom_bar() +
  scale_fill_brewer(palette = "Set1") +
  ylab("MAP") +
  xlab("") +
  #scale_x_discrete(position = "top") +
  #rotate() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
  #facet_grid(cols = vars(`Comparison`), switch = "y", scales = "free")
g

ggsave("LtoP.png", plot = g, width = 6, height = 10, units = "in", dpi = 300)

nrow(tmpChange[tmpChange$pXg_Events != "Coding" & tmpChange$Laumont_Event == "Coding", ])

#lOnly[lOnly$pXg_Events == "Coding;asRNA", ]$pXg_Events <- "asRNA"
tmpChange[tmpChange$pXg_Events == "PC;asRNA", ]$pXg_Events <- "asRNA"
tmpChange[tmpChange$pXg_IsCanonical == "TRUE", ]$pXg_Events <- "Coding"

networkLink <- tmpChange[tmpChange$Class == "Unique", ] %>% mutate(Laumont_Event = paste(Laumont_Event, " ")) %>% 
  group_by(source = Laumont_Event, target = pXg_Events) %>% 
  summarise(value = n()) %>%
  as.data.frame()

networkNode <- data.frame(name=c(as.character(networkLink$source), 
                                 as.character(networkLink$target)) %>% unique())

networkLink$IDsource <- match(networkLink$source, networkNode$name)-1 
networkLink$IDtarget <- match(networkLink$target, networkNode$name)-1

networkPlot <- sankeyNetwork(Links = networkLink, Nodes = networkNode,
                             Source = "IDsource", Target = "IDtarget",
                             Value = "value", NodeID = "name", sinksRight=FALSE, fontSize = 15, 
                             nodeWidth = 30, nodePadding = 15, height = 350, width = 600)

networkPlot

saveNetwork(networkPlot, "network.html")

g <- ggplot(data = tmpChange, aes(x=MascotScore, y=ALCScore, colour=Category)) +
  theme_bw() +
  staticThemeRightTop +
  #staticThemeNone +
  geom_point() +
  scale_fill_brewer(palette = "Set1") +
  ylab("ALC score") +
  xlab("Mascot score")

g <- ggExtra::ggMarginal(g, type = "boxplot", groupFill = TRUE, groupColour = T)
g


ggsave("Comparison_Diff_Equi.png", plot = g, width = 8, height = 8, units = "in", dpi = 300)



