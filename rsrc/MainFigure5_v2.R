
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

staticThemeRightTop <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), 
                             axis.text.y = element_text(size=20, color="black"),
                             legend.position= c(0.8, 0.1),
                             legend.text = element_text(size=15, color = "black"),
                             legend.title = element_text(size=0, color = "black"))

g <- ggplot(data = tmp, aes(x=Mascot_Score, y=`pXg_ALC (%)`, color = Category)) +
  theme_bw() +
  staticThemeRightTop +
  #staticThemeNone +
  geom_point() +
  scale_fill_brewer(palette = "Set1") +
  ylab("ALC score") +
  xlab("Mascot score")

g <- ggExtra::ggMarginal(g, type = "boxplot", groupColour = T, groupFill = T)
g

mean(tmp$`pXg_ALC (%)`)
mean(tmp$Mascot_Score)

ggsave("Comparison_score.png", plot = g, width = 8, height = 8, units = "in", dpi = 300)

classFlow <- tmp
classFlow[classFlow$pXg_IsCanonical == "TRUE", ]$pXg_Events <- "Coding"
classFlow[classFlow$pXg_Events == "PC;asRNA", ]$pXg_Events <- "asRNA"
classFlow[classFlow$pXg_Events == "5`-UTR|FS", ]$pXg_Events <- "FS|5`-UTR"
classFlow$Laumont_GenomicLoci <- paste("chr", classFlow$Laumont_Chromosome, ":", classFlow$Laumont_Start, "-", classFlow$Laumont_Stop, sep = "")

a <- classFlow[classFlow$Category == "Different sequence" & classFlow$Laumont_Event == classFlow$pXg_Events, 
               c("Mascot_Scan", "Laumont_Event", "pXg_Events","Mascot_Peptide", "Laumont_Peptide", "pXg_InferredPeptide", 
                 "Laumont_Ensembl_gene_id", "pXg_GeneIDs", "pXg_IsCanonical", "Laumont_GenomicLoci", "pXg_GenomicLoci")]

nrow(classFlow[classFlow$Category == "Equivalent sequence", ])
nrow(classFlow[classFlow$Category != "Equivalent sequence", ])
classFlow$Mascot_Scan
classFlow$pXg_GenomicLoci
classFlow$Mascot_Peptide

#networkLink <- classFlow[classFlow$Category == "Different sequence" & classFlow$Laumont_Event != classFlow$pXg_Events, ] %>% 
networkLink <- classFlow[classFlow$Category == "Different sequence", ] %>% 
  mutate(Laumont_Event = paste(Laumont_Event, " ")) %>% 
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

saveNetwork(networkPlot, "network_diffSeq.html")

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



