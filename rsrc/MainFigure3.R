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

sampleC <- c(`Subject 1` = redC, `Subject 2`=blueC, `Subject 3`=greenC, `Subject 4`=purpleC)

staticThemeLeftTop <- theme(text = element_text(size=25, color = "black")
                            , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                            legend.justification = c("right", "top"),
                            legend.position= c(.18, .98),
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


setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/")

dataS1 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL1_FDR10")
dataS2 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL2_FDR10")
dataS3 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL3_FDR10")
dataS4 <- read_excel(path = "pXgResults.xlsx", sheet = "B-LCL4_FDR10")

dataS1 <- dataS1[, -c(38:43)]
dataS2 <- dataS2[, -c(38:43)]
dataS3 <- dataS3[, -c(38:43)]
dataS4 <- dataS4[, -c(38:43)]

dataS1$Subject <- "Subject 1"
dataS2$Subject <- "Subject 2"
dataS3$Subject <- "Subject 3"
dataS4$Subject <- "Subject 4"

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
subData$Type <- "WB"
subData[subData$BestScore < 2, ]$Type <- "WB"
subData[subData$BestScore < 0.5, ]$Type <- "SB"
subData$Class <- "Canonical"
subData[subData$IsCanonical == "FALSE", ]$Class <- "Noncanonical"
nrow(subData[subData$Class == "Noncanonical", ])
nrow(subData[subData$Class == "Canonical", ])

nrow(subData[subData$Subject == "Subject 1" & subData$Class == "Canonical", ])
nrow(subData[subData$Subject == "Subject 2" & subData$Class == "Canonical", ])
nrow(subData[subData$Subject == "Subject 3" & subData$Class == "Canonical", ])
nrow(subData[subData$Subject == "Subject 4" & subData$Class == "Canonical", ])

nrow(subData[subData$Subject == "Subject 1" & subData$Class == "Noncanonical", ])
nrow(subData[subData$Subject == "Subject 2" & subData$Class == "Noncanonical", ])
nrow(subData[subData$Subject == "Subject 3" & subData$Class == "Noncanonical", ])
nrow(subData[subData$Subject == "Subject 4" & subData$Class == "Noncanonical", ])


subData$Subject <- factor(subData$Subject, levels = c("Subject 1", "Subject 2", "Subject 3", "Subject 4"))

### Unique selection
subDataUnique <- subData
subDataUnique <- subDataUnique[subDataUnique$Class == "Noncanonical", ]

subDataUniqueSelection <- subDataUnique[subDataUnique$EventCount == 1, ]
subDataUniqueSelection <- subDataUniqueSelection[subDataUniqueSelection$GeneIDCount <= 1, ]
subDataUniqueSelection <- subDataUniqueSelection[subDataUniqueSelection$GenomicLociCount <= 1, ]
### shared
subDataSharedSelection <- subDataUnique[subDataUnique$EventCount > 1 | subDataUnique$GeneIDCount > 1 | subDataUnique$GenomicLociCount > 1, ]

nrow(subDataUnique)
nrow(subDataUniqueSelection)
nrow(subDataSharedSelection)

nrow(subDataUniqueSelection[subDataUniqueSelection$Subject == "Subject 1", ])
nrow(subDataUniqueSelection[subDataUniqueSelection$Subject == "Subject 2", ])
nrow(subDataUniqueSelection[subDataUniqueSelection$Subject == "Subject 3", ])
nrow(subDataUniqueSelection[subDataUniqueSelection$Subject == "Subject 4", ])

tmp <- subDataUniqueSelection[subDataUniqueSelection$Subject == "Subject 1" | subDataUniqueSelection$Subject == "Subject 2" | subDataUniqueSelection$Subject == "Subject 3", ]
a1 <- nrow(tmp)
tmp <- tmp[!duplicated(tmp[,c('InferredPeptide')]), ]
a2 <- nrow(tmp)
a3 <- a1 - a2
print(a2)
print(a3)

byEventCount <- subDataSharedSelection[subDataSharedSelection$EventCount > 1, ]
byGeneIDCount <- subDataSharedSelection[subDataSharedSelection$GeneIDCount > 1, ]
byGenomicLociCount <- subDataSharedSelection[subDataSharedSelection$GenomicLociCount > 1, ]

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

t.test(subDataUnique[subDataUnique$Log10Abundance !=0 & subDataUnique$Subject == "Subject 1" & subDataUnique$Class == "Canonical", ]$Log2Abundance, subDataUnique[subDataUnique$Log10Abundance !=0 & subDataUnique$Subject == "Subject 1" & subDataUnique$Class == "Noncanonical", ]$Log2Abundance)
abPlot <- ggplot(data=subDataUnique[subDataUnique$Log10Abundance != 0 & subDataUnique$Class == "Noncanonical", ], aes(x=Events, y=Log10Abundance)) +
  theme_bw() +
  scale_fill_brewer(palette="Set3") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  staticThemeNone +
  labs(y= TeX("$Log_{10}$(Abundance+1)"), x =  NULL) +
  coord_flip()
  #facet_grid(rows = vars(Events))
abPlot
ggsave("quant.event.png", plot = abPlot, width = 12, height = 10, units = "in", dpi = 300)
#nrow(subDataUnique[subDataUnique$Class == "Noncanonical" & subDataUnique$BestScore < 2, ])
#tmp <- subDataUnique[!duplicated(subDataUnique[,c('InferredPeptide')]), ]
#nrow(tmp[tmp$Class == "Noncanonical"&tmp$BestScore < 2, ])

## Rename Events
## With/without mutations separately.


subDataUniqueSelection[subDataUniqueSelection$Events == "5'-UTR;sense", ]$Events <- "5'-UTR"
subDataUniqueSelection[subDataUniqueSelection$Events == "3'-UTR;sense", ]$Events <- "3'-UTR"
subDataUniqueSelection[subDataUniqueSelection$Events == "noncoding;sense", ]$Events <- "ncRNA"
subDataUniqueSelection[subDataUniqueSelection$Events == "frameshift;sense", ]$Events <- "frameshift"
subDataUniqueSelection[subDataUniqueSelection$Events == "intronic;sense", ]$Events <- "intron"
subDataUniqueSelection[subDataUniqueSelection$Events == "intergenic;sense", ]$Events <- "intergenic"
subDataUniqueSelection[subDataUniqueSelection$Events == "proteincoding;antisense", ]$Events <- "PC;asRNA"
subDataUniqueSelection[subDataUniqueSelection$Events == "proteincoding;sense;alternativesplicing", ]$Events <- "PC;AS"
subDataUniqueSelection[subDataUniqueSelection$Events == "noncoding;sense;alternativesplicing", ]$Events <- "ncRNA;AS"
subDataUniqueSelection[subDataUniqueSelection$Events == "noncoding;antisense;alternativesplicing", ]$Events <- "ncRNA;asRNA;AS"
subDataUniqueSelection[subDataUniqueSelection$Events == "5'-UTR;antisense", ]$Events <- "5'-UTR;asRNA"
subDataUniqueSelection[subDataUniqueSelection$Events == "3'-UTR;antisense", ]$Events <- "3'-UTR;asRNA"
subDataUniqueSelection[subDataUniqueSelection$Events == "intronic;antisense", ]$Events <- "intron;asRNA"
subDataUniqueSelection[subDataUniqueSelection$Events == "noncoding;antisense", ]$Events <- "ncRNA;asRNA"
subDataUniqueSelection[subDataUniqueSelection$Events == "unknown", ]$Events <- "unknown"
subDataUniqueSelection[subDataUniqueSelection$Events == "proteincoding;sense", ]$Events <- "PC"

#subDataUniqueSelection[subDataUniqueSelection$Mutations != "-", ]$Events <- paste(subDataUniqueSelection[subDataUniqueSelection$Mutations != "-", ]$Events, "mutations", sep=";")
subDataUniqueSelectionWOMut <- subDataUniqueSelection[subDataUniqueSelection$Mutations != "-", ]
tmp <- subDataUniqueSelection
tmp <- tmp[tmp$Subject == 'Subject 1' | tmp$Subject == 'Subject 2', ]

tmp <- tmp[!duplicated(tmp[,c('InferredPeptide')]), ]
nrow(tmp[tmp$Events == "PC", ])
tmp$GeneNames
tmp <- tmp[!duplicated(tmp[, c("GeneNames")]), ]
nrow(tmp)
subDataUniqueSelectionWOMut$Events <- factor(subDataUniqueSelectionWOMut$Events, 
                               levels = c("PC", "5'-UTR", "ncRNA", "frameshift", "3'-UTR", "intron", "intergenic", 
                                          "PC;AS", "3'-UTR;asRNA", "PC;asRNA", "ncRNA;AS",  
                                          "intron;asRNA", "ncRNA;asRNA","unknown" , 
                                          "ncRNA;asRNA;AS", "5'-UTR;asRNA"))

## mutated
topRankedPlot <- ggplot(data=subDataUniqueSelectionWOMut[subDataUniqueSelectionWOMut$Class == "Noncanonical", ], aes(x=Events, fill=Subject)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  geom_bar(aes(x= (..count..), y=reorder(Events, desc(Events))), position=position_dodge2(preserve = "single", padding = 0, reverse = T), width = 0.5) +
  staticThemeRightBottom +
  labs(y="Event", x = "A number of unique ncMAPs")
  #coord_flip() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

topRankedPlot
ggsave("ncCategoryMut.png", plot = topRankedPlot, width = 12, height = 8, units = "in", dpi = 300)

