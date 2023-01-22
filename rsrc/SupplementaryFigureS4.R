install.packages("ggpmisc")
install.packages("heatmaply")

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
library(enrichplot)

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

staticThemeRightTop <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                             legend.justification = c("right", "top"),
                             legend.text = element_text(size=20, color = "black"))

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


setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pNovo3/")

# Subject 1: 40.623265
# Subject 2: 43.07716996
# Subject 3: 53.244961
# Subject 4: 35.841886

subject1_score <- 40.623265
subject2_score <- 43.07716996
subject3_score <- 53.244961
subject4_score <- 35.841886

tmp <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 1")
tmp <- tmp[!duplicated(tmp$InferredPeptide), ]
nrow(tmp[tmp$IsCanonical == "TRUE" & tmp$`MHC-I` != "NB", ])
nrow(tmp[tmp$IsCanonical == "TRUE", ])
nrow(tmp[tmp$IsCanonical == "FALSE" & tmp$`MHC-I` != "NB" & tmp$`Score` >= subject1_score, ])
nrow(tmp[tmp$IsCanonical == "FALSE" & tmp$`Score` >= subject1_score, ])

dataS1 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 1")
dataS2 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 2")
dataS3 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 3")
dataS4 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 4")
dataS1 <- dataS1[, -c(24:29)]
dataS2 <- dataS2[, -c(24:29)]
dataS3 <- dataS3[, -c(24:29)]
dataS4 <- dataS4[, -c(24:29)]

dataS1 <- dataS1[dataS1$Score >= subject1_score | dataS1$IsCanonical == TRUE, ]
dataS2 <- dataS2[dataS2$Score >= subject2_score | dataS2$IsCanonical == TRUE, ]
dataS3 <- dataS3[dataS3$Score >= subject3_score | dataS3$IsCanonical == TRUE, ]
dataS4 <- dataS4[dataS4$Score >= subject4_score | dataS4$IsCanonical == TRUE, ]

nrow(dataS1)
nrow(dataS2)
nrow(dataS3)
nrow(dataS4)

dataS1$Sample <- "Subject 1"
dataS2$Sample <- "Subject 2"
dataS3$Sample <- "Subject 3"
dataS4$Sample <- "Subject 4"

data <- rbind(dataS1, dataS2, dataS3, dataS4)

data$Class <- "Canonical"
data[data$IsCanonical == "TRUE", ]$Class <- "Canonical"
data[data$IsCanonical == "FALSE", ]$Class <- "Noncanonical"
data$Type <- "NB"
data[data$BestScore < 2, ]$Type <- "WB"
data[data$BestScore < 0.5, ]$Type <- "SB"
data_rank <- data.frame(Rank=numeric(), Type=character(), Class=character(), Subject=character(), PSM=numeric())

## Canonical ratio
nrow(data[data$Sample == "Subject 1" & data$Type != "NB" & data$Class == "Canonical", ])/nrow(data[data$Sample == "Subject 1" & data$Class == "Canonical", ])
nrow(data[data$Sample == "Subject 2" & data$Type != "NB" & data$Class == "Canonical", ])/nrow(data[data$Sample == "Subject 2" & data$Class == "Canonical", ])
nrow(data[data$Sample == "Subject 3" & data$Type != "NB" & data$Class == "Canonical", ])/nrow(data[data$Sample == "Subject 3" & data$Class == "Canonical", ])
nrow(data[data$Sample == "Subject 4" & data$Type != "NB" & data$Class == "Canonical", ])/nrow(data[data$Sample == "Subject 4" & data$Class == "Canonical", ])

## Noncanonical ratio
nrow(data[data$Sample == "Subject 1" & data$Type != "NB" & data$Class == "Noncanonical", ])/nrow(data[data$Sample == "Subject 1" & data$Class == "Noncanonical", ])
nrow(data[data$Sample == "Subject 2" & data$Type != "NB" & data$Class == "Noncanonical", ])/nrow(data[data$Sample == "Subject 2" & data$Class == "Noncanonical", ])
nrow(data[data$Sample == "Subject 3" & data$Type != "NB" & data$Class == "Noncanonical", ])/nrow(data[data$Sample == "Subject 3" & data$Class == "Noncanonical", ])
nrow(data[data$Sample == "Subject 4" & data$Type != "NB" & data$Class == "Noncanonical", ])/nrow(data[data$Sample == "Subject 4" & data$Class == "Noncanonical", ])


for(idx in c(1: nrow(data))) {
  sub_data <- data[idx, ]
  rank <- sub_data$Rank
  type <- sub_data$Type
  class_ <- sub_data$Class
  subject <- sub_data$Sample
  
  if(nrow(data_rank[data_rank$Rank == rank & data_rank$Type == type & data_rank$Class == class_ & data_rank$Subject == subject, ]) == 0){
    data_rank[nrow(data_rank)+1, ] = c(rank, type, class_, subject, 1)
    
  } else {
    data_rank[data_rank$Rank == rank & data_rank$Type == type & data_rank$Class == class_ & data_rank$Subject == subject, ]$PSM <- as.numeric(data_rank[data_rank$Rank == rank & data_rank$Type == type & data_rank$Class == class_ & data_rank$Subject == subject, ]$PSM)+1
  }
}

data_rank$Rank <- factor(data_rank$Rank, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
data_rank$PSM <- as.numeric(data_rank$PSM)
data_rank$Type <- factor(data_rank$Type, levels = c("SB", "WB", "NB"))
data_rank$isTop <- "Rank N"
data_rank[data_rank$Rank == "1", ]$isTop <- "Rank 1"

inc <- sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Canonical"  & data_rank$Subject == "Subject 1", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 1" & data_rank$Class == "Canonical", ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Canonical"  & data_rank$Subject == "Subject 2", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 2" & data_rank$Class == "Canonical", ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Canonical"  & data_rank$Subject == "Subject 3", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 3" & data_rank$Class == "Canonical", ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Canonical"  & data_rank$Subject == "Subject 4", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 4" & data_rank$Class == "Canonical", ]$PSM)
inc/4

inc <- sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Noncanonical"  & data_rank$Subject == "Subject 1", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 1" & data_rank$Class == "Noncanonical" , ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Noncanonical"  & data_rank$Subject == "Subject 2", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 2" & data_rank$Class == "Noncanonical" , ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Noncanonical"  & data_rank$Subject == "Subject 3", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 3" & data_rank$Class == "Noncanonical" , ]$PSM)
inc <- inc + sum(data_rank[data_rank$isTop != "Rank 1" & data_rank$Type != "NB" & data_rank$Class == "Noncanonical"  & data_rank$Subject == "Subject 4", ]$PSM)/sum(data_rank[data_rank$isTop == "Rank 1" & data_rank$Type != "NB" & data_rank$Subject == "Subject 4" & data_rank$Class == "Noncanonical" , ]$PSM)
inc/4

g <- ggplot(data = data_rank[data_rank$Type != "NB", ], aes(x=Subject, y=PSM, fill=isTop)) +
  theme_bw() +
  geom_bar(stat = "identity", width = 0.8)+
  scale_fill_manual(values = c(blueC, greenC)) +
  ylab("PSM") +
  xlab(NULL) +
  staticThemeRight +
  facet_grid(rows = vars(`Class`), scales = "free") +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +
  labs(fill = "Rank")
g

ggsave("PSMs.png", plot = g, width = 7, height = 5.5, units = "in", dpi = 300)


#figure2 <- ggarrange(binderPlot, idPlot, ncol = 2, nrow = 1, common.legend = T)
#ggsave("Figure2.png", plot = figure2, width = 12, height = 6, units = "in", dpi = 300)

## Event catalog
dataS1 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 1")
dataS2 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 2")
dataS3 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 3")
dataS4 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 4")

dataS1 <- dataS1[, -c(24:29)]
dataS2 <- dataS2[, -c(24:29)]
dataS3 <- dataS3[, -c(24:29)]
dataS4 <- dataS4[, -c(24:29)]

dataS1 <- dataS1[dataS1$Score >= subject1_score | dataS1$IsCanonical == TRUE, ]
dataS2 <- dataS2[dataS2$Score >= subject2_score | dataS2$IsCanonical == TRUE, ]
dataS3 <- dataS3[dataS3$Score >= subject3_score | dataS3$IsCanonical == TRUE, ]
dataS4 <- dataS4[dataS4$Score >= subject4_score | dataS4$IsCanonical == TRUE, ]

dataS1$Subject <- "Subject 1"
dataS2$Subject <- "Subject 2"
dataS3$Subject <- "Subject 3"
dataS4$Subject <- "Subject 4"

allData <- rbind(dataS1, dataS2, dataS3, dataS4)

subData <- allData[allData$BestScore < 2, ]
subData$Type <- "WB"
subData[subData$BestScore < 2, ]$Type <- "WB"
subData[subData$BestScore < 0.5, ]$Type <- "SB"
subData$Class <- "Canonical"
subData[subData$IsCanonical == "FALSE", ]$Class <- "Noncanonical"
subData$Subject <- factor(subData$Subject, levels = c("Subject 1", "Subject 2", "Subject 3", "Subject 4"))

### Unique selection
subDataUnique <- subData[subData$Class == "Noncanonical", ]

subDataUnique <- subDataUnique[subDataUnique$EventCount == 1, ]
subDataUnique <- subDataUnique[subDataUnique$GeneIDCount <= 1, ]
subDataUnique <- subDataUnique[subDataUnique$GenomicLociCount <= 1, ]

tmp <- subDataUnique[subDataUnique$GeneIDCount > 1 | subDataUnique$EventCount > 1 | subDataUnique$GenomicLociCount > 1, ]
nrow(tmp[!duplicated(tmp$InferredPeptide), ])

tmp <- tmp[!duplicated(tmp$InferredPeptide), ]

# meta characteristics
subDataUnique$AS <- FALSE
subDataUnique$asRNA <- FALSE

subDataUnique[subDataUnique$Events == "PC;sense;AS", ]$AS <- TRUE
subDataUnique[subDataUnique$Events == "ncRNA;asRNA;AS", ]$AS <- TRUE
subDataUnique[subDataUnique$Events == "ncRNA;sense;AS", ]$AS <- TRUE
subDataUnique[subDataUnique$Events == "IR;sense;AS", ]$AS <- TRUE
subDataUnique[subDataUnique$Events == "IR;asRNA;AS", ]$AS <- TRUE

subDataUnique[subDataUnique$Events == "PC;asRNA", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "5`-UTR;asRNA", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "3`-UTR;asRNA", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "IR;asRNA", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "IR;asRNA;AS", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "ncRNA;asRNA;AS", ]$asRNA <- TRUE
subDataUnique[subDataUnique$Events == "ncRNA;asRNA", ]$asRNA <- TRUE


subDataUnique[subDataUnique$Events == "5`-UTR;sense", ]$Events <- "5`-UTR"
subDataUnique[subDataUnique$Events == "3`-UTR;sense", ]$Events <- "3`-UTR"
subDataUnique[subDataUnique$Events == "ncRNA;sense", ]$Events <- "ncRNA"
subDataUnique[subDataUnique$Events == "FS;sense", ]$Events <- "FS"
subDataUnique[subDataUnique$Events == "IR;sense", ]$Events <- "IR"
subDataUnique[subDataUnique$Events == "IR;sense;AS", ]$Events <- "IR"
subDataUnique[subDataUnique$Events == "IR;asRNA;AS", ]$Events <- "IR"
subDataUnique[subDataUnique$Events == "IGR;sense", ]$Events <- "IGR"
subDataUnique[subDataUnique$Events == "PC;asRNA", ]$Events <- "Coding"
subDataUnique[subDataUnique$Events == "PC;sense;AS", ]$Events <- "Coding"
subDataUnique[subDataUnique$Events == "ncRNA;sense;AS", ]$Events <- "ncRNA"
subDataUnique[subDataUnique$Events == "ncRNA;asRNA;AS", ]$Events <- "ncRNA"
subDataUnique[subDataUnique$Events == "5`-UTR;asRNA", ]$Events <- "5`-UTR"
subDataUnique[subDataUnique$Events == "3`-UTR;asRNA", ]$Events <- "3`-UTR"
subDataUnique[subDataUnique$Events == "IR;asRNA", ]$Events <- "IR"
subDataUnique[subDataUnique$Events == "ncRNA;asRNA", ]$Events <- "ncRNA"
subDataUnique[subDataUnique$Events == "unknown", ]$Events <- "Unknown"
subDataUnique[subDataUnique$Events == "PC;sense", ]$Events <- "Coding"

subDataUnique$Mutation <- TRUE
subDataUnique[subDataUnique$Mutations == "-", ]$Mutation <- FALSE

subDataUnique$Events <- factor(subDataUnique$Events, 
                               levels = c("Coding", "5`-UTR", "ncRNA", "FS", "3`-UTR", "IR", "IGR", 
                                          "Unknown"))

subDataUnique$InferredPeptide
subDataUnique$Reads

peptideLevel <- ddply(subDataUnique, .(InferredPeptide, Subject, asRNA, AS, Mutation, GenomicLoci, Events, Reads), nrow)
peptideLevel$`Subject 1` <- 0
peptideLevel$`Subject 2` <- 0
peptideLevel$`Subject 3` <- 0
peptideLevel$`Subject 4` <- 0

peptideLevel[peptideLevel$Subject == "Subject 1",]$`Subject 1` <- peptideLevel[peptideLevel$Subject == "Subject 1",]$V1
peptideLevel[peptideLevel$Subject == "Subject 2",]$`Subject 2` <- peptideLevel[peptideLevel$Subject == "Subject 2",]$V1
peptideLevel[peptideLevel$Subject == "Subject 3",]$`Subject 3` <- peptideLevel[peptideLevel$Subject == "Subject 3",]$V1
peptideLevel[peptideLevel$Subject == "Subject 4",]$`Subject 4` <- peptideLevel[peptideLevel$Subject == "Subject 4",]$V1

#write.table(peptideLevel, file = "ncMAPs.tsv", sep = "\t")

peptideLevel <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "ncMAP")
peptideLevel$Events <- factor(peptideLevel$Events, 
                              levels = c("Coding", "5`-UTR", "FS", "IR", "ncRNA", "3`-UTR", "IGR", "Unknown"))
peptideLevel <- peptideLevel[order(peptideLevel$Events, -peptideLevel$PSM1, -peptideLevel$PSM2, -peptideLevel$PSM3, -peptideLevel$PSM4),]

peptideLevel$Event <- 8
peptideLevel[peptideLevel$Events == "5`-UTR", ]$Event <- 7
peptideLevel[peptideLevel$Events == "FS", ]$Event <- 6
peptideLevel[peptideLevel$Events == "IR", ]$Event <- 5
peptideLevel[peptideLevel$Events == "ncRNA", ]$Event <- 4
peptideLevel[peptideLevel$Events == "3`-UTR", ]$Event <- 3
peptideLevel[peptideLevel$Events == "IGR", ]$Event <- 2
peptideLevel[peptideLevel$Events == "Unknown", ]$Event <- 1

peptideLevel[peptideLevel$PSM1 > 2, ]$PSM1 <- 3
peptideLevel[peptideLevel$PSM2 > 2, ]$PSM2 <- 3
peptideLevel[peptideLevel$PSM3 > 2, ]$PSM3 <- 3
peptideLevel[peptideLevel$PSM4 > 2, ]$PSM4 <- 3

peptideLevel$Read1 <- log2(peptideLevel$Read1+1)
peptideLevel$Read2 <- log2(peptideLevel$Read2+1)
peptideLevel$Read3 <- log2(peptideLevel$Read3+1)
peptideLevel$Read4 <- log2(peptideLevel$Read4+1)

nrow(peptideLevel[peptideLevel$Event == 8, ])
nrow(peptideLevel[peptideLevel$Event == 7, ])
nrow(peptideLevel[peptideLevel$Event == 6, ])
nrow(peptideLevel[peptideLevel$Event == 5, ])
nrow(peptideLevel[peptideLevel$Event == 4, ])
nrow(peptideLevel[peptideLevel$Event == 3, ])
nrow(peptideLevel[peptideLevel$Event == 2, ])
nrow(peptideLevel[peptideLevel$Event == 1, ])
nrow(peptideLevel[peptideLevel$asRNA1 == 1 | peptideLevel$asRNA2 == 1 | peptideLevel$asRNA3 == 1 | peptideLevel$asRNA4 == 1, ])

nrow(peptideLevel[peptideLevel$Event == 8 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 7 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 6 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 5 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 4 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 3 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 2 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 1 & peptideLevel$`Observed subjects` > 1, ])

pdf(file = "Events.pdf", width = 10, height = 6)
Heatmap(t(peptideLevel[, c("Event")]), width = 0.5, height = unit(0.3, "cm"), col = brewer.pal(n=8, name="Set1"),
        show_row_names = F, 
        cluster_columns = F, show_column_dend = F, heatmap_legend_param = list(
          title = "Event", labels = c("Coding", "5`-UTR", "FS", "IR", "ncRNA", "3`-UTR", "IGR", 
                                      "Unknown"),
          color_bar = "discrete", direction = "horizontal", ncol = 2, nrow = 4
        )) %v%
  Heatmap(t(peptideLevel[, c("PSM1", "PSM2", "PSM3", "PSM4")]), width = 0.5, height = unit(2.2, "cm"), 
          row_title = "PSM", row_title_rot = 90, show_row_names = F, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Blues", n=4), 
          cluster_columns = F, show_column_dend = F, cluster_rows = F, heatmap_legend_param = list(
            title = "Number of PSMs", at = c(0, 1, 2, 3), 
            labels = c("0", "1", "2", ">2"), legend_height = unit(1, "cm"),
            color_bar = "discrete", nrow = 1
          )) %v%
  Heatmap(t(peptideLevel[, c("Read1","Read2","Read3","Read4")]), width = 0.5, height = unit(2.2, "cm"),
          show_row_names = F,
          row_title = "Read", row_title_rot = 90, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Reds", n=8),
          cluster_columns = F, show_column_dend = F, cluster_rows = F, heatmap_legend_param = list(
            title = "log2(number of reads+1)",legend_height = unit(1, "cm"), nrow = 1, direction = "horizontal"
          )) %v%
  Heatmap(t(peptideLevel[, c("Mutation1","Mutation2","Mutation3","Mutation4")]), width = 0.5, height = unit(2.2, "cm"),
          show_row_names = F,
          cluster_columns = F, show_column_dend = F, cluster_rows = F, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          row_title = "Mutation", row_title_rot = 90,
          col = brewer.pal(name="Blues", n=5),
          heatmap_legend_param = list(
            title = "Number of mutations", at = c(0, 1, 2, 3, 4), 
            labels = c("0", "1", "2", "3", "4"), legend_height = unit(1, "cm"), 
            color_bar = "discrete", nrow = 1
          ))  %v%
  Heatmap(t(peptideLevel[, c("AS1","AS2","AS3","AS4")]), width = 0.5, height = unit(1.3, "cm"),
          show_row_names = F,
          row_title = "AS", row_title_rot = 90, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Reds", n=3),
          cluster_columns = F, show_column_dend = F, cluster_rows = F, heatmap_legend_param = list(
            title = "Is AS", at = c(0, 1), 
            labels = c("N", "Y"), legend_height = unit(1, "cm"),
            color_bar = "discrete", nrow = 1
          )) %v%
  Heatmap(t(peptideLevel[, c("asRNA1","asRNA2","asRNA3","asRNA4")]), width = 0.5, height = unit(1.3, "cm"),
          show_row_names = F,
          row_title = "asRNA", row_title_rot = 90, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Blues", n=3),
          cluster_columns = F, show_column_dend = F, cluster_rows = F,
          heatmap_legend_param = list(
            title = "Is asRNA", at = c(0, 1), 
            labels = c("N", "Y"), legend_height = unit(1, "cm"),
            color_bar = "discrete", nrow = 1
          )) 

dev.off()

## Compare with pXg
#1. Load pXg results
pXgdataS1 <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXgResults_p001_FDR.xlsx", sheet = "Subject 1")
pXgdataS2 <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXgResults_p001_FDR.xlsx", sheet = "Subject 2")
pXgdataS3 <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXgResults_p001_FDR.xlsx", sheet = "Subject 3")
pXgdataS4 <- read_excel(path = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXgResults_p001_FDR.xlsx", sheet = "Subject 4")

pXgdataS1 <- pXgdataS1[, -c(38:43)]
pXgdataS2 <- pXgdataS2[, -c(38:43)]
pXgdataS3 <- pXgdataS3[, -c(38:43)]
pXgdataS4 <- pXgdataS4[, -c(38:43)]

pXgdataS1$Subject <- "Subject 1"
pXgdataS2$Subject <- "Subject 2"
pXgdataS3$Subject <- "Subject 3"
pXgdataS4$Subject <- "Subject 4"

pXgallData <- rbind(pXgdataS1, pXgdataS2, pXgdataS3, pXgdataS4)

pXgsubData <- pXgallData[pXgallData$BestScore < 2, ]
pXgsubData$Type <- "WB"
pXgsubData[pXgsubData$BestScore < 2, ]$Type <- "WB"
pXgsubData[pXgsubData$BestScore < 0.5, ]$Type <- "SB"
pXgsubData$Class <- "Canonical"
pXgsubData[pXgsubData$IsCanonical == "FALSE", ]$Class <- "Noncanonical"

pXgsubDataUnique <- pXgsubData[pXgsubData$EventCount == 1, ]
pXgsubDataUnique <- pXgsubDataUnique[pXgsubDataUnique$GeneIDCount <= 1, ]
pXgsubDataUnique <- pXgsubDataUnique[pXgsubDataUnique$GenomicLociCount <= 1, ]

#2. Load pNovo3 results
pNovodataS1 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 1")
pNovodataS2 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 2")
pNovodataS3 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 3")
pNovodataS4 <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "Subject 4")

pNovodataS1 <- pNovodataS1[, -c(24:29)]
pNovodataS2 <- pNovodataS2[, -c(24:29)]
pNovodataS3 <- pNovodataS3[, -c(24:29)]
pNovodataS4 <- pNovodataS4[, -c(24:29)]

pNovodataS1 <- pNovodataS1[pNovodataS1$Score >= subject1_score | pNovodataS1$IsCanonical == TRUE, ]
pNovodataS2 <- pNovodataS2[pNovodataS2$Score >= subject2_score | pNovodataS2$IsCanonical == TRUE, ]
pNovodataS3 <- pNovodataS3[pNovodataS3$Score >= subject3_score | pNovodataS3$IsCanonical == TRUE, ]
pNovodataS4 <- pNovodataS4[pNovodataS4$Score >= subject4_score | pNovodataS4$IsCanonical == TRUE, ]

pNovodataS1$Subject <- "Subject 1"
pNovodataS2$Subject <- "Subject 2"
pNovodataS3$Subject <- "Subject 3"
pNovodataS4$Subject <- "Subject 4"

pNovoallData <- rbind(pNovodataS1, pNovodataS2, pNovodataS3, pNovodataS4)

pNovosubData <- pNovoallData[pNovoallData$BestScore < 2, ]
pNovosubData$Type <- "WB"
pNovosubData[pNovosubData$BestScore < 2, ]$Type <- "WB"
pNovosubData[pNovosubData$BestScore < 0.5, ]$Type <- "SB"
pNovosubData$Class <- "Canonical"
pNovosubData[pNovosubData$IsCanonical == "FALSE", ]$Class <- "Noncanonical"

pNovosubDataUnique <- pNovosubData[pNovosubData$EventCount == 1, ]
pNovosubDataUnique <- pNovosubDataUnique[pNovosubDataUnique$GeneIDCount <= 1, ]
pNovosubDataUnique <- pNovosubDataUnique[pNovosubDataUnique$GenomicLociCount <= 1, ]

# unique peptides in canonical peptides
pNovoPeptides <- pNovosubDataUnique[pNovosubDataUnique$Class == "Canonical", c("InferredPeptide", "Class", "Subject") ]
pXgPeptides <- pXgsubDataUnique[pXgsubDataUnique$Class == "Canonical", c("InferredPeptide", "Class", "Subject")  ]

subjectList <- c("Subject 1", "Subject 2", "Subject 3", "Subject 4")

compData <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(compData) <- c("PEAKS","pNovo3", "Common", "Subject", "Class")

for(subject in subjectList) {
  pNovoS <- pNovoPeptides[pNovoPeptides$Subject == subject, ]
  pXgS <- pXgPeptides[pXgPeptides$Subject == subject, ]
  
  pNovoS <- pNovoS[!duplicated(pNovoS$InferredPeptide), ]
  pXgS <- pXgS[!duplicated(pXgS$InferredPeptide), ]
  
  unionPeptides <- rbind(pNovoS, pXgS)
  unionPeptides <- unionPeptides[!duplicated(unionPeptides$InferredPeptide),]
  
  commonCount <- ((nrow(pNovoS) + nrow(pXgS))-nrow(unionPeptides))
  pXgCount <- nrow(pXgS) - commonCount
  pNovoCount <- nrow(pNovoS) - commonCount
  
  compData = rbind(compData, c(pXgCount, pNovoCount, commonCount, subject, "Canonical"))
}

# unique peptides in noncanonical peptides
pNovoPeptides <- pNovosubDataUnique[pNovosubDataUnique$Class == "Noncanonical", c("InferredPeptide", "Class", "Subject") ]
pXgPeptides <- pXgsubDataUnique[pXgsubDataUnique$Class == "Noncanonical", c("InferredPeptide", "Class", "Subject")  ]

for(subject in subjectList) {
  pNovoS <- pNovoPeptides[pNovoPeptides$Subject == subject, ]
  pXgS <- pXgPeptides[pXgPeptides$Subject == subject, ]
  
  pNovoS <- pNovoS[!duplicated(pNovoS$InferredPeptide), ]
  pXgS <- pXgS[!duplicated(pXgS$InferredPeptide), ]
  
  unionPeptides <- rbind(pNovoS, pXgS)
  unionPeptides <- unionPeptides[!duplicated(unionPeptides$InferredPeptide),]
  
  commonCount <- ((nrow(pNovoS) + nrow(pXgS))-nrow(unionPeptides))
  pXgCount <- nrow(pXgS) - commonCount
  pNovoCount <- nrow(pNovoS) - commonCount
  
  compData = rbind(compData, c(pXgCount, pNovoCount, commonCount, subject, "Noncanonical"))
}

colnames(compData) <- c("PEAKS","pNovo3", "Common", "Subject", "Class")
compData

a<-pXgsubDataUnique[pXgsubDataUnique$Subject == "Subject 2" & pXgsubDataUnique$Class == "Noncanonical", ]
nrow(a[!duplicated(a$InferredPeptide),])

## Check distinct sets
#ncMAP_PEAKS_pNovo3
PEAKSpNovo3Comp <- read_excel(path = "pXgResults_pNovo3_p001_FDR.xlsx", sheet = "ncMAP_PEAKS_pNovo3")

S1comp <- PEAKSpNovo3Comp[PEAKSpNovo3Comp$PSM1 > 0, ]
S2comp <- PEAKSpNovo3Comp[PEAKSpNovo3Comp$PSM2 > 0, ]
S3comp <- PEAKSpNovo3Comp[PEAKSpNovo3Comp$PSM3 > 0, ]
S4comp <- PEAKSpNovo3Comp[PEAKSpNovo3Comp$PSM4 > 0, ]

S1comp <- S1comp[(!duplicated(S1comp$InferredPeptide) & !duplicated(S1comp$InferredPeptide, fromLast = TRUE)),]
S2comp <- S2comp[(!duplicated(S2comp$InferredPeptide) & !duplicated(S2comp$InferredPeptide, fromLast = TRUE)),]
S3comp <- S3comp[(!duplicated(S3comp$InferredPeptide) & !duplicated(S3comp$InferredPeptide, fromLast = TRUE)),]
S4comp <- S4comp[(!duplicated(S4comp$InferredPeptide) & !duplicated(S4comp$InferredPeptide, fromLast = TRUE)),]

distinctPEAKS <- rbind(S1comp[S1comp$Search == "PEAKS", ], S2comp[S2comp$Search == "PEAKS", ], S3comp[S3comp$Search == "PEAKS", ], S4comp[S4comp$Search == "PEAKS", ])
distinctpNovo <- rbind(S1comp[S1comp$Search == "pNovo3", ], S2comp[S2comp$Search == "pNovo3", ], S3comp[S3comp$Search == "pNovo3", ], S4comp[S4comp$Search == "pNovo3", ])

peptideLevel <- distinctPEAKS

peptideLevel$Events <- factor(peptideLevel$Events, 
                              levels = c("Coding", "5`-UTR", "FS", "IR", "ncRNA", "3`-UTR", "IGR", "Unknown"))
peptideLevel <- peptideLevel[order(peptideLevel$Events, -peptideLevel$PSM1, -peptideLevel$PSM2, -peptideLevel$PSM3, -peptideLevel$PSM4),]

peptideLevel$Event <- 8
peptideLevel[peptideLevel$Events == "5`-UTR", ]$Event <- 7
peptideLevel[peptideLevel$Events == "FS", ]$Event <- 6
peptideLevel[peptideLevel$Events == "IR", ]$Event <- 5
peptideLevel[peptideLevel$Events == "ncRNA", ]$Event <- 4
peptideLevel[peptideLevel$Events == "3`-UTR", ]$Event <- 3
peptideLevel[peptideLevel$Events == "IGR", ]$Event <- 2
peptideLevel[peptideLevel$Events == "Unknown", ]$Event <- 1

peptideLevel[peptideLevel$PSM1 > 2, ]$PSM1 <- 3
peptideLevel[peptideLevel$PSM2 > 2, ]$PSM2 <- 3
peptideLevel[peptideLevel$PSM3 > 2, ]$PSM3 <- 3
peptideLevel[peptideLevel$PSM4 > 2, ]$PSM4 <- 3

peptideLevel$Read1 <- log2(peptideLevel$Read1+1)
peptideLevel$Read2 <- log2(peptideLevel$Read2+1)
peptideLevel$Read3 <- log2(peptideLevel$Read3+1)
peptideLevel$Read4 <- log2(peptideLevel$Read4+1)

nrow(peptideLevel[peptideLevel$Event == 8, ])
nrow(peptideLevel[peptideLevel$Event == 7, ])
nrow(peptideLevel[peptideLevel$Event == 6, ])
nrow(peptideLevel[peptideLevel$Event == 5, ])
nrow(peptideLevel[peptideLevel$Event == 4, ])
nrow(peptideLevel[peptideLevel$Event == 3, ])
nrow(peptideLevel[peptideLevel$Event == 2, ])
nrow(peptideLevel[peptideLevel$Event == 1, ])
nrow(peptideLevel[peptideLevel$asRNA1 == 1 | peptideLevel$asRNA2 == 1 | peptideLevel$asRNA3 == 1 | peptideLevel$asRNA4 == 1, ])

nrow(peptideLevel[peptideLevel$Event == 8 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 7 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 6 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 5 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 4 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 3 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 2 & peptideLevel$`Observed subjects` > 1, ])
nrow(peptideLevel[peptideLevel$Event == 1 & peptideLevel$`Observed subjects` > 1, ])

pdf(file = "Events_Distinct_pNovo3.pdf", width = 10, height = 6)
Heatmap(t(peptideLevel[, c("Event")]), width = 0.5, height = unit(0.3, "cm"), col = brewer.pal(n=8, name="Set1"),
        show_row_names = F, 
        cluster_columns = F, show_column_dend = F, heatmap_legend_param = list(
          title = "Event", labels = c("Coding", "5`-UTR", "FS", "IR", "ncRNA", "3`-UTR", "IGR", 
                                      "Unknown"),
          color_bar = "discrete", direction = "horizontal", ncol = 2, nrow = 4
        )) %v%
  Heatmap(t(peptideLevel[, c("PSM1", "PSM2", "PSM3", "PSM4")]), width = 0.5, height = unit(2.2, "cm"), 
          row_title = "PSM", row_title_rot = 90, show_row_names = F, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Blues", n=4), 
          cluster_columns = F, show_column_dend = F, cluster_rows = F, heatmap_legend_param = list(
            title = "Number of PSMs", at = c(0, 1, 2, 3), 
            labels = c("0", "1", "2", ">2"), legend_height = unit(1, "cm"),
            color_bar = "discrete", nrow = 1
          )) %v%
  Heatmap(t(peptideLevel[, c("Read1","Read2","Read3","Read4")]), width = 0.5, height = unit(2.2, "cm"),
          show_row_names = F,
          row_title = "Read", row_title_rot = 90, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Reds", n=8),
          cluster_columns = F, show_column_dend = F, cluster_rows = F, heatmap_legend_param = list(
            title = "log2(number of reads+1)",legend_height = unit(1, "cm"), nrow = 1, direction = "horizontal"
          )) %v%
  Heatmap(t(peptideLevel[, c("Mutation1","Mutation2","Mutation3","Mutation4")]), width = 0.5, height = unit(2.2, "cm"),
          show_row_names = F,
          cluster_columns = F, show_column_dend = F, cluster_rows = F, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          row_title = "Mutation", row_title_rot = 90,
          col = brewer.pal(name="Blues", n=5),
          heatmap_legend_param = list(
            title = "Number of mutations", at = c(0, 1, 2, 3, 4), 
            labels = c("0", "1", "2", "3", "4"), legend_height = unit(1, "cm"), 
            color_bar = "discrete", nrow = 1
          ))  %v%
  Heatmap(t(peptideLevel[, c("AS1","AS2","AS3","AS4")]), width = 0.5, height = unit(1.3, "cm"),
          show_row_names = F,
          row_title = "AS", row_title_rot = 90, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Reds", n=3),
          cluster_columns = F, show_column_dend = F, cluster_rows = F, heatmap_legend_param = list(
            title = "Is AS", at = c(0, 1), 
            labels = c("N", "Y"), legend_height = unit(1, "cm"),
            color_bar = "discrete", nrow = 1
          )) %v%
  Heatmap(t(peptideLevel[, c("asRNA1","asRNA2","asRNA3","asRNA4")]), width = 0.5, height = unit(1.3, "cm"),
          show_row_names = F,
          row_title = "asRNA", row_title_rot = 90, row_title_gp = gpar(fontsize = 11, fontface="bold"),
          col = brewer.pal(name="Blues", n=3),
          cluster_columns = F, show_column_dend = F, cluster_rows = F,
          heatmap_legend_param = list(
            title = "Is asRNA", at = c(0, 1), 
            labels = c("N", "Y"), legend_height = unit(1, "cm"),
            color_bar = "discrete", nrow = 1
          )) 

dev.off()





