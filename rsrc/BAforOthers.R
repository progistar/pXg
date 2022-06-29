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

setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/PreviousNCMAPs")
dataPC9 <- read_excel(path = "Qi_MCP2021.xlsx", sheet = "PC9 De Novo", skip = 1)
dataH1975 <- read_excel(path = "Qi_MCP2021.xlsx", sheet = "H1975 De Novo", skip = 1)
data3784Mel <- read_excel(path = "Qi_MCP2021.xlsx", sheet = "3784Mel De Novo", skip = 1)
data3785Mel <- read_excel(path = "Qi_MCP2021.xlsx", sheet = "3795Mel De Novo", skip = 1)
dataRA007 <- read_excel(path = "Qi_MCP2021.xlsx", sheet = "RA007 De Novo", skip = 1)
dataHCT116 <- read_excel(path = "Chen_DB_HCT116_NetMHCpan.xlsx", sheet = "1", skip = 1)
dataJurkat <- read_excel(path = "Chen_DB_Jurkat_NetMHCpan.xlsx", sheet = "1", skip = 1)

dataPC9 <- dataPC9[, c(1,2, length(dataPC9))]
dataH1975 <- dataH1975[, c(1,2, length(dataH1975))]
data3784Mel <- data3784Mel[, c(1,2, length(data3784Mel))]
data3785Mel <- data3785Mel[, c(1,2, length(data3785Mel))]
dataRA007 <- dataRA007[, c(1,2, length(dataRA007))]
dataHCT116 <- dataHCT116[, c(2, 3, length(dataHCT116))]
dataJurkat <- dataJurkat[, c(2, 3, length(dataJurkat))]

dataPC9$Sample <- "PC9 (Qi et al., MCP2021)"
dataH1975$Sample <- "H1975 (Qi et al., MCP2021)"
data3784Mel$Sample <- "3784Mel (Qi et al., MCP2021)"
data3785Mel$Sample <- "3785Mel (Qi et al., MCP2021)"
dataRA007$Sample <- "RA007 (Qi et al., MCP2021)"
dataHCT116$Sample <- "HCT116 (Chen et al., JASMS2021)"
dataJurkat$Sample <- "Jurkat E6-1 (Chen et al., JASMS2021)"

dataPC9 <- dataPC9[!duplicated(dataPC9[,1]), ]
dataH1975 <- dataH1975[!duplicated(dataH1975[,1]), ]
data3784Mel <- data3784Mel[!duplicated(data3784Mel[,1]), ]
data3785Mel <- data3785Mel[!duplicated(data3785Mel[,1]), ]
dataRA007 <- dataRA007[!duplicated(dataRA007[,1]), ]
dataHCT116 <- dataHCT116[!duplicated(dataHCT116[,1]), ]
dataJurkat <- dataJurkat[!duplicated(dataJurkat[,1]), ]

defNames <- c("Peptide", "Length", "BestScore", "Sample")
names(dataPC9) <- defNames
names(dataH1975) <- defNames
names(data3784Mel) <- defNames
names(data3785Mel) <- defNames
names(dataRA007) <- defNames
names(dataHCT116) <- defNames
names(dataJurkat) <- defNames

QiData <- rbind(dataPC9, dataH1975, data3784Mel, data3785Mel, dataRA007)
ChenData <- rbind(dataHCT116, dataJurkat)

## select 
QiData$Length <- as.character(QiData$Length)
QiData$Length <- factor(QiData$Length, levels=c("8", "9", "10", "11", "12", "13", "14", "15"))
QiData$BestScore <- log2(QiData$BestScore+1.5)

ChenData$Length <- as.character(ChenData$Length)
ChenData$Length <- factor(ChenData$Length, levels=c("8", "9", "10", "11", "12", "13", "14", "15"))
ChenData$BestScore <- log2(ChenData$BestScore+1.5)


subData <- QiData[QiData$Sample == "PC9 (Qi et al., MCP2021)", ]
#subData <- ChenData[ChenData$Sample == "Jurkat E6-1 (Chen et al., JASMS2021)", ]

mapPlot <- ggplot(data=subData, aes(x=Length, y=BestScore, fill=Length)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  ggtitle(substitute(paste("HLA binding affinity for ", "noncanonical ", "peptides (PC9)", sep=""))) +
  scale_y_continuous(breaks=seq(from=0, to=7, by= 1)) +
  labs(y= TeX("$Log_{2}$(1.5+%Rank)"), x = "Peptide length") +
  annotate("text", label = nrow(subData[subData$Length == 8,]), x =1, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 9,]), x =2, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 10,]), x =3, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 11,]), x =4, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 12,]), x =5, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 13,]), x =6, y = 7, size = 8, family="serif", colour = "black") +
  annotate("text", label = nrow(subData[subData$Length == 14,]), x =7, y = 7, size = 8, family="serif", colour = "black") +
  #annotate("text", label = nrow(subData[subData$Length == 15,]), x =8, y = 7, size = 8, family="serif", colour = "black")
  geom_hline(yintercept=log2(1.5+2.0), linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=log2(1.5+0.5), linetype="dashed", color = "black", size=1)
  

mapPlot

ggsave("BA.DeNovo.Qi.PC9.png", plot = mapPlot, width = 10, height = 8, units = "in", dpi = 600)
.