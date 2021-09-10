library(ggplot2)

allDeNovoCSVFilePath <- "C:\\Users\\progi\\Desktop\\Projects\\pXg\\ProteomeTools\\all de novo candidates.tsv.annot"

data <- read.csv(file = allDeNovoCSVFilePath, header = T, sep = "\t",  as.is = as.double())

minPeptLen <- min(data$MaxQuantPeptLength)
maxPeptLen <- max(data$MaxQuantPeptLength)

## If you want to set free, set agglomeratePeptLen to maxPeptLen
agglomeratePeptLen <- 20


peptLenRange  <- c(minPeptLen:maxPeptLen)
counts        <-  c(minPeptLen:agglomeratePeptLen) * 0
types         <-  c(minPeptLen:agglomeratePeptLen) * 0

for(peptLen in peptLenRange) {
  index <- peptLen - minPeptLen + 1
  
  ## Index limitation
  if(peptLen > agglomeratePeptLen) index <- agglomeratePeptLen - minPeptLen + 1

    counts[index] <- nrow(data[data$MaxQuantPeptLength == peptLen, ]) + counts[index]
  types[index] <- "PSMs in MaxQuant"
}

MaxQuantByPeptLen <- data.frame(PeptideLength = c(minPeptLen:agglomeratePeptLen),
                                PSMs =  counts,
                                Type = types)

counts        <-  c(minPeptLen:agglomeratePeptLen) * 0
types         <-  c(minPeptLen:agglomeratePeptLen) * 0
for(peptLen in peptLenRange) {
  index <- peptLen - minPeptLen + 1
  
  ## Index limitation
  if(peptLen > agglomeratePeptLen) index <- agglomeratePeptLen - minPeptLen + 1
  
  counts[index] <- nrow(data[data$MaxQuantPeptLength == peptLen & data$Aggrement == "true", ]) + counts[index]
  types[index] <- "Agreements in N-ranked"
}

nRankedByPeptLen <- data.frame(PeptideLength = c(minPeptLen:agglomeratePeptLen),
                               PSMs = counts,
                               Type = types)

counts        <-  c(minPeptLen:agglomeratePeptLen) * 0
types         <-  c(minPeptLen:agglomeratePeptLen) * 0
for(peptLen in peptLenRange) {
  index <- peptLen - minPeptLen + 1
  
  ## Index limitation
  if(peptLen > agglomeratePeptLen) index <- agglomeratePeptLen - minPeptLen + 1
  
  counts[index] <- nrow(data[data$MaxQuantPeptLength == peptLen & data$Aggrement == "true" & data$DeltaRank == 0, ]) + counts[index]
  types[index] <- "Agreements in Top-ranked"
}

topRankedByPeptLen <- data.frame(PeptideLength = c(minPeptLen:agglomeratePeptLen),
                               PSMs = counts,
                               Type = types)

dataByPeptLen <- rbind(MaxQuantByPeptLen, nRankedByPeptLen)
dataByPeptLen <- rbind(dataByPeptLen, topRankedByPeptLen)

sum(dataByPeptLen[dataByPeptLen$Type == "PSMs in MaxQuant", ]$PSMs)


topRankedPlot <- ggplot(data=dataByPeptLen[dataByPeptLen$Type != "Agreements in N-ranked", ], aes(x=PeptideLength, y=PSMs, fill=Type)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_continuous(breaks=minPeptLen:agglomeratePeptLen) +
  ggtitle("Agreements in PSM-level between MaxQuant and Top-ranked PEAKS") +
  theme(text = element_text(size=15), plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 15, color = "darkblue")) +
  labs(x="Peptide length")


nRankedPlot <- ggplot(data=dataByPeptLen[dataByPeptLen$Type != "Agreements in Top-ranked", ], aes(x=PeptideLength, y=PSMs, fill=Type)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_continuous(breaks=minPeptLen:agglomeratePeptLen) +
  ggtitle("Agreements in PSM-level between MaxQuant and N-ranked PEAKS") +
  theme(text = element_text(size=15), plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 15, color = "darkblue")) +
  labs(x="Peptide length")

cowplot::plot_grid(topRankedPlot, nRankedPlot, labels = "AUTO", nrow = 2, align = "v")

## for recall
topRankedByPeptLen$PSMs / MaxQuantByPeptLen$PSMs 
nRankedByPeptLen$PSMs / MaxQuantByPeptLen$PSMs 


## score distribution

#### Score distribution
#data$MaxQuantPeptLength == 7
#data$Aggrement == "true"
#data$DeltaRank > 0
# show data$DeltaScore

data$DeltaRank <- as.double(as.character(data$DeltaRank))
data$DeltaScore <- as.double(as.character(data$DeltaScore))
data$MaxQuantPeptLength <- as.double(as.character(data$MaxQuantPeptLength))
subData <- data[data$MaxQuantPeptLength > 7 & data$MaxQuantPeptLength < 14 & data$DeltaRank > 0 & data$Aggrement == "true", ]
subData$MaxQuantPeptLength <- factor(subData$MaxQuantPeptLength, levels = c('8','9','10','11','12','13'))

topRankedBoxPlot <- ggplot(data=subData, aes(x=MaxQuantPeptLength, y=DeltaScore, fill = as.character(MaxQuantPeptLength))) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_violin(position = "dodge") +
  ggtitle("Agreements in PSM-level between MaxQuant and N-ranked PEAKS") +
  theme(text = element_text(size=15), plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 15, color = "darkblue"),
        legend.position = "none") +
  labs(x="Peptide length", y = "Delta de novo score")

topRankedBoxPlot

nrow(subData)
nrow(subData[subData$DeltaScore < 5, ])

sum(subData$DeltaRank)
sum(subData[subData$DeltaScore < 5, ]$DeltaRank)

