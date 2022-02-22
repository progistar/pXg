library(ggplot2)
library(RColorBrewer)

setwd("C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\1.M_pXg")
data <- read.csv(file = "PeptideAnnotationM_5ppm_002.rep1.rank10.pXg.pval.dist", header = T, sep="\t", as.is = as.double())
data$PeptideLength <- as.character(data$PeptideLength)
data$Mock <- log2(data$Mock+1)
data$Experiment <- log2(data$Experiment+1)
data$PeptideLength <- factor(data$PeptideLength, levels = c('8','9','10','11','12','13', '14', '15'))
data <- data[data$ReadCount <= 20, ]

distPlot <- ggplot(data = data, aes(x=ReadCount, y=Mock, fill=PeptideLength)) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  geom_line(size = 1, aes(color=PeptideLength)) +
  xlab("Mock Reads") +
  ylab("Log2(number of peptides + 1)") +
  scale_x_continuous(breaks=seq(from=0, to=20, by = 1)) +
  geom_point()

distPlot



dataAll <- read.csv(file = "subjectM.5ppm.002.rep1.all.pXg", header = T, sep="\t", as.is = as.double())
dataMax <- read.csv(file = "subjectM.5ppm.002.rep1.max.pXg", header = T, sep="\t", as.is = as.double())

dataAll <- dataAll[!duplicated(dataAll[,c('Fraction', 'Source.File', 'Scan')]),]
dataMax <- dataMax[!duplicated(dataMax[,c('Fraction', 'Source.File', 'Scan')]),]

dataAll <- dataAll[!duplicated(dataAll[,c('InferredPeptide')]),]
dataMax <- dataMax[!duplicated(dataMax[,c('InferredPeptide')]),]

dataAll$Methods <- "Use all"
dataMax$Methods <- "Use max-one"

data <- rbind(dataAll, dataMax)

topRankedPlot <- ggplot(data=data, aes(x=Length, fill=Methods)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  geom_histogram(binwidth = 0.5, position=position_dodge()) +
  scale_x_continuous(breaks=seq(from=8, to=15, by = 1)) +
  theme(text = element_text(size=15), plot.title = element_text(family = "serif", face = "bold", hjust = 0, size = 15, color = "darkblue")) +
  labs(x="Peptide length", y="Peptides")
topRankedPlot

