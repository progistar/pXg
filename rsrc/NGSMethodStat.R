library(ggplot2)
library(RColorBrewer)

setwd("C:\\Users\\progi\\Desktop\\Projects\\pXg\\MockTest")
data <- read.csv(file = "subjectM.5ppm.002.rep1.pXg.ngs.max.stat", header = T, sep="\t", as.is = as.double())

data <- data[data$ReadCount <= 20, ]

distPlot <- ggplot(data = data, aes(x=ReadCount, y=Target, group=(PeptideLength))) +
  theme_bw() +
  scale_x_continuous(breaks=seq(from=1, to=20, by = 1)) +
  geom_line(size = 0.5) +
  geom_point()

distPlot


all
sub
all = sum(data[data$PeptideLength == 8, ]$Decoy)
sub = sum(data[data$PeptideLength == 8 & data$ReadCount >= 20, ]$Decoy)
p <- sub/all
p
