library(ggplot2)
library(RColorBrewer)

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


setwd("C:/Users/progi/Desktop/Projects/pXg/Laumont_NatCommun2016/Results")

scores <- read.csv(file = "ScoreDist.txt", header = T, sep = "\t", as.is = as.double())

tdPlot <- ggplot(data=scores, aes(x=Score, fill=Class)) +
  scale_fill_manual(values=c(blueC, redC)) +
  theme_bw() +
  geom_bar(aes(y= (..count..), x=Score), position=position_dodge(), width = 0.5) +
  #geom_text(stat='count', aes(label = scales::comma(..count.., accuracy = 1)), vjust=-0.5, hjust=c(1.0,-0.5,1.0,-0.5,1.0,-0.5,1.0,-0.5), family = "serif", size = 6) +
  theme(text = element_text(size=25)) +
  labs(y="PSMs", x = "ALC")

tdPlot
