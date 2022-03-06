install.packages("ggplot2")

library(ggplot2)
library(RColorBrewer)
library(readxl)

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
setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/3.withCalibrationAddScanNum")

pResult <- read_excel(path = "Laumont_Results.xlsx", sheet = "Deamidation")
pResult <- data.frame(pResult)

g <- ggplot(data = pResult, aes(x=DeltaALCScore, fill=Class)) +
  theme_bw() +
  geom_histogram(position = "dodge", bins = 15) +
  theme(text = element_text(size=25)) +
  ggtitle("Delta score distribution of deamidated spectra") +
  labs(y="Spectra", x = "Delta ALC score")
ggsave("subjectS1.deltaDeamiScore.png", plot = g, width = 10, height = 8, units = "in", dpi = 300)

nrow(pResult[pResult$Class == "same" & pResult$DeltaALCScore == 1, ])

