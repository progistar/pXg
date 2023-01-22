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

staticThemeRightTop <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                             legend.justification = c("right", "top"),
                             legend.position= c(.95, .98),
                             legend.text = element_text(size=20, color = "black"),
                             legend.box.background = element_rect(linetype = 1, size = 1))
setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXg_010/")

s1 <- read.csv(file = "S1.RAW.PEAKS.pxg.fdr.dist", header = T, sep = "\t", as.is = as.double())
s2 <- read.csv(file = "S2.RAW.PEAKS.pxg.fdr.dist", header = T, sep = "\t", as.is = as.double())
s3 <- read.csv(file = "S3.RAW.PEAKS.pxg.fdr.dist", header = T, sep = "\t", as.is = as.double())
s4 <- read.csv(file = "S4.RAW.PEAKS.pxg.fdr.dist", header = T, sep = "\t", as.is = as.double())

s1$Subject <- "Subject 1"
s2$Subject <- "Subject 2"
s3$Subject <- "Subject 3"
s4$Subject <- "Subject 4"

samples <- rbind(s1, s2, s3, s4)

scores <- samples
targets <- scores[, c("Score", "cTargetCount", "Subject")]
decoys <- scores[, c("Score", "cDecoyCount", "Subject")]
colnames(targets) <- c("Score", "PSMs", "Subject")
colnames(decoys) <- c("Score", "PSMs", "Subject")
targets$Label <- "Target"
decoys$Label <- "Decoy"
canonical <- rbind(targets, decoys)

targets <- scores[, c("Score", "ncTargetCount", "Subject")]
decoys <- scores[, c("Score", "ncDecoyCount", "Subject")]
colnames(targets) <- c("Score", "PSMs", "Subject")
colnames(decoys) <- c("Score", "PSMs", "Subject")
targets$Label <- "Target"
decoys$Label <- "Decoy"
noncanonical <- rbind(targets, decoys)

canonical$Class <- "Canonical"
noncanonical$Class <- "Nonanonical"

scoreDist <- rbind(canonical, noncanonical)
scoreDist$Label <- factor(scoreDist$Label, levels = c("Target", "Decoy"))

tdPlot <- ggplot(data=scoreDist, aes(x=Score, y=PSMs, fill=Label)) +
  scale_fill_manual(values=c(blueC, redC)) +
  theme_bw() +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(y="PSM", x = "Average local confidence (ALC)") +
  facet_grid(rows = vars(`Class`), cols = vars(`Subject`), scales = "free") +
  staticThemeRightTop

plot(tdPlot)

ggsave("scoreDist.png", plot = tdPlot, width = 20, height = 8, units = "in", dpi = 300)

