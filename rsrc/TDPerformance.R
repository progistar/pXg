
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

setwd("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Performances/MockReadMethods")

staticThemeLeftTop <- theme(text = element_text(size=25, color = "black")
                            , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                            legend.justification = c("right", "top"),
                            legend.position= c(.22, .98),
                            legend.text = element_text(size=20, color = "black"),
                            legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeRightTop <- theme(text = element_text(size=25, color = "black")
                             , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"),
                             legend.justification = c("right", "top"),
                             legend.position= c(.98, .98),
                             legend.text = element_text(size=20, color = "black"),
                             legend.box.background = element_rect(linetype = 1, size = 1))

staticThemeNone <- theme(text = element_text(size=25, color = "black")
                         , axis.text.x = element_text(size=20, color="black"), axis.text.y = element_text(size=20, color="black"))

PSDG_data <-  read_excel(path = "MockMethodComparison.xlsx", sheet = "six-frame-mean-std-group")
PSDG_data <- as.data.frame(PSDG_data)
PSDG_data[PSDG_data$RNACutoff == T, ]$RNACutoff <- "with read cutoff"
PSDG_data[PSDG_data$RNACutoff == F, ]$RNACutoff <- "w/o read cutoff"

PSDG_plot <- ggplot(PSDG_data, aes(x=Rank, y=Mean, fill=Label)) +
  geom_bar(position=position_dodge2(preserve = "single", padding = 0), stat="identity",
           colour='black') +
  geom_errorbar(aes(ymin=Mean-Std, ymax=Mean+Std), width=.2, 
                position=position_dodge(.9)) +
  scale_y_continuous(breaks = seq(from=0, to=0.5, by=0.05), limits=c(0,0.5), expand = expansion(mult = c(0, .05))) +
  theme_bw() +
  staticThemeNone +
  scale_fill_manual(values=c(blueC, redC, greenC), name ="Distance") +
  labs(y="Mean of decoy PSM ratio", x = "Rank") +
  facet_grid(. ~ RNACutoff)

PSDG_plot

ggsave("PSDG_group.png", plot = PSDG_plot, width = 11, height = 8, units = "in", dpi = 300)

## real spectra
Real_data <-  read_excel(path = "MockMethodComparison.xlsx", sheet = "six-frame-real")
Real_data <- as.data.frame(Real_data)
Real_data[Real_data$RNACutoff == T, ]$RNACutoff <- "with read cutoff"
Real_data[Real_data$RNACutoff == F, ]$RNACutoff <- "w/o read cutoff"

Real_plot <- ggplot(Real_data, aes(x=Rank, y=DecoyRatio)) +
  geom_bar(position=position_dodge2(preserve = "single", padding = 0), stat="identity",
           colour='black') +
  scale_x_continuous(breaks = seq(from=1, to=10, by = 1)) +
  scale_y_continuous(breaks = seq(from=0, to=0.5, by=0.05), limits=c(0,0.5), expand = expansion(mult = c(0, .05))) +
  theme_bw() +
  staticThemeNone +
  labs(y="Decoy PSM ratio", x = "Rank") +
  facet_grid(. ~ RNACutoff)

Real_plot

ggsave("Real.png", plot = Real_plot, width = 11, height = 8, units = "in", dpi = 300)


y = c(0.027124655, 0.141690009, 0.250847458, 0.301310044, 0.363005051, 0.39019337, 0.409701493, 0.386115445,0.396788991, 0.392912173)

log(y, 10)
