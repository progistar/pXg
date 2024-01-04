
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

setwd("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy")

generate_scatter_plot <- function(data) {
  data <- data[!is.na(data$`E-value`), ]
  scatterPhredPlot <- ggplot(data=data, aes(x=InferredPeptide, y=-log(`E-value`,10), group=Origin, color=Origin)) +
    theme_bw() +
    scale_color_brewer(palette = "Set3") +
    geom_point(color='black', shape=21, size=3, aes(fill=Origin)) + 
    scale_fill_brewer(palette="Set3") +
    theme(text = element_text(size=20), strip.background = element_blank(), legend.title = element_blank(),
          plot.margin = margin(0.1, 0.1, 0, -0.3, "in"),
          axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
    #scale_x_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
    #scale_y_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
    staticThemeRightTop +
    labs(y= "", x = "")
  
  return (scatterPhredPlot)
}

unknown <- read_excel(path = "UnmappedReads.xlsx", sheet = "Unknown")


scatter_plot <- generate_scatter_plot(unknown)
scatter_plot


ggsave("FigureS2.png", plot = scatter_plot, width = 20, height = 6, units = "in", dpi = 600)

