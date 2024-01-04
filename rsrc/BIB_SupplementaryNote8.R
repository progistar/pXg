
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
  scatterPhredPlot <- ggplot(data=data, aes(x=Log2Reads, y=SA, group=InferredPeptide, color=InferredPeptide)) +
    theme_bw() +
    scale_color_brewer(palette = "Set3") +
    geom_point(color='black', shape=21, size=3, aes(fill=InferredPeptide)) + 
    scale_fill_brewer(palette="Set3") +
    theme(text = element_text(size=20), strip.background = element_blank(), legend.title = element_blank(),
          plot.margin = margin(0.1, 0.1, 0, -0.3, "in")) +
    #scale_x_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
    #scale_y_continuous(breaks = c(25,50,75,100), limits = c(0,100)) +
    staticThemeRightTop +
    labs(y= "", x = "") +
    facet_grid(cols = vars(`Events`), scales = "free")
  
  return (scatterPhredPlot)
}

fdr_5_res_noncanonical <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "Noncanonical")
retro_viruses <- fdr_5_res_noncanonical[fdr_5_res_noncanonical$GeneNames == "Squirrel monkey retrovirus", ]
retro_viruses[retro_viruses$InferredPeptide == "SLYPNLATL", ]$Events <- "Gag protein"
retro_viruses[retro_viruses$InferredPeptide == "SLLDGVKTLV", ]$Events <- "Gag protein"
retro_viruses[retro_viruses$InferredPeptide == "SLLDGVKTL", ]$Events <- "Gag protein"
retro_viruses[retro_viruses$InferredPeptide == "LLDGVKTLV", ]$Events <- "Gag protein"
retro_viruses[retro_viruses$InferredPeptide == "LLPPAGVMA", ]$Events <- "Gag protein"

retro_viruses[retro_viruses$InferredPeptide == "KLFSGILDTGA", ]$Events <- "Intergenic region"
retro_viruses[retro_viruses$InferredPeptide == "LTYEKTLAA", ]$Events <- "Intergenic region"

retro_viruses[retro_viruses$InferredPeptide == "VLAHQPFNL", ]$Events <- "Pol protein"

retro_viruses[retro_viruses$InferredPeptide == "ALDISNPSL", ]$Events <- "Envelope protein"
retro_viruses[retro_viruses$InferredPeptide == "QLINDVQAL", ]$Events <- "Envelope protein"
retro_viruses[retro_viruses$InferredPeptide == "SLLGKPIQI", ]$Events <- "Envelope protein"

nrow(retro_viruses[retro_viruses$Events == "Gag protein", ])
nrow(retro_viruses[retro_viruses$Events == "Intergenic region", ])
nrow(retro_viruses[retro_viruses$Events == "Pol protein", ])
nrow(retro_viruses[retro_viruses$Events == "Envelope protein", ])

scatter_plot <- generate_scatter_plot(retro_viruses)
scatter_plot
ggsave("FigureS2_SA_Log2Reads_ScatterPlot.png", plot = scatter_plot, width = 18, height = 6, units = "in", dpi = 600)
nrow(retro_viruses)
retro_viruses <- retro_viruses[order(-retro_viruses$SA),]
viruse_peptide_level <- retro_viruses[!duplicated(retro_viruses$InferredPeptide), ]


