
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
fdr_5_res_noncanonical <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "Noncanonical")
fdr_5_res_canonical <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "Canonical")

#THP1-1 + THP1-2 + THP1-3 = THP1
fdr_5_res_canonical[fdr_5_res_canonical$Sample == "THP1-1"|
                    fdr_5_res_canonical$Sample == "THP1-2"|
                    fdr_5_res_canonical$Sample == "THP1-3",]$Sample <- "THP1"

fdr_5_res_noncanonical[fdr_5_res_noncanonical$Sample == "THP1-1"|
                         fdr_5_res_noncanonical$Sample == "THP1-2"|
                         fdr_5_res_noncanonical$Sample == "THP1-3",]$Sample <- "THP1"

get_z_score_norm <- function (data, sample) {
  
  data <- data[data$Sample == sample, ]
  pivot_data <- data
  
  pivot_data <- pivot_data[!duplicated(pivot_data$InferredPeptide), ]
  mean_rna <- mean(pivot_data$Log2Reads)
  sd_rna <- sd(pivot_data$Log2Reads)
  
  data$NormLog2Reads <- (data$Log2Reads - mean_rna) / sd_rna
  
  return(data[!duplicated(data$InferredPeptide),])
}

fdr_5_res <- rbind(fdr_5_res_noncanonical, fdr_5_res_canonical)

THP1 <- get_z_score_norm(fdr_5_res, "THP1")
B_LCL1 <- get_z_score_norm(fdr_5_res, "B-LCL1")
B_LCL2 <- get_z_score_norm(fdr_5_res, "B-LCL2")
B_LCL3 <- get_z_score_norm(fdr_5_res, "B-LCL3")
B_LCL4 <- get_z_score_norm(fdr_5_res, "B-LCL4")
DOHH2 <- get_z_score_norm(fdr_5_res, "DOHH2")
HBL1 <- get_z_score_norm(fdr_5_res, "HBL1")
SUDHL4 <- get_z_score_norm(fdr_5_res, "SUDHL4")

t.test(SUDHL4[SUDHL4$IsCanonical == T, ]$NormLog2Reads,
       SUDHL4[SUDHL4$IsCanonical == F, ]$NormLog2Reads)

t.test(B_LCL1[B_LCL1$IsCanonical == T, ]$NormLog2Reads,
       B_LCL1[B_LCL1$IsCanonical == F, ]$NormLog2Reads)

t.test(B_LCL2[B_LCL2$IsCanonical == T, ]$NormLog2Reads,
       B_LCL2[B_LCL2$IsCanonical == F, ]$NormLog2Reads)

t.test(B_LCL3[B_LCL3$IsCanonical == T, ]$NormLog2Reads,
       B_LCL3[B_LCL3$IsCanonical == F, ]$NormLog2Reads)

t.test(B_LCL4[B_LCL4$IsCanonical == T, ]$NormLog2Reads,
       B_LCL4[B_LCL4$IsCanonical == F, ]$NormLog2Reads)

t.test(DOHH2[DOHH2$IsCanonical == T, ]$NormLog2Reads,
       DOHH2[DOHH2$IsCanonical == F, ]$NormLog2Reads)

t.test(HBL1[HBL1$IsCanonical == T, ]$NormLog2Reads,
       HBL1[HBL1$IsCanonical == F, ]$NormLog2Reads)


### RNA plot ###
Total_sample <- rbind(THP1, B_LCL1,B_LCL2, B_LCL3, B_LCL4, DOHH2, HBL1, SUDHL4)
selectedColor <- c(brewer.pal(n = 8, name = "Set3")[5],
                   brewer.pal(n = 8, name = "Set3")[4])

Total_sample$Class <- "Canonical"
Total_sample[Total_sample$IsCanonical == F, ]$Class <- "Noncanonical"
Total_sample$Report <- "Known"
Total_sample[Total_sample$NumberOfMappedSources == 0, ]$Report <- "Novel"
Total_sample$Report <- factor(Total_sample$Report, levels = c("Known", "Novel"))

mdPlot <- ggplot(data=Total_sample, aes(x=factor(Class), y=NormLog2Reads, fill=Report)) +
  theme_bw() +
  scale_fill_manual(values = selectedColor) +
  #scale_fill_brewer(palette = "Set3") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  scale_y_continuous(n.breaks = 10, limits = c(-5, 6)) +
  staticThemeTop +
  theme(legend.title = element_blank(),
        plot.margin = margin(0.1, 0.1, -0.0, -0.3, "in")) +
  labs(y= "", x="")
mdPlot

t.test(Total_sample[Total_sample$Report == "Known" & Total_sample$Class != "Canonical", ]$NormLog2Reads,
       Total_sample[Total_sample$Report == "Novel" & Total_sample$Class != "Canonical", ]$NormLog2Reads)

ggsave("Figure5_RNAexp.png", plot = mdPlot, width = 4.5, height = 6, units = "in", dpi = 600)

## Mutant vs Wild
Total_sample$Mutant <- "Wild"
Total_sample[Total_sample$Mutations != "-", ]$Mutant <- "Mutant"
mdPlot <- ggplot(data=Total_sample[Total_sample$IsCanonical == F, ], aes(x=factor(Mutant), y=NormLog2Reads, fill=Report)) +
  theme_bw() +
  scale_fill_manual(values = selectedColor) +
  #scale_fill_brewer(palette = "Set3") +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  scale_y_continuous(n.breaks = 10, limits = c(-5, 6)) +
  staticThemeTop +
  theme(legend.title = element_blank(),
        plot.margin = margin(0.1, 0.1, -0.0, -0.3, "in")) +
  labs(y= "", x="")
mdPlot

t.test(Total_sample[Total_sample$Report == "Novel" & Total_sample$Class != "Canonical" & Total_sample$Mutant == "Wild", ]$NormLog2Reads,
       Total_sample[Total_sample$Report == "Known" & Total_sample$Class != "Canonical" & Total_sample$Mutant == "Wild", ]$NormLog2Reads)

t.test(Total_sample[Total_sample$Report == "Novel" & Total_sample$Class != "Canonical" & Total_sample$Mutant == "Mutant", ]$NormLog2Reads,
       Total_sample[Total_sample$Report == "Known" & Total_sample$Class != "Canonical" & Total_sample$Mutant == "Mutant", ]$NormLog2Reads)


write.table(Total_sample[Total_sample$IsCanonical == F, ], file = "Noncanonical_NormLog2Reads.tsv", append = FALSE, quote = F, sep = "\t",
            na = "NA", dec = ".", row.names = F,
            col.names = TRUE)

write.table(Total_sample[Total_sample$IsCanonical == T, ], file = "Canonical_NormLog2Reads.tsv", append = FALSE, quote = F, sep = "\t",
            na = "NA", dec = ".", row.names = F,
            col.names = TRUE)

### noncanonical ###
Total_sample$key <- paste(Total_sample$BestHLAType,Total_sample$InferredPeptide,sep="_")
Total_sample[!is.na(Total_sample$`Immunogenicity 3.0`), ]$`Immunogenicity 3.0` <- Total_sample[!is.na(Total_sample$`Immunogenicity 3.0`), ]$`Immunogenicity 3.0` - 0.5
unique_list <- Total_sample[!duplicated(Total_sample$key),]
#unique_list <- unique_list[unique_list$GeneNames == "Squirrel monkey retrovirus", ]

immune_pred <- data.frame(matrix(nrow = 0, ncol = 4))
deepHLApan <- data.frame(matrix(nrow = 0, ncol = 3))
deepNetBim <- data.frame(matrix(nrow = 0, ncol = 3))
immunegenicity <- data.frame(matrix(nrow = 0, ncol = 3))

deepHLApan <- data.frame(unique_list[!is.na(unique_list$DeepHLApan), ]$Class, 
                         unique_list[!is.na(unique_list$DeepHLApan), ]$DeepHLApan,
                         unique_list[!is.na(unique_list$DeepHLApan), ]$Events)
deepHLApan$Tool <- "DeepHLApan"
colnames(deepHLApan) <- c("Class","Immunogenicity", "Events", "Tool")

deepNetBim <- data.frame(unique_list[!is.na(unique_list$DeepNetBim), ]$Class, 
                unique_list[!is.na(unique_list$DeepNetBim), ]$DeepNetBim,
                unique_list[!is.na(unique_list$DeepNetBim), ]$Events)
deepNetBim$Tool <- "DeepNetBim"
colnames(deepNetBim) <-  c("Class","Immunogenicity", "Events", "Tool")

immunegenicity <- data.frame(unique_list[!is.na(unique_list$`Immunogenicity 3.0`), ]$Class, 
                    unique_list[!is.na(unique_list$`Immunogenicity 3.0`), ]$`Immunogenicity 3.0`,
                    unique_list[!is.na(unique_list$`Immunogenicity 3.0`), ]$Events)
immunegenicity$Tool <- "Immunogenicity 3.0"
colnames(immunegenicity) <- c("Class","Immunogenicity", "Events", "Tool")

immune_pred <- rbind(deepHLApan, deepNetBim)
immune_pred <- rbind(immune_pred, immunegenicity)
colnames(immune_pred) <- c("Class","Immunogenicity", "Events", "Tool")
immune_pred$Class <- factor(immune_pred$Class, levels = c("Canonical", "Noncanonical"))

mdPlot <- ggplot(data=immune_pred, aes(x=factor(Tool), y=Immunogenicity, fill=Class)) +
  theme_bw() +
  scale_fill_manual(values=c(brewer.pal(n = 10, name = "Set3")[5], brewer.pal(n = 10, name = "Set3")[4])) +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  scale_y_continuous(n.breaks = 6, limits = c(-1, 1.5)) +
  staticThemeTop +
  theme(legend.title = element_blank(),
        plot.margin = margin(0.1, 0.1, -0.0, -0.3, "in")) +
  labs(y= "", x="")
mdPlot



t.test(immune_pred[immune_pred$Class == "Canonical" & immune_pred$Tool == "DeepHLApan", ]$Immunogenicity,
       immune_pred[immune_pred$Class == "Noncanonical" & immune_pred$Tool == "DeepHLApan", ]$Immunogenicity)
nrow(immune_pred[immune_pred$Class == "Canonical" & immune_pred$Tool == "DeepHLApan", ])
nrow(immune_pred[immune_pred$Class == "Noncanonical" & immune_pred$Tool == "DeepHLApan", ])

t.test(immune_pred[immune_pred$Class == "Canonical" & immune_pred$Tool == "DeepNetBim", ]$Immunogenicity,
       immune_pred[immune_pred$Class == "Noncanonical" & immune_pred$Tool == "DeepNetBim", ]$Immunogenicity)
nrow(immune_pred[immune_pred$Class == "Canonical" & immune_pred$Tool == "DeepNetBim", ])
nrow(immune_pred[immune_pred$Class == "Noncanonical" & immune_pred$Tool == "DeepNetBim", ])

t.test(immune_pred[immune_pred$Class == "Canonical" & immune_pred$Tool == "Immunogenicity 3.0", ]$Immunogenicity,
       immune_pred[immune_pred$Class == "Noncanonical" & immune_pred$Tool == "Immunogenicity 3.0", ]$Immunogenicity)
nrow(immune_pred[immune_pred$Class == "Canonical" & immune_pred$Tool == "Immunogenicity 3.0", ])
nrow(immune_pred[immune_pred$Class == "Noncanonical" & immune_pred$Tool == "Immunogenicity 3.0", ])

ggsave("Figure6_Immunogenicity.png", plot = mdPlot, width = 8, height = 7, units = "in", dpi = 600)

canonical_list <- Total_sample[Total_sample$IsCanonical == T, ]
unique_canonical_list <- canonical_list[!duplicated(canonical_list$key), ]

boxplot(unique_canonical_list$DeepNetBim, unique_canonical_list$DeepHLApan,
        unique_canonical_list$`Immunogenicity 3.0`)

t_deepNetBim <- quantile(unique_canonical_list$DeepNetBim, probs = 0.99, na.rm = T)
t_deepHLApan <- quantile(unique_canonical_list$DeepHLApan, probs = 0.99, na.rm = T)
t_immunogen <- quantile(unique_canonical_list$`Immunogenicity 3.0`, probs = 0.99, na.rm = T)

t_deepNetBim
t_deepHLApan
t_immunogen

noncanonical_list <- Total_sample[Total_sample$IsCanonical == F, ]
unique_noncanonical_list <- noncanonical_list[!duplicated(noncanonical_list$key), ]
nrow(unique_noncanonical_list)

unique_noncanonical_list$pTSA <- "No"

## intersection
#unique_noncanonical_list[(!is.na(unique_noncanonical_list$DeepNetBim) & 
#                           unique_noncanonical_list$DeepNetBim >= t_deepNetBim) &
#                         (!is.na(unique_noncanonical_list$DeepHLApan) & 
#                            unique_noncanonical_list$DeepHLApan >= t_deepHLApan) &
#                         (!is.na(unique_noncanonical_list$`Immunogenicity 3.0`) & 
#                            unique_noncanonical_list$`Immunogenicity 3.0` >= t_immunogen), ]$pTSA <- "Yes"

## Union
unique_noncanonical_list[!is.na(unique_noncanonical_list$DeepNetBim) & 
                           unique_noncanonical_list$DeepNetBim >= t_deepNetBim, ]$pTSA <- "Yes"
nrow(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes", ])

unique_noncanonical_list[!is.na(unique_noncanonical_list$DeepHLApan) & 
                           unique_noncanonical_list$DeepHLApan >= t_deepHLApan, ]$pTSA <- "Yes"
nrow(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes", ])

unique_noncanonical_list[!is.na(unique_noncanonical_list$`Immunogenicity 3.0`) & 
                          unique_noncanonical_list$`Immunogenicity 3.0` >= t_immunogen, ]$pTSA <- "Yes"
nrow(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes", ])

unique_noncanonical_list[unique_noncanonical_list$`HLA-ligand` == "Yes" | 
                           unique_noncanonical_list$`IEAtlas normal` == "Yes",]$pTSA <- "No"

nrow(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes", ])
nrow(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes" & unique_noncanonical_list$Report == "Known", ])
nrow(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes" & unique_noncanonical_list$Report == "Novel", ])

tmp <- unique_noncanonical_list[!is.na(unique_noncanonical_list$DeepNetBim) & 
                                  unique_noncanonical_list$DeepNetBim >= t_deepNetBim &
                                  !is.na(unique_noncanonical_list$DeepHLApan) &
                                  unique_noncanonical_list$DeepHLApan >= t_deepHLApan & unique_noncanonical_list$pTSA == "Yes", ]

tmp <- unique_noncanonical_list[!is.na(unique_noncanonical_list$DeepNetBim) & 
                                  unique_noncanonical_list$DeepNetBim >= t_deepNetBim &
                                  !is.na(unique_noncanonical_list$`Immunogenicity 3.0`) & 
                                  unique_noncanonical_list$`Immunogenicity 3.0` >= t_immunogen & unique_noncanonical_list$pTSA == "Yes", ]

tmp <- unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes" & unique_noncanonical_list$`IEDB T Cell-` == "Yes" ,]

mdPlot <- ggplot(data=unique_noncanonical_list, aes(x=factor(pTSA), y=NormLog2Reads, fill=pTSA)) +
  theme_bw() +
  scale_fill_manual(values=c(brewer.pal(n = 10, name = "Set3")[5], brewer.pal(n = 10, name = "Set3")[4])) +
  geom_boxplot() +
  theme(text = element_text(size=20)) +
  scale_y_continuous(n.breaks = 6, limits = c(-4, 6)) +
  staticThemeTop +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0.1, 0.1, -0.0, -0.3, "in")) +
  labs(y= "", x="")
mdPlot

nrow(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes", ])
nrow(unique_noncanonical_list[unique_noncanonical_list$pTSA == "No", ])

shapiro.test(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes", ]$NormLog2Reads)
shapiro.test(unique_noncanonical_list[unique_noncanonical_list$pTSA == "No", ]$NormLog2Reads)

wilcox.test(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes", ]$NormLog2Reads,
       unique_noncanonical_list[unique_noncanonical_list$pTSA == "No", ]$NormLog2Reads)

t.test(unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes", ]$NormLog2Reads,
        unique_noncanonical_list[unique_noncanonical_list$pTSA == "No", ]$NormLog2Reads)

ggsave("Figure6_RNA_pTSA.png", plot = mdPlot, width = 4, height = 5, units = "in", dpi = 600)

pTSA <- unique_noncanonical_list[unique_noncanonical_list$pTSA == "Yes", ]
nrow(pTSA)
nrow(pTSA[pTSA$Memo == "Unique", ])

write.table(pTSA, file = "pTSA.tsv", append = FALSE, quote = F, sep = "\t",
            na = "NA", dec = ".", row.names = F,
            col.names = TRUE)
### Event of pTSA ###
noncanonical_peptide <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "pTSA_report_ratio")
noncanonical_peptide$Category <- factor(noncanonical_peptide$Category, levels = c(
  "Coding", "5`-UTR", "FS", "ncRNA", "IR", "IGR", "3`-UTR", "asRNA", "Unknown", "Coding+AS","IR+AS", "IGR+AS", "FS+AS", "5`-UTR+AS", "ncRNA+AS","Multiple"
))

#display.brewer.pal(n = 10, name = 'Set3')
noncanonical_peptide$Report <- factor(x = noncanonical_peptide$Report, levels = c("Novel", "Known"))
event_plot <- ggplot(data=noncanonical_peptide[noncanonical_peptide$`Mutation = 0` != 0, ], 
                     aes(x=Category, y=`Mutation = 0`, fill = Report)) +
  theme_bw() +
  scale_fill_manual(values=c(brewer.pal(n = 10, name = "Set3")[4], brewer.pal(n = 10, name = "Set3")[5])) +
  geom_bar(stat="identity") +
  staticThemeRightTopInner +
  labs(y= "", x = "") +
  scale_y_continuous(n.breaks = 4, limits = c(0, 20)) +
  theme(text = element_text(size=20), strip.background = element_blank(), 
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(0.1, 0.1, 0, -0.2, "in"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
event_plot

ggsave("Figure6_Event_pTSA_Mut0.png", plot = event_plot, width = 5, height = 5, units = "in", dpi = 600)

#### Event proportion of pTSA ######

pTSA_ratio <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "pTSA_ratio_AS")

pTSA_wild <- pTSA_ratio[pTSA_ratio$Wild_ratio != "NA", ]
pTSA_wild$Wild_ratio <- as.numeric(pTSA_wild$Wild_ratio)
pTSA_mutant <- pTSA_ratio[pTSA_ratio$Mutant_ratio != "NA", ]
pTSA_mutant$Mutant_ratio <- as.numeric(pTSA_mutant$Mutant_ratio)

pTSA_wild$Relative_ratio <- pTSA_wild$Wild_ratio / pTSA_mutant[pTSA_mutant$Category == "Coding", ]$Mutant_ratio
pTSA_mutant$Relative_ratio <- pTSA_mutant$Mutant_ratio / pTSA_mutant[pTSA_mutant$Category == "Coding", ]$Mutant_ratio

pTSA_wild$Category <- factor(pTSA_wild$Category, levels = c("Multiple", "IGR", "AS", "IGR+AS", "FS+AS", "Coding+AS",
                                                             "ncRNA", "5`-UTR", "IR", 
                                                            "FS", "3`-UTR", "asRNA", "5`-UTR+AS",
                                                            "IR+AS", "ncRNA+AS", "Unknown"))

pTSA_mutant$Category <- factor(pTSA_mutant$Category, levels = c("Multiple", "AS", "5`-UTR", "Coding", "FS",
                                                                "3`-UTR", "ncRNA", "IR", "IGR","asRNA",
                                                                "5`-UTR+AS", "IR+AS", "IGR+AS"))

event_plot <- ggplot(data=pTSA_wild[pTSA_wild$Wild_ratio != 0, ], 
                     aes(x=Category, y=Wild_ratio)) +
  theme_bw() +
  geom_bar(stat="identity") +
  staticThemeRightTopInner +
  labs(y= "", x = "") +
  scale_y_continuous(n.breaks = 5, limits = c(0, 0.4)) +
  theme(text = element_text(size=20), strip.background = element_blank(), 
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0.1, 0.1, 0, -0.3, "in"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

event_plot

ggsave("Figure6_pTSA_proportion_Mut0.png", plot = event_plot, width = 6, height = 5, units = "in", dpi = 600)

## Fisher's exact test

fisher.test(matrix(c(12, 116, 1, 5), nrow = 2, ncol = 2, byrow = T))
fisher.test(matrix(c(12, 116, 4, 21), nrow = 2, ncol = 2, byrow = T))
fisher.test(matrix(c(12, 116, 5, 28), nrow = 2, ncol = 2, byrow = T))
fisher.test(matrix(c(12, 116, 7, 71), nrow = 2, ncol = 2, byrow = T))
fisher.test(matrix(c(12, 116, 16, 185), nrow = 2, ncol = 2, byrow = T))
fisher.test(matrix(c(12, 116, 2, 46), nrow = 2, ncol = 2, byrow = T))
fisher.test(matrix(c(12, 116, 1, 36), nrow = 2, ncol = 2, byrow = T))


fisher.test(matrix(c(12, 116, 2, 4), nrow = 2, ncol = 2, byrow = T))
fisher.test(matrix(c(12, 116, 2, 5), nrow = 2, ncol = 2, byrow = T))
fisher.test(matrix(c(12, 116, 1, 11), nrow = 2, ncol = 2, byrow = T))



