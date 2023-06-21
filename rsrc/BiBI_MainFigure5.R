setwd("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy")

staticThemeRightTopInner <- theme(text = element_text(size=25, color = "black")
                                  , axis.text.x = element_text(size=20, color="black"), 
                                  axis.text.y = element_text(size=20, color="black"),
                                  legend.justification = c("right", "top"),
                                  legend.position= c(.98, .99),
                                  legend.text = element_text(size=20, color = "black"))

## load noncaonical FDR 5% ##
fdr_5_res_noncanonical <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "Noncanonical")
fdr_5_res_canonical <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "Canonical")

## Peptide proportion
canonical_peptide <- fdr_5_res_canonical
canonical_peptide <- canonical_peptide[!duplicated(canonical_peptide$InferredPeptide), ]
noncanonical_peptide <- fdr_5_res_noncanonical
noncanonical_peptide <- noncanonical_peptide[!duplicated(noncanonical_peptide$InferredPeptide), ]

peptide_matrix <- data.frame(matrix(nrow=0, ncol=12))
colnames(peptide_matrix) <- c("Not reported","Reported", "IEDB", "IEAtlas cancer", "IEAtlas normal", 
                              "HLA-ligand",
                     "Canonical (Scull et al.)", "Canonical (Cuevas et al.)", "Canonical (Laumont et al.)",
                     "Noncanonical (Scull et al.)", "Noncanonical (Cuevas et al.)", "Noncanonical (Laumont et al.)")

peptide_matrix <- data.frame(matrix(nrow=0, ncol=3))
colnames(peptide_matrix) <- c("Class","Label", "Count")
peptide_matrix[nrow(peptide_matrix)+1, ] <- c("cMAP", "Reported", nrow(canonical_peptide[canonical_peptide$NumberOfMappedSources > 0, ]))
peptide_matrix[nrow(peptide_matrix)+1, ] <- c("cMAP", "Not reported", nrow(canonical_peptide[canonical_peptide$NumberOfMappedSources == 0, ]))
peptide_matrix[nrow(peptide_matrix)+1, ] <- c("ncMAP", "Reported", nrow(noncanonical_peptide[noncanonical_peptide$NumberOfMappedSources > 0, ]))
peptide_matrix[nrow(peptide_matrix)+1, ] <- c("ncMAP", "Not reported", nrow(noncanonical_peptide[noncanonical_peptide$NumberOfMappedSources == 0, ]))
peptide_matrix$Count <- as.double(peptide_matrix$Count)

peptide_matrix[peptide_matrix$Class == "cMAP", ]$Count <- peptide_matrix[peptide_matrix$Class == "cMAP", ]$Count / nrow(canonical_peptide)
peptide_matrix[peptide_matrix$Class == "ncMAP", ]$Count <- peptide_matrix[peptide_matrix$Class == "ncMAP", ]$Count / nrow(noncanonical_peptide)
peptide_matrix

 
display.brewer.pal(n = 10, name = 'Set3')
display.brewer.pal(n = 10, name = 'Paired')

tdRatio <- ggplot(data=peptide_matrix, aes(x=Class, y=Count, fill = Label)) +
  theme_bw() +
  scale_fill_manual(values=c(brewer.pal(n = 10, name = "Paired")[8], brewer.pal(n = 10, name = "Paired")[4])) +
  geom_bar(stat="identity") +
  staticThemeRightTopInner +
  labs(y= "", x = "") +
  theme(text = element_text(size=25), strip.background = element_blank(), legend.title = element_blank(), legend.position = "top",
        plot.margin = margin(0.1, 0.1, 0, -0.3, "in"))

tdRatio

ggsave("Figure5_PrevMatch.png", plot = tdRatio, width = 4, height = 5, units = "in", dpi = 600)

### Event type ####

noncanonical_peptide <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "Report_ratio")

noncanonical_peptide$Category <- factor(noncanonical_peptide$Category, levels = c(
  "Coding","FS", "5`-UTR", "ncRNA", "IR", "3`-UTR", "asRNA", "IGR", "Unknown", "Coding+AS","IR+AS", "IGR+AS", "FS+AS", "5`-UTR+AS", "ncRNA+AS",
  "Multiple"
))

noncanonical_peptide$Report = factor(noncanonical_peptide$Report, levels = c("Novel", "Known"))

#display.brewer.pal(n = 10, name = 'Set3')
event_plot <- ggplot(data=noncanonical_peptide[noncanonical_peptide$`Mutation > 0` != 0, ], 
                  aes(x=Category, y=`Mutation > 0`, fill = Report)) +
  theme_bw() +
  scale_fill_manual(values=c(brewer.pal(n = 10, name = "Set3")[4], brewer.pal(n = 10, name = "Set3")[5])) +
  geom_bar(stat="identity") +
  staticThemeRightTopInner +
  labs(y= "", x = "") +
  theme(text = element_text(size=20), strip.background = element_blank(), 
        #legend.position = "none",
        legend.title = element_blank(), 
        plot.margin = margin(0.1, 0.1, 0, -0.3, "in"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

event_plot

ggsave("Figure5_Event_Mut1.png", plot = event_plot, width = 6, height = 5, units = "in", dpi = 600)

#### Categorize mutations

category_mutations <- read_excel(path = "BAAnalysis_FDR.xlsx", sheet = "Mutations_Category")
category_mutations <- category_mutations[category_mutations$Class != "All", ]
category_mutations$Category <- factor(category_mutations$Category, 
                                      levels = c("1 SNV", "2 SNVs", "3 SNVs", "7 SNVs", "1 INS",
                                                 "1 DEL", "1 SNV + 1 DEL", "3 SNVs + 1 INS"))
category_mutations$Type <- factor(category_mutations$Type, levels = c("Full", "Partial", "None"))

tdRatio <- ggplot(data=category_mutations, aes(x=Category, y=Count, fill = Type)) +
  theme_bw() +
  scale_fill_manual(values=c(brewer.pal(n = 10, name = "Paired")[2], 
                             brewer.pal(n = 10, name = "Paired")[1],
                             brewer.pal(n = 10, name = "Paired")[5])) +
  geom_bar(stat="identity") +
  staticThemeRightTopInner +
  labs(y= "", x = "") +
  theme(text = element_text(size=25), strip.background = element_blank(), 
        legend.title = element_blank(), legend.position = "top",
        plot.margin = margin(0.1, 0.1, 0, -0.3, "in")) +
  coord_flip() +
  facet_grid(cols = vars(`Class`), scales = "free") +
  scale_x_discrete(limits = rev(levels(category_mutations$Category)))

tdRatio

ggsave("Figure5_dbSNP+COSMIC.png", plot = tdRatio, width = 6, height = 4, units = "in", dpi = 600)

display.brewer.pal(n = 10, name = 'Set3')
display.brewer.pal(n = 10, name = 'Paired')










