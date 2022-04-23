BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(cluster)

setwd("C:\\Users\\progi\\Desktop\\ntMODa")
data <- read.csv(file = "KEGGPathway.txt", sep="\t", header=TRUE)
row.names(data) <- data$Term
matrix <- data[, 2:3]
matrix <- data.matrix(matrix)
#matrix <- t(matrix)

heatmap <- Heatmap(matrix = matrix, column_title="parent genes", row_title = "sets",
                   name = "expression", cluster_rows = as.dendrogram(diana(matrix)),
                   cluster_columns = as.dendrogram(agnes(t(matrix)))
                   ,column_names_gp = gpar(fontsize = 8)
)

library(circlize)
cm = colorRamp2(c(0, 5, 10), c("white", "yellow", "red4"))
heatmap <- Heatmap(matrix = matrix, name = "expression" ,column_names_gp = gpar(fontsize = 8), width = unit(3, "cm"),
                   rect_gp = gpar(col = "gray14"),
                   show_column_dend = FALSE, show_row_dend = FALSE, heatmap_legend_param = list(
                     title = "-log(p-value)", at = c(0,2,4,6,8,10), labels = c(0,2,4,6,8,10), legend_width = unit(2, "cm"), legend_direction = "horizontal", border = "gray14"), col = cm,
)

draw(heatmap, heatmap_legend_side = "top")

