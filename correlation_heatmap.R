library(pheatmap)
library(RColorBrewer)
library(readxl)

setwd("D:/")
corr_data <- read_excel("CORRELATION.xlsx", col_names = TRUE)
corr_data <- as.data.frame(corr_data)
rownames(corr_data) <- corr_data[, 1]
corr_data <- corr_data[, -1]
corr_data <- as.matrix(corr_data)

annotation <- data.frame(
  CellLine = ifelse(grepl("U87", colnames(corr_data)), "U87", "SH-SY5Y"),
  Condition = ifelse(grepl("siMATR3", colnames(corr_data)), "KD", "Control")
)
rownames(annotation) <- colnames(corr_data)

ann_colors <- list(
  CellLine = c("SH-SY5Y" = "#1a6faf", "U87" = "#f5c800"),
  Condition = c("Control" = "#888888", "KD" = "#E84040")
)

tiff("Correlation_heatmap_with_numbers.tiff", 
     width = 16, 
     height = 14, 
     units = "cm", 
     res = 300,
     compression = "lzw")

pheatmap(corr_data,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(0.6, 1, length.out = 101),
         annotation_col = annotation,
         annotation_row = annotation,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 4.5,
         fontsize = 7,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,
         border_color = "grey60",
         main = "Pearson correlation between samples",
         legend = TRUE,
         legend_breaks = c(0.6, 0.7, 0.8, 0.9, 1.0),
         legend_labels = c("0.6", "0.7", "0.8", "0.9", "1.0"))

dev.off()

print("✅ Correlation heatmap saved as Correlation_heatmap_with_numbers.tiff")
