library(ggplot2)
library(readxl)
library(pheatmap)

setwd("D:/")
fpkm <- read.csv("GEO submission-Sareh/MATR3_KD_mRNA_FPKM_matrix.csv")

fpkm_mat <- as.matrix(fpkm[, c("siLUC_1_SH.SY5Y_mRNA","siLUC_2_SH.SY5Y_mRNA",
                                 "siLUC_3_SH.SY5Y_mRNA","siMATR3_1_SH.SY5Y_mRNA",
                                 "siMATR3_2_SH.SY5Y_mRNA","siMATR3_3_SH.SY5Y_mRNA",
                                 "siLUC_1_U87_mRNA","siLUC_2_U87_mRNA",
                                 "siLUC_3_U87_mRNA","siMATR3_1_U87_mRNA",
                                 "siMATR3_2_U87_mRNA","siMATR3_3_U87_mRNA")])

 log2 transform
fpkm_log <- log2(fpkm_mat + 1)

fpkm_log <- fpkm_log[apply(fpkm_log, 1, var) > 0.5,]

 PCA
pca_result <- prcomp(t(fpkm_log), scale.=TRUE)
var_explained <- round(summary(pca_result)$importance[2,] * 100, 2)

pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Sample = c("siLUC_1_SH","siLUC_2_SH","siLUC_3_SH",
             "siMATR3_1_SH","siMATR3_2_SH","siMATR3_3_SH",
             "siLUC_1_U87","siLUC_2_U87","siLUC_3_U87",
             "siMATR3_1_U87","siMATR3_2_U87","siMATR3_3_U87"),
  Group = c(rep("LUC_SH",3), rep("M3_SH",3),
            rep("LUC_U87",3), rep("M3_U87",3))
)

pca_df$Group <- factor(pca_df$Group,
                        levels=c("LUC_SH","M3_SH","LUC_U87","M3_U87"))

p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Group, label=Sample)) +
  geom_point(size=3.5, alpha=0.95) +
  ggrepel::geom_text_repel(size=2.8, show.legend=FALSE,
                            box.padding=0.4, max.overlaps=20) +
  scale_color_manual(
    values=c("LUC_SH"="#E84040","M3_SH"="#3AB54A",
             "LUC_U87"="#1a6faf","M3_U87"="#CC79A7"),
    labels=c("LUC SH"="siLUC SH-SY5Y","M3_SH"="siMATR3 SH-SY5Y",
             "LUC_U87"="siLUC U87","M3_U87"="siMATR3 U87")
  ) +
  geom_hline(yintercept=0, linetype="dashed", color="grey70", linewidth=0.4) +
  geom_vline(xintercept=0, linetype="dashed", color="grey70", linewidth=0.4) +
  labs(x=paste0("PC1 (",var_explained[1],"%)"),
       y=paste0("PC2 (",var_explained[2],"%)"),
       color="group") +
  theme_classic() +
  theme(
    plot.background=element_rect(fill="white", color=NA),
    panel.background=element_rect(fill="white"),
    axis.title=element_text(size=11, color="black"),
    axis.text=element_text(size=10, color="black"),
    legend.title=element_text(size=10),
    legend.text=element_text(size=9),
    legend.background=element_rect(fill="white", color="grey80"),
    plot.margin=margin(10,20,10,10)
  )

ggsave("PCA_plot.tiff", plot=p_pca, width=14, height=10, units="cm", dpi=300)
