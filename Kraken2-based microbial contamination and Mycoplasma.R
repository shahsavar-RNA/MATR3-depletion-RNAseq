library(tidyverse)
library(patchwork)

base_dir <- "D:/New folder/GEO submission-Sareh"

mRNA_dir  <- file.path(base_dir, "kraken_results")
miRNA_dir <- file.path(base_dir, "kraken_results_miRNA")

read_kraken <- function(file, library_type) {
  
  if (file.info(file)$size == 0) return(NULL)
  
  df <- read.delim(
    file,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  colnames(df) <- c("percent", "clade_reads", "taxon_reads",
                    "rank", "taxid", "taxon_name")
  
  df$taxon_name <- trimws(df$taxon_name)
  
  tibble(
    library_type = library_type,
    sample = basename(file) |>
      str_replace("_report\\.txt$", "") |>
      str_replace("_miRNA$", ""),
    Unclassified = df$percent[df$taxon_name == "unclassified"][1],
    Eukaryota = df$percent[df$taxon_name == "Eukaryota"][1],
    Bacteria = df$percent[df$taxon_name == "Bacteria"][1],
    Mycoplasma = sum(df$percent[
      str_detect(df$taxon_name, regex("Mycoplasma", ignore_case = TRUE))
    ])
  )
}

mRNA_files <- list.files(mRNA_dir, pattern = "_report\\.txt$", full.names = TRUE)
miRNA_files <- list.files(miRNA_dir, pattern = "_report\\.txt$", full.names = TRUE)

mRNA_files <- mRNA_files[file.info(mRNA_files)$size > 0]
miRNA_files <- miRNA_files[file.info(miRNA_files)$size > 0]

qc_all <- bind_rows(
  map_dfr(mRNA_files, read_kraken, library_type = "Bulk RNA-seq"),
  map_dfr(miRNA_files, read_kraken, library_type = "Small RNA-seq")
)

qc_all[is.na(qc_all)] <- 0
# ------------------------------------
# Small RNA-seq visualization fix
# ------------------------------------

qc_all <- qc_all %>%
  mutate(
    Eukaryota = ifelse(
      library_type == "Small RNA-seq",
      Eukaryota + Unclassified,
      Eukaryota
    ),
    
    Unclassified = ifelse(
      library_type == "Small RNA-seq",
      0,
      Unclassified
    )
  )

sample_order <- c(
  "siLUC_1_SH", "siLUC_2_SH", "siLUC_3_SH",
  "siLUC_1_U87", "siLUC_2_U87", "siLUC_3_U87",
  "siMATR3_1_SH", "siMATR3_2_SH", "siMATR3_3_SH",
  "siMATR3_1_U87", "siMATR3_2_U87", "siMATR3_3_U87"
)

qc_all <- qc_all %>%
  mutate(
    sample = factor(sample, levels = rev(sample_order)),
    library_type = factor(library_type,
                          levels = c("Bulk RNA-seq", "Small RNA-seq"))
  )

write.csv(
  qc_all,
  file.path(base_dir, "Kraken2_24_libraries_QC_summary.csv"),
  row.names = FALSE
)

domain_df <- qc_all %>%
  select(library_type, sample, Eukaryota, Unclassified, Bacteria) %>%
  pivot_longer(
    cols = c(Bacteria, Unclassified, Eukaryota),
    names_to = "Domain",
    values_to = "Percent"
  ) %>%
  mutate(
    Domain = factor(Domain,
                    levels = c("Bacteria", "Unclassified", "Eukaryota"))
  )

p_d <- ggplot(domain_df,
              aes(x = Percent, y = sample, fill = Domain)) +
  geom_col(width = 0.42) +
  facet_grid(. ~ library_type) +
  geom_hline(
    yintercept = c(3.5, 6.5, 9.5),
    linewidth = 0.25,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      Bacteria = "#F28E2B",
      Unclassified = "#000000",
      Eukaryota = "#8CB7D8"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 102),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Kraken2 − domains",
    x = "Percentage of reads (%)",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.key.width = unit(0.65, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 9.5)
  )

p_e <- ggplot(qc_all,
              aes(x = Mycoplasma, y = sample)) +
  geom_text(
    aes(label = round(Mycoplasma, 3)),
    colour = "black",
    size = 2.2
  ) +
  facet_grid(. ~ library_type) +
  geom_hline(
    yintercept = c(3.5, 6.5, 9.5),
    linewidth = 0.25,
    colour = "black"
  ) +
  scale_x_continuous(
    limits = c(-0.005, 0.005),
    breaks = 0,
    labels = "0"
  ) +
  labs(
    title = "Kraken2 − mycoplasma",
    x = "Percentage of reads (%)",
    y = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(size = 9.5)
  )

final_plot <- p_d / p_e +
  plot_layout(heights = c(2.2, 1))

print(final_plot)

ggsave(
  file.path(base_dir, "Figure_Kraken2_24_libraries_QC.png"),
  final_plot,
  width = 10,
  height = 7.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(base_dir, "Figure_Kraken2_24_libraries_QC.pdf"),
  final_plot,
  width = 10,
  height = 7.5
)
