# MATR3-depletion-RNAseq

[[DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20744648.svg)](https://doi.org/10.5281/zenodo.20744648)
[[GEO](https://img.shields.io/badge/GEO-GSE336029-blue)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336029)
[[License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Custom R scripts for the integrated bulk and small RNA sequencing analysis of MATR3 depletion in SH-SY5Y (neuroblastoma) and U87 (glioblastoma) cells, accompanying the Data Descriptor:

Integrated bulk and small RNA sequencing datasets of MATR3 depletion in SH-SY5Y and U87 neural-derived cellular models.* Submitted to *Scientific Data*.

Overview

MATR3 is an RNA-binding protein implicated in RNA metabolism, alternative splicing, and neurodegenerative disease. This repository contains the scripts used to generate every main-text figure from raw/processed pipeline outputs (read counts, FPKM matrices, DESeq2/rMATS results, and small RNA differential expression tables).

Raw and processed sequencing data are deposited in NCBI GEO under accession **[GSE336029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336029)**.

Quick Start

```bash
Clone the repository
git clone https://github.com/yourusername/MATR3-depletion-RNAseq.git
cd MATR3-depletion-RNAseq

Run all figures (if you have the required input files)
Rscript scripts/run_all_figures.R
```

Repository Structure

```
MATR3-depletion-RNAseq/
├── scripts/
│   ├── Fig1_workflow.R                   # Experimental workflow diagram
│   ├── Fig2_knockdown_validation.R       # RT-qPCR + Western blot validation panels
│   ├── Fig3_QC.R                         # RNA-sequencing quality control (QC) summary
│   ├── Fig4_PCA_Correlation_combined.R   # PCA + sample-sample Pearson correlation heatmap
│   ├── Fig5_DEG_combined.R               # Volcano plots (SH-SY5Y, U87) + Venn + shared-DEG heatmap
│   ├── Fig6_Splicing_combined.R          # Alternative splicing event bar chart + Venn diagram
│   └── Fig7_miRNA_combined.R             # miRNA volcano plots + let-7c-3p expression bar chart
├── figures/                              # Rendered output (.tiff/.pdf/.png) for each script above
└── README.md
```

## Requirements

All scripts are written in R (≥ 4.5.1) and rely on the following packages:

```r
install.packages(c(
  "ggplot2", "readxl", "dplyr", "tibble", "tidyr",
  "VennDiagram", "pheatmap", "patchwork", "ggplotify",
  "RColorBrewer", "ggforce", "DiagrammeR", "DiagrammeRsvg", "rsvg"
))


## Data Availability

Raw and processed sequencing data are available at NCBI GEO under accession [GSE336029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336029).

Processed data files used as input for these scripts (count matrices, FPKM tables, DESeq2 results, rMATS outputs, and miRNA differential expression tables) are provided as Supplementary Files with the manuscript and are also available for download from the GEO Series record.

## Usage

Each script in `scripts/` is self-contained: it reads its required input file(s) (paths are set near the top of each script — edit `setwd()` and the file names to match your local copy of the processed data tables), performs the analysis/plotting, and saves a publication-ready `.tiff`/`.pdf` (and `.png`/`.svg` for Fig. 1) to the working directory.

To regenerate a figure:

```r
source("scripts/Fig5_DEG_combined.R")
```

**Input files required by each script** (gene count matrices, FPKM tables, DESeq2/rMATS output tables, and the small RNA differential expression tables) are derived from the raw sequencing data deposited in GEO (accession GSE336029) following the processing pipeline described in the Methods section of the manuscript (fastp → HISAT2/Bowtie → featureCounts/miRDeep2 → DESeq2/rMATS).

## Figure-to-Script Mapping

| Figure | Script | Description |
|--------|--------|-------------|
| Fig. 1 | `Fig1_workflow.R` | Overview of the experimental and analytical workflow |
| Fig. 2 | `Fig2_knockdown_validation.R` | RT-qPCR and Western blot validation of MATR3 knockdown |
| Fig. 3 | `Fig3_QC.R` | **RNA-sequencing quality control (QC) summary** (mRNA + miRNA QC metrics, mapping statistics, Kraken2 taxonomic classification) |
| Fig. 4 | `Fig4_PCA_Correlation_combined.R` | PCA and pairwise sample correlation of bulk RNA-seq libraries |
| Fig. 5 | `Fig5_DEG_combined.R` | Differential expression analysis (volcano plots, Venn, shared-DEG heatmap) |
| Fig. 6 | `Fig6_Splicing_combined.R` | Alternative splicing analysis (event counts, shared splicing genes) |
| Fig. 7 | `Fig7_miRNA_combined.R` | MicroRNA differential expression and let-7c-3p validation |

## Code Availability

The version of the code used to generate the figures in the published manuscript is archived on Zenodo: [https://doi.org/10.5281/zenodo.20744648](https://doi.org/10.5281/zenodo.20744648)

## Citation

If you use these scripts or the associated dataset, please cite the Data Descriptor (citation to be added upon publication) and the GEO accession GSE336029.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

Sareh Shahsavar — for questions about the scripts or dataset, please open an issue on this repository.
