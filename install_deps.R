#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install CRAN packages
cran_packages <- c(
  "dplyr", "tidyr", "RColorBrewer", "ggplot2", "stringr",
  "nnls", "ggpubr", "remotes"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, ask = FALSE)
  }
}

# Install Bioconductor packages
bioc_packages <- c(
  "Seurat", "SeuratObject", "Signac", "GenomicRanges",
  "WGCNA", "ggtree", "ape", "Matrix", "matrixStats",
  "plyranges", "sparseMatrixStats", "easylift"
)

BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

cat("Dependencies installed successfully!\n")
