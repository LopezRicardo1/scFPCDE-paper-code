# Load required libraries
library(readr)
library(fda)
library(scFPCDE)
library(parallel)
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggvenn)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(UpSetR)
library(cowplot)

theme_set(theme_cowplot())

# Set the repository root manually if the script is not run
# from within the RStudio project
# setwd("/path/to/scFPCDE-paper-code")

# Load data from the repository data directory
seu <- readRDS("data/scPure2_HB6_UMAP3D.rds") |>
  UpdateSeuratObject()

cds <- readRDS("data/cds_HB6.rds")

# Extract the principal-graph cell projection distance matrix
mst <- t(
  cds@principal_graph_aux$UMAP$pr_graph_cell_proj_dist
)