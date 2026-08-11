## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE, dpi = 400)


## -----------------------------------------------------------------------------
# Load required libraries
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)
theme_set(theme_classic())


## -----------------------------------------------------------------------------
# Load Seurat object from file
seu <- readRDS("scPure2_HB6_UMAP3D.rds") %>% 
  UpdateSeuratObject() 


## -----------------------------------------------------------------------------
# 1 Convert to cell_data_set object
cds <- as.cell_data_set(seu)
coldat = colData(cds)
fData(cds)$gene_short_name <- rownames(fData(cds))
# 2 Cluster cells (using clustering info from seurat's UMAP)
# Use the Seurat object clustering information
# Assign partitions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition
# Assign the cluster info
list_cluster <- seu@active.ident
cds@clusters$UMAP$clusters <- list_cluster
# Assign UMAP coordinate - cell embedding
cds@int_colData@listData$reducedDims$UMAP <- seu@reductions$umap@cell.embeddings


## -----------------------------------------------------------------------------
# 3 Learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE, verbose = F)


## -----------------------------------------------------------------------------
colData(cds)$Cell_Type = coldat$cluster
p12 = (plot_cells(cds,1,3,
           color_cells_by = 'Cell_Type',
           label_cell_groups = F, 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 4) +
  theme(legend.position = "right")) +
  theme(aspect.ratio = 1)

p23 =  (plot_cells(cds,2,3,
           color_cells_by = 'Cell_Type',
           label_cell_groups = F, 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 4) +
  theme(legend.position = "none")) +
  theme(aspect.ratio = 1)

umap_clplot = (p12+p23)


## -----------------------------------------------------------------------------
# 4 Order the cells in pseudotime
get_earliest_principal_node <- function(cds, cluster="Trans"){
  cell_ids <- which(colData(cds)[, "cluster"] == cluster)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))


## ----include=FALSE------------------------------------------------------------
umap_ptplot13 = plot_cells(cds,1,3, 
           color_cells_by = 'pseudotime',
           label_cell_groups = F, show_trajectory_graph = F,
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)  +
    theme(legend.position = "right", aspect.ratio = 1)

umap_ptplot23 =    plot_cells(cds,2,3, 
           color_cells_by = 'pseudotime',
           label_cell_groups = F, show_trajectory_graph = F,
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)+
    theme(legend.position = "none", aspect.ratio = 1)

umap_ptplot = umap_ptplot13 + umap_ptplot23


## ----umap_ptime,  echo=FALSE, fig.height=8------------------------------------
umap_clplot/
umap_ptplot 


## -----------------------------------------------------------------------------
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))


## ----violin_ptime, echo=FALSE, fig.height=7, fig.width=8----------------------
vplot = ggplot(data.pseudo, aes(monocle3_pseudotime, 
                        reorder(cluster, monocle3_pseudotime, median), 
                        fill = cluster)) +
  geom_jitter(size = 0.1, alpha = .25)+
  geom_violin(draw_quantiles = c(0.5), 
              scale = "width") + 
  labs(y = " ", x = "pseudotime", 
       title = "Ordered Violin of Pseudotime by Cluster")+ 
  theme(aspect.ratio = 1)
vplot


## -----------------------------------------------------------------------------
# saveRDS(cds, file = "csd_HB6.rds")

