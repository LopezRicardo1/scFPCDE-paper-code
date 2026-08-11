## ============================================================
## Full two-trajectory pipeline
##   - traj1: Trans / Naive / C-mem1
##   - traj2: M-mem1 / DN2 / DN3
##   - each trajectory has its own filtering and tuning settings
##   - all results saved in one list: traj_results
## ============================================================
suppressPackageStartupMessages({
  library(readr)
  library(fda)
  library(scFPCDE)
  library(parallel)
  library(PseudotimeDE)
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
})

# Load Seurat object from file
# seu = readRDS("my directory/scPure2_HB6_UMAP3D.rds")%>% 
#   UpdateSeuratObject() 
# cds = readRDS("my directory/cds_HB6.rds")
# mst = cds@principal_graph_aux$UMAP$pr_graph_cell_proj_dist %>% t()

## ------------------------------------------------------------
## Per-trajectory settings
## ------------------------------------------------------------
traj_defs <- list(
  traj1 = list(
    clusters = c("Trans", "Naive", "C-mem1"),
    analysis_tag = "trans_naive_cmem1",
    ptime_global_lo = 0.01,
    ptime_global_hi = 0.99,
    ptime_cluster_lo = 0.01,
    ptime_cluster_hi = 0.99,
    topvarper_tune = 0.10,
    topvarper_run = 0.10,
    topvarper_deg = 1.00,
    L = 3,
    nbasis_range = 30,
    r_pen_range_all = seq(10, 300, len = 21),
    r_pen_range_local = seq(5, 12, 0.21),
    n_perm_main = 500,
    n_perm_deg = 5,
    mc_cores_ptde = 12
  ),
  traj2 = list(
    clusters = c("M-mem1", "DN2", "DN3"),
    analysis_tag = "mmem1_dn2_dn3",
    ptime_global_lo = 0.15,
    ptime_global_hi = 1.00,
    ptime_cluster_lo = 0.005,
    ptime_cluster_hi = 0.995,
    topvarper_tune = 0.10,
    topvarper_run = 0.10,
    topvarper_deg = 1.00,
    L = 3,
    nbasis_range = 30,
    r_pen_range_all = seq(10, 300, len = 21),
    r_pen_range_local = seq(8, 12, 0.21),
    n_perm_main = 500,
    n_perm_deg = 5,
    mc_cores_ptde = 12
  )
)

## ------------------------------------------------------------
## Storage
## ------------------------------------------------------------
traj_results <- list()
set.seed(123)
## ------------------------------------------------------------
## Run both analyses
## ------------------------------------------------------------
for (nm in names(traj_defs)) {
  
  cfg <- traj_defs[[nm]]
  message("Running ", nm, " ...")
  
  ## ----------------------------------------------------------
  ## 1) Subset cds by clusters
  ## ----------------------------------------------------------
  cds_sub <- cds[, colData(cds)$cluster %in% cfg$clusters]
  
  ptime.filter <- cds_sub %>%
    colData() %>%
    as.data.frame() %>%
    mutate(cell_barcode = colnames(cds_sub)) %>%
    filter(
      monocle3_pseudotime > quantile(monocle3_pseudotime, cfg$ptime_global_lo),
      monocle3_pseudotime < quantile(monocle3_pseudotime, cfg$ptime_global_hi)
    ) %>%
    group_by(cluster) %>%
    filter(
      monocle3_pseudotime > quantile(monocle3_pseudotime, cfg$ptime_cluster_lo),
      monocle3_pseudotime < quantile(monocle3_pseudotime, cfg$ptime_cluster_hi)
    ) %>%
    arrange(monocle3_pseudotime) %>%
    pull(cell_barcode)
  
  cds_sub <- cds_sub[, ptime.filter]
  data.pseudo <- as.data.frame(colData(cds_sub))
  
  ## ----------------------------------------------------------
  ## 2) Matching Seurat subset
  ## ----------------------------------------------------------
  seu_sub <- seu[, names(pseudotime(cds_sub))]
  
  seu_sub <- subset(
    seu_sub,
    subset = nFeature_RNA > quantile(nFeature_RNA, 0.01) &
      nFeature_RNA < quantile(nFeature_RNA, 0.99) &
      percent.mt < quantile(percent.mt, 0.99) &
      percent.mt > quantile(percent.mt, 0.01)
  )
  
  seu_sub <- FindVariableFeatures(seu_sub, selection.method = "vst", nfeatures = 4000)
  top10 <- head(VariableFeatures(seu_sub), 10)
  
  ## ----------------------------------------------------------
  ## 3) Build y, tt, z
  ## ----------------------------------------------------------
  topvargenes <- c(VariableFeatures(seu_sub), "ACTB")
  
  y <- logcounts(cds_sub) %>% as.matrix() %>% t
  tt <- pseudotime(cds_sub)
  tt.cluster <- clusters(cds_sub)
  tt.cluster <- droplevels(tt.cluster)
  
  tt.order <- order(tt)
  y <- y[tt.order, , drop = FALSE]
  tt <- tt[tt.order]
  tt.cluster <- tt.cluster[tt.order]
  
  y <- y[, topvargenes, drop = FALSE]
  z <- ifelse(y == 0, 0, 1)
  if(nm=="traj1")
    gene.low.counts <- scFPCDE_filter_genes(y, 0.25)
  if(nm=="traj2")
    gene.low.counts <- scFPCDE_filter_genes(y, 0.5)
  
  y <- y[, gene.low.counts, drop = FALSE]
  z <- ifelse(y == 0, 0, 1)
  
  y <- t(t(y) - colMeans(y))
  
  gene.min.count <- names(which.min(colSums(z)))
  gene.max.count <- names(which.max(colSums(z)))
  
  ## ----------------------------------------------------------
  ## 4) Tune FPCA
  ## ----------------------------------------------------------
  set.seed(123)
  
  fpc.gcv.res.all <- scFPCDE_tune_fpca(
    y, tt,
    L = cfg$L,
    nbasis_range = cfg$nbasis_range,
    r_pen_range = cfg$r_pen_range_all,
    topvarper = 1
  )
  
  fpc.gcv.res <- scFPCDE_tune_fpca(
    y, tt,
    L = cfg$L,
    nbasis_range = cfg$nbasis_range,
    r_pen_range = cfg$r_pen_range_local,
    topvarper = cfg$topvarper_tune
  )
  
  ## ----------------------------------------------------------
  ## 5) Run scFPC-DE
  ## ----------------------------------------------------------
  set.seed(123)
  
  nbasis.star <- fpc.gcv.res$best_nbasis
  rpen.star <- fpc.gcv.res$best_r_pen
  
  fpca.res <- scFPCDE_run(
    y, tt,
    L = cfg$L,
    use_FPC_F = TRUE,
    r_pen = rpen.star,
    nbasis = nbasis.star,
    n_perm = cfg$n_perm_main,
    topvarper = cfg$topvarper_run,fpc_varmax = T
  )
  
  deg.names.p <- fpca.res$D_test_result %>%
    filter(q_value < 0.05) %>%
    arrange(p_value) %>%
    pull(ID)
  
  deg.names <- fpca.res$D_test_result %>%
    filter(q_value < 0.05) %>%
    arrange(desc(D_obs)) %>%
    pull(ID)
  
  all.genes <- fpca.res$D_test_result$ID
  
  fpca.res.deg <- scFPCDE_run(
    y[, deg.names, drop = FALSE], tt,
    L = cfg$L,
    use_FPC_F = TRUE,
    r_pen = rpen.star,
    nbasis = nbasis.star,
    n_perm = cfg$n_perm_deg,
    topvarper = cfg$topvarper_deg
  )
  
  ## ----------------------------------------------------------
  ## 6) PseudotimeDE
  ## ----------------------------------------------------------
  set.seed(123)
  
  cds_sub2 <- cds_sub[colnames(y), names(tt)]
  
  count.matrix <- BiocGenerics::counts(cds_sub2) %>% as.matrix()
  logcounts.matrix <- logcounts(cds_sub2) %>% as.matrix()
  
  sce <- SingleCellExperiment(
    list(
      counts = count.matrix,
      logcounts = logcounts.matrix
    ),
    colData = DataFrame(
      cell_type = clusters(cds_sub2),
      pseudotime = tt
    ),
    rowData = DataFrame(
      gene_short_name = rownames(count.matrix)
    ),
    metadata = list(
      study = cfg$analysis_tag
    )
  )
  
  ori_tbl <- tibble(
    cell = colnames(sce),
    pseudotime = colData(sce)$pseudotime
  )
  
  res.gauss <- runPseudotimeDE(
    gene.vec = rownames(sce),
    ori.tbl = ori_tbl,
    sub.tbl = NULL,
    mc.cores = cfg$mc_cores_ptde,
    mat = sce,
    assay.use = "counts",
    model = "nb"
  )
  
  gauss.pvals <- res.gauss$fix.pv
  gauss.qvals <- p.adjust(res.gauss$fix.pv, "BH")
  
  ptimeDE.gauss.deg <- res.gauss %>%
    filter(gauss.qvals < 0.05) %>%
    arrange(fix.pv) %>%
    pull(gene)
  
  ## ----------------------------------------------------------
  ## 7) Save all objects for this trajectory
  ## ----------------------------------------------------------
  traj_results[[nm]] <- list(
    config = cfg,
    cds_sub = cds_sub,
    cds_sub2 = cds_sub2,
    seu_sub = seu_sub,
    data.pseudo = data.pseudo,
    top10 = top10,
    topvargenes = topvargenes,
    y = y,
    z = z,
    tt = tt,
    tt.cluster = tt.cluster,
    tt.order = tt.order,
    gene.low.counts = gene.low.counts,
    gene.min.count = gene.min.count,
    gene.max.count = gene.max.count,
    fpc.gcv.res.all = fpc.gcv.res.all,
    fpc.gcv.res = fpc.gcv.res,
    nbasis.star = nbasis.star,
    rpen.star = rpen.star,
    fpca.res = fpca.res,
    fpca.res.deg = fpca.res.deg,
    deg.names = deg.names,
    deg.names.p = deg.names.p,
    all.genes = all.genes,
    sce = sce,
    res.gauss = res.gauss,
    gauss.pvals = gauss.pvals,
    gauss.qvals = gauss.qvals,
    ptimeDE.gauss.deg = ptimeDE.gauss.deg,
    n_cells = ncol(cds_sub),
    n_genes = ncol(y)
  )
}

## ------------------------------------------------------------
## 8) Quick summary table
## ------------------------------------------------------------
traj_summary <- data.frame(
  trajectory = names(traj_results),
  clusters = sapply(traj_results, function(x) paste(x$config$clusters, collapse = ", ")),
  n_cells = sapply(traj_results, function(x) x$n_cells),
  n_HVGs = 4000,
  n_filter_genes = sapply(traj_results, function(x) x$n_genes),
  # best_r_pen = sapply(traj_results, function(x) x$rpen.star),
  #  best_nbasis = sapply(traj_results, function(x) x$nbasis.star),
  n_deg_scFPCDE = sapply(traj_results, function(x) length(x$deg.names)),
  n_deg_ptimeDE = sapply(traj_results, function(x) length(x$ptimeDE.gauss.deg)),
  stringsAsFactors = FALSE
)

traj_summary %>% knitr::kable()

## ------------------------------------------------------------
## 9) Convenience aliases
## ------------------------------------------------------------
res_traj1 <- traj_results$traj1
res_traj2 <- traj_results$traj2

################################################################################
# write(
#   traj_results$traj1$fpca.res.deg$D_test_result %>%
#     arrange(desc(D_obs)) %>%
#     pull(ID) %>% head(100),
#   file = "traj1_scFPCDE_deg_names.txt"
# )
# 
# write(
#   traj_results$traj2$fpca.res.deg$D_test_result %>%
#     arrange(desc(D_obs)) %>%
#     pull(ID) %>% head(100),
#   file = "traj2_scFPCDE_deg_names.txt"
# )
# 
# write(
#   traj_results$traj1$res.gauss %>%
#     arrange(fix.pv) %>%
#     pull(gene) %>% head(100),
#   file = "traj1_pseudotimeDE_deg_names.txt"
# )
# 
# write(
#   traj_results$traj2$res.gauss %>%
#     arrange(fix.pv) %>%
#     pull(gene) %>% head(100),
#   file = "traj2_pseudotimeDE_deg_names.txt"
# )

################################################################################
