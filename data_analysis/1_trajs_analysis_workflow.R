## ============================================================
## Full two-trajectory pipeline
##   - traj1: Trans / Naive / C-mem1
##   - traj2: M-mem1 / DN2 / DN3
##   - each trajectory has its own filtering and tuning settings
##   - all results saved in one list: traj_results
## ============================================================

## Load either the local full objects or the packaged preprocessed trajectories.
## Sourcing this script directly therefore has the same behavior as running
## 0_load_data.R first.
if (!exists("analysis_data_source", inherits = FALSE)) {
  load_script <- if (file.exists(file.path("data_analysis", "0_load_data.R"))) {
    file.path("data_analysis", "0_load_data.R")
  } else {
    "0_load_data.R"
  }
  source(load_script)
}

if (!analysis_data_source %in% c(
  "local_full_objects",
  "package_preprocessed"
)) {
  stop("Unknown HB6 analysis_data_source.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(scFPCDE)
  library(PseudotimeDE)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(Matrix)
  library(dplyr)
})

if (analysis_data_source == "local_full_objects") {
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(SeuratWrappers)
    library(monocle3)
  })
}

make_packaged_sce <- function(trajectory, cfg) {
  stopifnot(
    identical(dimnames(trajectory$counts), dimnames(trajectory$yt)),
    identical(rownames(trajectory$yt), names(trajectory$tt)),
    identical(rownames(trajectory$yt), names(trajectory$clusters))
  )

  count_matrix <- Matrix::Matrix(t(trajectory$counts), sparse = TRUE)
  logcount_matrix <- Matrix::Matrix(t(trajectory$yt), sparse = TRUE)
  cell_data <- S4Vectors::DataFrame(
    cell_type = trajectory$clusters,
    cluster = trajectory$clusters,
    monocle3_pseudotime = trajectory$tt,
    pseudotime = trajectory$tt,
    row.names = trajectory$cell_id
  )
  gene_data <- S4Vectors::DataFrame(
    gene_short_name = trajectory$gene_id,
    row.names = trajectory$gene_id
  )

  SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = count_matrix,
      logcounts = logcount_matrix
    ),
    colData = cell_data,
    rowData = gene_data,
    metadata = list(
      study = cfg$analysis_tag,
      source = "scFPCDE::scFPCDE_hb6"
    )
  )
}

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
    ncores_scfpcde = 2,
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
    ncores_scfpcde = 2,
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
  
  if (analysis_data_source == "local_full_objects") {
    ## --------------------------------------------------------
    ## 1a) Reproduce cell and gene selection from full objects
    ## --------------------------------------------------------
    cds_sub <- cds[, colData(cds)$cluster %in% cfg$clusters]

    ptime.filter <- cds_sub %>%
      colData() %>%
      as.data.frame() %>%
      mutate(cell_barcode = colnames(cds_sub)) %>%
      filter(
        monocle3_pseudotime > quantile(
          monocle3_pseudotime,
          cfg$ptime_global_lo
        ),
        monocle3_pseudotime < quantile(
          monocle3_pseudotime,
          cfg$ptime_global_hi
        )
      ) %>%
      group_by(cluster) %>%
      filter(
        monocle3_pseudotime > quantile(
          monocle3_pseudotime,
          cfg$ptime_cluster_lo
        ),
        monocle3_pseudotime < quantile(
          monocle3_pseudotime,
          cfg$ptime_cluster_hi
        )
      ) %>%
      arrange(monocle3_pseudotime) %>%
      pull(cell_barcode)

    cds_sub <- cds_sub[, ptime.filter]
    data.pseudo <- as.data.frame(colData(cds_sub))

    seu_sub <- seu[, names(pseudotime(cds_sub))]
    seu_sub <- subset(
      seu_sub,
      subset = nFeature_RNA > quantile(nFeature_RNA, 0.01) &
        nFeature_RNA < quantile(nFeature_RNA, 0.99) &
        percent.mt < quantile(percent.mt, 0.99) &
        percent.mt > quantile(percent.mt, 0.01)
    )
    seu_sub <- FindVariableFeatures(
      seu_sub,
      selection.method = "vst",
      nfeatures = 4000
    )
    top10 <- head(VariableFeatures(seu_sub), 10)
    topvargenes <- c(VariableFeatures(seu_sub), "ACTB")

    y <- logcounts(cds_sub) %>% as.matrix() %>% t
    tt <- pseudotime(cds_sub)
    tt.cluster <- droplevels(clusters(cds_sub))

    tt.order <- order(tt)
    y <- y[tt.order, , drop = FALSE]
    tt <- tt[tt.order]
    tt.cluster <- tt.cluster[tt.order]

    y <- y[, topvargenes, drop = FALSE]
    filter_quantile <- if (nm == "traj1") 0.25 else 0.50
    gene.low.counts <- scFPCDE_filter_genes(y, filter_quantile)
    y <- y[, gene.low.counts, drop = FALSE]
  } else {
    ## --------------------------------------------------------
    ## 1b) Use the exact preprocessed package trajectory
    ## --------------------------------------------------------
    packaged <- hb6_package_data[[nm]]
    if (is.null(packaged)) {
      stop("Packaged trajectory was not found: ", nm, call. = FALSE)
    }

    cds_sub <- make_packaged_sce(packaged, cfg)
    data.pseudo <- as.data.frame(colData(cds_sub))
    seu_sub <- NULL

    y <- packaged$yt
    tt <- packaged$tt
    tt.cluster <- droplevels(packaged$clusters)
    tt.order <- order(tt)
    y <- y[tt.order, , drop = FALSE]
    tt <- tt[tt.order]
    tt.cluster <- tt.cluster[tt.order]

    topvargenes <- colnames(y)
    top10 <- head(topvargenes, 10)
    gene.low.counts <- colnames(y)

    message(
      "Using packaged ", nm, ": ", nrow(y), " cells x ",
      ncol(y), " already-selected genes."
    )
  }

  ## The FPCA model uses centered log-expression in both input modes.
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
    topvarper = 1,
    ncores = cfg$ncores_scfpcde
  )
  
  fpc.gcv.res <- scFPCDE_tune_fpca(
    y, tt,
    L = cfg$L,
    nbasis_range = cfg$nbasis_range,
    r_pen_range = cfg$r_pen_range_local,
    topvarper = cfg$topvarper_tune,
    ncores = cfg$ncores_scfpcde
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
    topvarper = cfg$topvarper_run,
    ncores = cfg$ncores_scfpcde,
    fpc_varmax = TRUE
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
    topvarper = cfg$topvarper_deg,
    ncores = cfg$ncores_scfpcde
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
      cell_type = tt.cluster,
      cluster = tt.cluster,
      pseudotime = tt,
      row.names = names(tt)
    ),
    rowData = DataFrame(
      gene_short_name = rownames(count.matrix),
      row.names = rownames(count.matrix)
    ),
    metadata = list(
      study = cfg$analysis_tag
    )
  )

  ## The package fallback uses this standard SingleCellExperiment for all
  ## downstream count, logcount, and cluster access. The local mode preserves
  ## the original Monocle object for exact reconstruction.
  if (analysis_data_source == "package_preprocessed") {
    cds_sub2 <- sce
  }
  
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
    data_source = analysis_data_source,
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
    n_cells = nrow(y),
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
  n_input_genes = sapply(traj_results, function(x) length(x$topvargenes)),
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
