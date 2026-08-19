## ============================================================
## Gene functional curve plots
##   scFPC-DE     : top K genes by D_obs
##   PseudotimeDE : top K genes by p-value
##
## Behavior:
##   1) Plots each figure in RStudio Plots pane
##   2) Saves each figure as TIFF with LZW compression
##
## Uses FULL FPCA object for both methods:
##   traj_results$traj*$fpca.res
##
## Saved TIFFs:
##   traj1_scfpcde_top10_curves.tiff
##   traj1_pseudotimeDE_top10_curves.tiff
##   traj2_scfpcde_top10_curves.tiff
##   traj2_pseudotimeDE_top10_curves.tiff
##
## TIFF output:
##   width  = 15 inches
##   height = 8 inches
##   layout = 2 rows x 5 columns
##   dpi    = 300
##
## Fixed cluster-color mapping:
##   traj1: Trans = 2, Naive = 3, C-mem1 = 4
##   traj2: M-mem1 = 5, DN3 = 6, DN2 = 7
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

graphics.off()

## ============================================================
## User-defined settings
## ============================================================

K_top <- 10

curve_tiff_width <- 15
curve_tiff_height <- 8
curve_tiff_dpi <- 300

curve_nrow_panel <- 2
curve_ncol_panel <- 5

## trajectory 1: pseudotime progression
cluster_order_traj1 <- c(
  "Trans",
  "Naive",
  "C-mem1"
)

cluster_cols_traj1 <- c(
  "Trans"  = 2,
  "Naive"  = 3,
  "C-mem1" = 4
)

## trajectory 2: pseudotime progression
cluster_order_traj2 <- c(
  "M-mem1",
  "DN3",
  "DN2"
)

cluster_cols_traj2 <- c(
  "M-mem1" = 5,
  "DN3"    = 6,
  "DN2"    = 7
)

## ============================================================
## White plotting defaults
## ============================================================

par(
  bg = "white",
  fg = "black",
  col.axis = "black",
  col.lab = "black",
  col.main = "black",
  col.sub = "black"
)

get_analysis_clusters <- function(cds_obj) {
  cell_data <- as.data.frame(SummarizedExperiment::colData(cds_obj))
  cluster_column <- c("cluster", "cell_type")
  cluster_column <- cluster_column[cluster_column %in% names(cell_data)]
  if (!length(cluster_column)) {
    stop("cds_obj must contain a cluster or cell_type column.", call. = FALSE)
  }

  cluster_values <- as.character(cell_data[[cluster_column[[1L]]]])
  names(cluster_values) <- rownames(cell_data)
  cluster_values
}

## ============================================================
## Helper: keep genes present in full FPCA object
## ============================================================

keep_genes_in_full_fpca <- function(fpca.obj,
                                    gene_vec) {
  
  gene_vec <- as.character(
    gene_vec
  )
  
  gene_ids_xt <- colnames(
    fpca.obj$fpca_result$xt_hat
  )
  
  gene_ids_y <- colnames(
    fpca.obj$fpca_result$fda_splines$y
  )
  
  if (is.null(gene_ids_xt)) {
    stop("xt_hat has no column names")
  }
  
  if (is.null(gene_ids_y)) {
    stop("fda_splines$y has no column names")
  }
  
  keep <- toupper(gene_vec) %in% toupper(gene_ids_xt) &
    toupper(gene_vec) %in% toupper(gene_ids_y)
  
  gene_vec[
    keep
  ]
}

## ============================================================
## Helper: top scFPC-DE genes by D_obs
## ============================================================

get_top_scfpcde_by_Dobs <- function(res,
                                    top_k = K_top,
                                    q_cut = 0.05) {
  
  if (
    is.null(res$fpca.res) ||
    is.null(res$fpca.res$D_test_result)
  ) {
    stop("res$fpca.res$D_test_result was not found")
  }
  
  dres <- as.data.frame(
    res$fpca.res$D_test_result
  )
  
  if (!"ID" %in% names(dres)) {
    stop("D_test_result must contain ID")
  }
  
  if (!"D_obs" %in% names(dres)) {
    stop("D_test_result must contain D_obs")
  }
  
  if (!"q_value" %in% names(dres)) {
    stop("D_test_result must contain q_value")
  }
  
  dres %>%
    dplyr::filter(
      is.finite(q_value),
      q_value < q_cut,
      is.finite(D_obs)
    ) %>%
    dplyr::arrange(
      dplyr::desc(D_obs)
    ) %>%
    dplyr::pull(
      ID
    ) %>%
    unique() %>%
    keep_genes_in_full_fpca(
      fpca.obj = res$fpca.res,
      gene_vec = .
    ) %>%
    head(
      top_k
    )
}

## ============================================================
## Helper: top PseudotimeDE genes by p-value
## ============================================================

get_top_pseudotimeDE_by_pvalue <- function(res,
                                           top_k = K_top) {
  
  if (is.null(res$res.gauss)) {
    stop("res$res.gauss was not found")
  }
  
  pt_tbl <- as.data.frame(
    res$res.gauss
  )
  
  if (!"gene" %in% names(pt_tbl)) {
    stop("res$res.gauss must contain gene")
  }
  
  if (!"fix.pv" %in% names(pt_tbl)) {
    stop("res$res.gauss must contain fix.pv")
  }
  
  pt_tbl %>%
    dplyr::filter(
      is.finite(fix.pv)
    ) %>%
    dplyr::arrange(
      fix.pv
    ) %>%
    dplyr::pull(
      gene
    ) %>%
    unique() %>%
    keep_genes_in_full_fpca(
      fpca.obj = res$fpca.res,
      gene_vec = .
    ) %>%
    head(
      top_k
    )
}

## ============================================================
## Base R panel plot of selected gene curves on absolute scale
## ============================================================

plot_gene_panels_base_absolute <- function(fpca.obj,
                                           label_vec,
                                           cds_obj,
                                           nrow_panel = 2,
                                           ncol_panel = 5,
                                           point_cex = 0.55,
                                           point_alpha = 0.75,
                                           curve_lwd = 2,
                                           zero_tick_lwd = 1,
                                           cluster_cols = NULL,
                                           cluster_order = NULL,
                                           main_title = NULL) {
  
  stopifnot(
    is.list(fpca.obj)
  )
  
  stopifnot(
    !is.null(fpca.obj$fpca_result$xt_hat)
  )
  
  stopifnot(
    !is.null(fpca.obj$fpca_result$fda_splines$y)
  )
  
  stopifnot(
    !is.null(fpca.obj$fpca_result$fda_splines$argvals)
  )
  
  xt_hat <- fpca.obj$fpca_result$xt_hat
  y_ctr  <- fpca.obj$fpca_result$fda_splines$y
  tt_raw <- fpca.obj$fpca_result$fda_splines$argvals[, 1]
  
  gene_ids_xt <- colnames(
    xt_hat
  )
  
  gene_ids_y <- colnames(
    y_ctr
  )
  
  if (is.null(gene_ids_xt)) {
    stop("xt_hat has no column names")
  }
  
  if (is.null(gene_ids_y)) {
    stop("fda_splines$y has no column names")
  }
  
  label_vec <- as.character(
    label_vec
  )
  
  idx_xt <- match(
    toupper(label_vec),
    toupper(gene_ids_xt)
  )
  
  idx_y <- match(
    toupper(label_vec),
    toupper(gene_ids_y)
  )
  
  keep <- !is.na(idx_xt) &
    !is.na(idx_y)
  
  if (!any(keep)) {
    stop("None of the selected genes were found in both xt_hat and fda_splines$y")
  }
  
  label_vec <- label_vec[
    keep
  ]
  
  idx_xt <- idx_xt[
    keep
  ]
  
  idx_y <- idx_y[
    keep
  ]
  
  panel_n <- nrow_panel * ncol_panel
  
  keep_n <- seq_len(
    min(
      panel_n,
      length(label_vec)
    )
  )
  
  label_vec <- label_vec[
    keep_n
  ]
  
  idx_xt <- idx_xt[
    keep_n
  ]
  
  idx_y <- idx_y[
    keep_n
  ]
  
  genes_use <- gene_ids_y[
    idx_y
  ]
  
  ## ----------------------------------------------------------
  ## Cluster information with fixed named color mapping
  ## ----------------------------------------------------------
  
  clu_raw <- get_analysis_clusters(cds_obj)
  
  if (!is.null(cluster_order)) {
    
    clu <- factor(
      clu_raw,
      levels = cluster_order
    )
    
  } else {
    
    clu <- factor(
      clu_raw
    )
  }
  
  if (length(clu) != nrow(y_ctr)) {
    stop("Number of cells in cds_obj clusters does not match number of rows in centered matrix")
  }
  
  if (anyNA(clu)) {
    
    missing_clu <- sort(
      unique(
        clu_raw[
          is.na(clu)
        ]
      )
    )
    
    stop(
      "Some clusters in cds_obj were not included in cluster_order: ",
      paste(
        missing_clu,
        collapse = ", "
      )
    )
  }
  
  clu_levels <- levels(
    clu
  )
  
  if (is.null(cluster_cols)) {
    
    cluster_cols <- setNames(
      seq_along(clu_levels) + 1,
      clu_levels
    )
  }
  
  if (is.null(names(cluster_cols))) {
    
    if (length(cluster_cols) < length(clu_levels)) {
      stop("cluster_cols must have at least one color per cluster level")
    }
    
    cluster_cols <- setNames(
      cluster_cols[
        seq_along(clu_levels)
      ],
      clu_levels
    )
  }
  
  missing_cols <- setdiff(
    clu_levels,
    names(cluster_cols)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      "Missing colors for cluster levels: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }
  
  cluster_cols_use <- cluster_cols[
    clu_levels
  ]
  
  point_cols_use <- grDevices::adjustcolor(
    unname(cluster_cols_use),
    alpha.f = point_alpha
  )
  
  names(point_cols_use) <- clu_levels
  
  cell_cols <- point_cols_use[
    as.character(clu)
  ]
  
  if (anyNA(cell_cols)) {
    stop("Cell colors contain NA values. Check cluster_order and cluster_cols names.")
  }
  
  ## ----------------------------------------------------------
  ## Expression matrix from cds_obj
  ## ----------------------------------------------------------
  
  expr_mat <- SingleCellExperiment::logcounts(cds_obj)
  
  cell_ids_cds <- colnames(
    expr_mat
  )
  
  gene_ids_cds <- rownames(
    expr_mat
  )
  
  cell_ids_fpca <- rownames(
    y_ctr
  )
  
  if (
    is.null(cell_ids_fpca) ||
    is.null(cell_ids_cds) ||
    is.null(gene_ids_cds)
  ) {
    stop("Could not match cells/genes between fpca object and cds_obj")
  }
  
  cell_match <- match(
    cell_ids_fpca,
    cell_ids_cds
  )
  
  if (anyNA(cell_match)) {
    stop("Some FPCA cells were not found in cds_obj expression matrix")
  }
  
  gene_match <- match(
    toupper(genes_use),
    toupper(gene_ids_cds)
  )
  
  if (anyNA(gene_match)) {
    stop("Some selected genes were not found in cds_obj expression matrix")
  }
  
  expr_sub_sparse <- expr_mat[
    gene_match,
    cell_match,
    drop = FALSE
  ]
  
  expr_sub <- as.matrix(
    expr_sub_sparse
  )
  
  gene_means <- Matrix::rowMeans(
    expr_sub_sparse
  )
  
  names(gene_means) <- genes_use
  
  ## ----------------------------------------------------------
  ## Layout
  ## ----------------------------------------------------------
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  layout_mat <- rbind(
    matrix(
      seq_len(panel_n),
      nrow = nrow_panel,
      byrow = TRUE
    ),
    rep(
      panel_n + 1,
      ncol_panel
    )
  )
  
  layout(
    layout_mat,
    heights = c(
      rep(1, nrow_panel),
      0.22
    )
  )
  
  if (!is.null(main_title)) {
    
    par(
      mar = c(3.2, 3.4, 2.3, 0.8),
      oma = c(0, 0, 2.6, 0),
      bg = "white",
      fg = "black",
      col.axis = "black",
      col.lab = "black",
      col.main = "black",
      col.sub = "black"
    )
    
  } else {
    
    par(
      mar = c(3.2, 3.4, 2.3, 0.8),
      oma = c(0, 0, 0, 0),
      bg = "white",
      fg = "black",
      col.axis = "black",
      col.lab = "black",
      col.main = "black",
      col.sub = "black"
    )
  }
  
  ## ----------------------------------------------------------
  ## Plot each gene
  ## ----------------------------------------------------------
  
  for (j in seq_len(panel_n)) {
    
    if (j <= length(genes_use)) {
      
      g_nm <- genes_use[j]
      g_mu <- gene_means[g_nm]
      
      yj <- y_ctr[, idx_y[j]] + g_mu
      fj <- xt_hat[, idx_xt[j]] + g_mu
      
      raw_abs <- expr_sub[j, ]
      
      y_rng <- range(
        c(
          yj,
          fj
        ),
        finite = TRUE
      )
      
      y_pad <- 0.06 * diff(
        y_rng
      )
      
      if (!is.finite(y_pad) || y_pad <= 0) {
        y_pad <- 0.5
      }
      
      y_rng <- y_rng + c(
        -y_pad,
        y_pad
      )
      
      plot(
        tt_raw,
        yj,
        type = "n",
        xlab = "Pseudotime",
        ylab = "Expression",
        main = g_nm,
        ylim = y_rng,
        axes = FALSE,
        ann = FALSE,
        xaxs = "r",
        yaxs = "r"
      )
      
      points(
        tt_raw,
        yj,
        col = cell_cols,
        pch = 16,
        cex = point_cex
      )
      
      lines(
        tt_raw,
        fj,
        col = "black",
        lwd = curve_lwd
      )
      
      axis(
        side = 1,
        cex.axis = 0.75
      )
      
      axis(
        side = 2,
        cex.axis = 0.75
      )
      
      box()
      
      title(
        main = g_nm,
        xlab = "Pseudotime",
        ylab = "Expression",
        cex.main = 1.0,
        cex.lab = 0.85
      )
      
      zero_idx <- which(
        raw_abs == 0
      )
      
      if (length(zero_idx) > 0) {
        
        rug(
          tt_raw[zero_idx],
          side = 1,
          ticksize = 0.035,
          lwd = zero_tick_lwd,
          col = "black"
        )
      }
      
    } else {
      
      plot.new()
    }
  }
  
  ## ----------------------------------------------------------
  ## Legend
  ## ----------------------------------------------------------
  
  par(
    mar = c(0, 0, 0, 0),
    bg = "white"
  )
  
  plot.new()
  
  legend(
    "center",
    legend = c(
      clu_levels,
      "Fitted curve",
      "Zero-count ticks"
    ),
    col = c(
      cluster_cols_use,
      "black",
      "black"
    ),
    pch = c(
      rep(16, length(clu_levels)),
      NA,
      124
    ),
    lty = c(
      rep(NA, length(clu_levels)),
      1,
      NA
    ),
    lwd = c(
      rep(NA, length(clu_levels)),
      curve_lwd,
      NA
    ),
    pt.cex = c(
      rep(1.2, length(clu_levels)),
      NA,
      1.3
    ),
    horiz = TRUE,
    bty = "n",
    xpd = NA,
    text.col = "black",
    cex = 0.95
  )
  
  if (!is.null(main_title)) {
    
    mtext(
      main_title,
      outer = TRUE,
      side = 3,
      line = 0.7,
      font = 2,
      cex = 1.1
    )
  }
  
  invisible(
    list(
      genes = genes_use,
      gene_means = gene_means,
      cluster_levels = clu_levels,
      cluster_cols = cluster_cols_use,
      idx_xt = idx_xt,
      idx_y = idx_y
    )
  )
}

## ============================================================
## Build top K gene vectors
## ============================================================

topK_scfpcde_curves_traj1 <- get_top_scfpcde_by_Dobs(
  res = traj_results$traj1,
  top_k = K_top
)

topK_ptimeDE_curves_traj1 <- get_top_pseudotimeDE_by_pvalue(
  res = traj_results$traj1,
  top_k = K_top
)

topK_scfpcde_curves_traj2 <- get_top_scfpcde_by_Dobs(
  res = traj_results$traj2,
  top_k = K_top
)

topK_ptimeDE_curves_traj2 <- get_top_pseudotimeDE_by_pvalue(
  res = traj_results$traj2,
  top_k = K_top
)

print(
  topK_scfpcde_curves_traj1
)

print(
  topK_ptimeDE_curves_traj1
)

print(
  topK_scfpcde_curves_traj2
)

print(
  topK_ptimeDE_curves_traj2
)

## ============================================================
## Plot in RStudio Plots pane: 2 x 5 panels
## ============================================================

plot_gene_panels_base_absolute(
  fpca.obj = traj_results$traj1$fpca.res,
  label_vec = topK_scfpcde_curves_traj1,
  cds_obj = traj_results$traj1$cds_sub2,
  nrow_panel = curve_nrow_panel,
  ncol_panel = curve_ncol_panel,
  point_cex = 0.65,
  point_alpha = 0.75,
  curve_lwd = 2,
  zero_tick_lwd = 1,
  cluster_cols = cluster_cols_traj1,
  cluster_order = cluster_order_traj1,
  main_title = "Trajectory 1: scFPC-DE top 10 functional curves"
)

plot_gene_panels_base_absolute(
  fpca.obj = traj_results$traj1$fpca.res,
  label_vec = topK_ptimeDE_curves_traj1,
  cds_obj = traj_results$traj1$cds_sub2,
  nrow_panel = curve_nrow_panel,
  ncol_panel = curve_ncol_panel,
  point_cex = 0.65,
  point_alpha = 0.75,
  curve_lwd = 2,
  zero_tick_lwd = 1,
  cluster_cols = cluster_cols_traj1,
  cluster_order = cluster_order_traj1,
  main_title = "Trajectory 1: PseudotimeDE top 10 functional curves"
)

plot_gene_panels_base_absolute(
  fpca.obj = traj_results$traj2$fpca.res,
  label_vec = topK_scfpcde_curves_traj2,
  cds_obj = traj_results$traj2$cds_sub2,
  nrow_panel = curve_nrow_panel,
  ncol_panel = curve_ncol_panel,
  point_cex = 0.65,
  point_alpha = 0.75,
  curve_lwd = 2,
  zero_tick_lwd = 1,
  cluster_cols = cluster_cols_traj2,
  cluster_order = cluster_order_traj2,
  main_title = "Trajectory 2: scFPC-DE top 10 functional curves"
)

plot_gene_panels_base_absolute(
  fpca.obj = traj_results$traj2$fpca.res,
  label_vec = topK_ptimeDE_curves_traj2,
  cds_obj = traj_results$traj2$cds_sub2,
  nrow_panel = curve_nrow_panel,
  ncol_panel = curve_ncol_panel,
  point_cex = 0.65,
  point_alpha = 0.75,
  curve_lwd = 2,
  zero_tick_lwd = 1,
  cluster_cols = cluster_cols_traj2,
  cluster_order = cluster_order_traj2,
  main_title = "Trajectory 2: PseudotimeDE top 10 functional curves"
)

## ============================================================
## Save TIFF helper: 2 x 5 panels, 15 x 8 inches
## ============================================================

save_gene_curve_tiff <- function(filename,
                                 fpca.obj,
                                 label_vec,
                                 cds_obj,
                                 cluster_cols,
                                 cluster_order,
                                 main_title,
                                 width = curve_tiff_width,
                                 height = curve_tiff_height,
                                 dpi = curve_tiff_dpi,
                                 point_cex = 0.55,
                                 point_alpha = 0.75,
                                 curve_lwd = 2,
                                 zero_tick_lwd = 1) {
  
  grDevices::tiff(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    compression = "lzw",
    type = "cairo",
    bg = "white"
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  par(
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    col.sub = "black"
  )
  
  plot_gene_panels_base_absolute(
    fpca.obj = fpca.obj,
    label_vec = label_vec,
    cds_obj = cds_obj,
    nrow_panel = curve_nrow_panel,
    ncol_panel = curve_ncol_panel,
    point_cex = point_cex,
    point_alpha = point_alpha,
    curve_lwd = curve_lwd,
    zero_tick_lwd = zero_tick_lwd,
    cluster_cols = cluster_cols,
    cluster_order = cluster_order,
    main_title = main_title
  )
  
  invisible(
    filename
  )
}

## ============================================================
## Save all four TIFF figures
## ============================================================

save_gene_curve_tiff(
  filename = "traj1_scfpcde_top10_curves.tiff",
  fpca.obj = traj_results$traj1$fpca.res,
  label_vec = topK_scfpcde_curves_traj1,
  cds_obj = traj_results$traj1$cds_sub2,
  cluster_cols = cluster_cols_traj1,
  cluster_order = cluster_order_traj1,
  main_title = "Trajectory 1: scFPC-DE top 10 functional curves",
  width = curve_tiff_width,
  height = curve_tiff_height,
  dpi = curve_tiff_dpi
)

save_gene_curve_tiff(
  filename = "traj1_pseudotimeDE_top10_curves.tiff",
  fpca.obj = traj_results$traj1$fpca.res,
  label_vec = topK_ptimeDE_curves_traj1,
  cds_obj = traj_results$traj1$cds_sub2,
  cluster_cols = cluster_cols_traj1,
  cluster_order = cluster_order_traj1,
  main_title = "Trajectory 1: PseudotimeDE top 10 functional curves",
  width = curve_tiff_width,
  height = curve_tiff_height,
  dpi = curve_tiff_dpi
)

save_gene_curve_tiff(
  filename = "traj2_scfpcde_top10_curves.tiff",
  fpca.obj = traj_results$traj2$fpca.res,
  label_vec = topK_scfpcde_curves_traj2,
  cds_obj = traj_results$traj2$cds_sub2,
  cluster_cols = cluster_cols_traj2,
  cluster_order = cluster_order_traj2,
  main_title = "Trajectory 2: scFPC-DE top 10 functional curves",
  width = curve_tiff_width,
  height = curve_tiff_height,
  dpi = curve_tiff_dpi
)

save_gene_curve_tiff(
  filename = "traj2_pseudotimeDE_top10_curves.tiff",
  fpca.obj = traj_results$traj2$fpca.res,
  label_vec = topK_ptimeDE_curves_traj2,
  cds_obj = traj_results$traj2$cds_sub2,
  cluster_cols = cluster_cols_traj2,
  cluster_order = cluster_order_traj2,
  main_title = "Trajectory 2: PseudotimeDE top 10 functional curves",
  width = curve_tiff_width,
  height = curve_tiff_height,
  dpi = curve_tiff_dpi
)

## ============================================================
## Confirm saved files
## ============================================================

saved_curve_files <- c(
  "traj1_scfpcde_top10_curves.tiff",
  "traj1_pseudotimeDE_top10_curves.tiff",
  "traj2_scfpcde_top10_curves.tiff",
  "traj2_pseudotimeDE_top10_curves.tiff"
)

cat("\nCurrent working directory:\n")
print(
  getwd()
)

cat("\nSaved TIFF files:\n")
print(
  normalizePath(
    saved_curve_files,
    mustWork = FALSE
  )
)

cat("\nFile exists check:\n")
print(
  file.exists(saved_curve_files)
)
