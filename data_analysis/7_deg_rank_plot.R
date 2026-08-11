## ============================================================
## Functional curves for top 10 genes from each method
##
## For each trajectory:
##   - top 10 scFPC-DE genes ranked by scFPC-DE p-value
##   - top 10 PseudotimeDE genes ranked by PseudotimeDE p-value
##
## Behavior:
##   1) Displays four independent figures in the RStudio Plots pane
##   2) Saves four TIFF figures
##
## Outputs:
##   traj1_scFPCDE_top10_functional_curves.tiff
##   traj1_PseudotimeDE_top10_functional_curves.tiff
##   traj2_scFPCDE_top10_functional_curves.tiff
##   traj2_PseudotimeDE_top10_functional_curves.tiff
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(monocle3)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

## ------------------------------------------------------------
## User settings
## ------------------------------------------------------------

K_curve <- 10

method_cols_curve <- c(
  "scFPC-DE" = "red3",
  "PseudotimeDE" = "blue3"
)

cluster_order_traj1 <- c("Trans", "Naive", "C-mem1")

cluster_cols_traj1 <- c(
  "Trans"  = 2,
  "Naive"  = 3,
  "C-mem1" = 4
)

cluster_order_traj2 <- c("M-mem1", "DN3", "DN2")

cluster_cols_traj2 <- c(
  "M-mem1" = 5,
  "DN3"    = 6,
  "DN2"    = 7
)

## ============================================================
## Keep only genes present in the full FPCA object
## ============================================================

keep_genes_in_fpca <- function(res, gene_vec) {
  
  gene_vec <- unique(as.character(gene_vec))
  
  gene_ids_fit <- colnames(
    res$fpca.res$fpca_result$xt_hat
  )
  
  gene_ids_y <- colnames(
    res$fpca.res$fpca_result$fda_splines$y
  )
  
  if (is.null(gene_ids_fit)) {
    stop("xt_hat has no gene names")
  }
  
  if (is.null(gene_ids_y)) {
    stop("fda_splines$y has no gene names")
  }
  
  keep <- toupper(gene_vec) %in% toupper(gene_ids_fit) &
    toupper(gene_vec) %in% toupper(gene_ids_y)
  
  gene_vec[keep]
}

## ============================================================
## Extract top genes from each method
## ============================================================

get_top_method_genes <- function(res,
                                 method = c("scFPC-DE", "PseudotimeDE"),
                                 top_k = 10) {
  
  method <- match.arg(method)
  
  if (method == "scFPC-DE") {
    
    gene_vec <- res$fpca.res$D_test_result %>%
      as.data.frame() %>%
      arrange(p_value, desc(D_obs)) %>%
      pull(ID)
    
  } else {
    
    gene_vec <- res$res.gauss %>%
      as.data.frame() %>%
      arrange(fix.pv) %>%
      pull(gene)
  }
  
  gene_vec <- keep_genes_in_fpca(
    res = res,
    gene_vec = gene_vec
  )
  
  head(gene_vec, top_k)
}

## ============================================================
## Extract aligned FPCA, expression, cluster, and pseudotime data
## ============================================================

prepare_curve_data <- function(res,
                               gene_vec,
                               cluster_order,
                               cluster_cols) {
  
  fpca_result <- res$fpca.res$fpca_result
  
  xt_hat <- as.matrix(
    fpca_result$xt_hat
  )
  
  y_centered <- as.matrix(
    fpca_result$fda_splines$y
  )
  
  tt <- fpca_result$fda_splines$argvals
  
  if (is.matrix(tt) || is.data.frame(tt)) {
    tt <- tt[, 1]
  }
  
  tt <- as.numeric(tt)
  
  if (length(tt) != nrow(y_centered)) {
    stop("Pseudotime length does not match the FPCA expression matrix")
  }
  
  if (nrow(xt_hat) != nrow(y_centered)) {
    stop("xt_hat and fda_splines$y have different numbers of rows")
  }
  
  ## ----------------------------------------------------------
  ## Match selected genes
  ## ----------------------------------------------------------
  
  idx_fit <- match(
    toupper(gene_vec),
    toupper(colnames(xt_hat))
  )
  
  idx_y <- match(
    toupper(gene_vec),
    toupper(colnames(y_centered))
  )
  
  keep <- !is.na(idx_fit) & !is.na(idx_y)
  
  gene_vec <- gene_vec[keep]
  idx_fit <- idx_fit[keep]
  idx_y <- idx_y[keep]
  
  if (length(gene_vec) == 0) {
    stop("None of the selected genes were found in the FPCA object")
  }
  
  gene_names <- colnames(y_centered)[idx_y]
  
  ## ----------------------------------------------------------
  ## Match cells between FPCA object and cell_data_set
  ## ----------------------------------------------------------
  
  cds_obj <- res$cds_sub2
  
  cell_ids_fpca <- rownames(y_centered)
  
  if (is.null(cell_ids_fpca)) {
    stop("The FPCA matrix has no cell names")
  }
  
  cell_ids_cds <- colnames(cds_obj)
  
  cell_match <- match(
    cell_ids_fpca,
    cell_ids_cds
  )
  
  if (anyNA(cell_match)) {
    stop("Some FPCA cells were not found in cds_sub2")
  }
  
  ## ----------------------------------------------------------
  ## Extract log-expression and raw counts
  ## ----------------------------------------------------------
  
  log_expr <- SingleCellExperiment::logcounts(cds_obj)
  raw_counts <- SummarizedExperiment::assay(cds_obj, "counts")
  
  gene_match_log <- match(
    toupper(gene_names),
    toupper(rownames(log_expr))
  )
  
  gene_match_count <- match(
    toupper(gene_names),
    toupper(rownames(raw_counts))
  )
  
  if (anyNA(gene_match_log)) {
    stop("Some genes were not found in logcounts(cds_sub2)")
  }
  
  if (anyNA(gene_match_count)) {
    stop("Some genes were not found in counts(cds_sub2)")
  }
  
  log_expr_sub <- as.matrix(
    log_expr[
      gene_match_log,
      cell_match,
      drop = FALSE
    ]
  )
  
  count_sub <- as.matrix(
    raw_counts[
      gene_match_count,
      cell_match,
      drop = FALSE
    ]
  )
  
  gene_means <- Matrix::rowMeans(log_expr_sub)
  names(gene_means) <- gene_names
  
  ## ----------------------------------------------------------
  ## Align cluster labels to FPCA cell order
  ## ----------------------------------------------------------
  
  cluster_all <- as.character(
    monocle3::clusters(cds_obj)
  )
  
  names(cluster_all) <- colnames(cds_obj)
  
  cluster_fpca <- cluster_all[cell_ids_fpca]
  
  clu <- factor(
    cluster_fpca,
    levels = cluster_order
  )
  
  if (anyNA(clu)) {
    missing_clusters <- unique(cluster_fpca[is.na(clu)])
    
    stop(
      "Clusters missing from cluster_order: ",
      paste(missing_clusters, collapse = ", ")
    )
  }
  
  missing_cols <- setdiff(
    levels(clu),
    names(cluster_cols)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing colors for clusters: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  ## ----------------------------------------------------------
  ## Order everything by pseudotime
  ## ----------------------------------------------------------
  
  ord <- order(tt)
  
  list(
    genes = gene_names,
    idx_fit = idx_fit,
    idx_y = idx_y,
    xt_hat = xt_hat[ord, , drop = FALSE],
    y_centered = y_centered[ord, , drop = FALSE],
    tt = tt[ord],
    clusters = clu[ord],
    cluster_cols = cluster_cols[levels(clu)],
    gene_means = gene_means,
    counts = count_sub[, ord, drop = FALSE]
  )
}

## ============================================================
## Plot top functional curves for one method and trajectory
## ============================================================

plot_top_method_curves <- function(res,
                                   gene_vec,
                                   method_name,
                                   trajectory_title,
                                   cluster_order,
                                   cluster_cols,
                                   nrow_panel = 2,
                                   ncol_panel = 5,
                                   point_cex = 0.45,
                                   point_alpha = 0.65,
                                   curve_lwd = 2.5,
                                   zero_tick_lwd = 1) {
  
  curve_data <- prepare_curve_data(
    res = res,
    gene_vec = gene_vec,
    cluster_order = cluster_order,
    cluster_cols = cluster_cols
  )
  
  genes_use <- curve_data$genes
  panel_n <- nrow_panel * ncol_panel
  
  if (length(genes_use) > panel_n) {
    genes_use <- genes_use[seq_len(panel_n)]
  }
  
  point_cols <- grDevices::adjustcolor(
    unname(curve_data$cluster_cols),
    alpha.f = point_alpha
  )
  
  names(point_cols) <- names(curve_data$cluster_cols)
  
  cell_cols <- point_cols[
    as.character(curve_data$clusters)
  ]
  
  old_par <- par(c(
    "mar", "oma", "mfrow", "bg", "fg",
    "col.axis", "col.lab", "col.main", "xpd"
  ))
  
  on.exit(par(old_par), add = TRUE)
  
  layout_matrix <- rbind(
    matrix(
      seq_len(panel_n),
      nrow = nrow_panel,
      byrow = TRUE
    ),
    rep(panel_n + 1, ncol_panel)
  )
  
  layout(
    layout_matrix,
    heights = c(rep(1, nrow_panel), 0.20)
  )
  
  par(
    mar = c(3.2, 3.5, 2.5, 0.8),
    oma = c(0, 0, 3.2, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black"
  )
  
  for (j in seq_len(panel_n)) {
    
    if (j <= length(genes_use)) {
      
      gene_name <- genes_use[j]
      
      original_index <- match(
        gene_name,
        curve_data$genes
      )
      
      idx_fit <- curve_data$idx_fit[original_index]
      idx_y <- curve_data$idx_y[original_index]
      
      gene_mean <- curve_data$gene_means[gene_name]
      
      observed_expression <-
        curve_data$y_centered[, idx_y] + gene_mean
      
      fitted_expression <-
        curve_data$xt_hat[, idx_fit] + gene_mean
      
      raw_counts <- curve_data$counts[
        original_index,
      ]
      
      y_range <- range(
        c(observed_expression, fitted_expression),
        finite = TRUE
      )
      
      plot(
        curve_data$tt,
        observed_expression,
        pch = 16,
        cex = point_cex,
        col = cell_cols,
        xlab = "Pseudotime",
        ylab = "Expression",
        main = gene_name,
        col.main = method_cols_curve[method_name],
        font.main = 2,
        ylim = y_range
      )
      
      lines(
        curve_data$tt,
        fitted_expression,
        col = "black",
        lwd = curve_lwd
      )
      
      zero_index <- which(raw_counts == 0)
      
      if (length(zero_index) > 0) {
        rug(
          curve_data$tt[zero_index],
          side = 1,
          ticksize = 0.03,
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
  
  par(mar = c(0, 0, 0, 0))
  plot.new()
  
  legend(
    "center",
    legend = c(
      names(curve_data$cluster_cols),
      "Fitted curve",
      paste0(method_name, " top-ranked gene"),
      "Zero-count ticks"
    ),
    col = c(
      curve_data$cluster_cols,
      "black",
      method_cols_curve[method_name],
      "black"
    ),
    pch = c(
      rep(16, length(curve_data$cluster_cols)),
      NA,
      15,
      124
    ),
    lty = c(
      rep(NA, length(curve_data$cluster_cols)),
      1,
      NA,
      NA
    ),
    lwd = c(
      rep(NA, length(curve_data$cluster_cols)),
      curve_lwd,
      NA,
      NA
    ),
    pt.cex = c(
      rep(1.2, length(curve_data$cluster_cols)),
      NA,
      1.2,
      1.3
    ),
    horiz = TRUE,
    bty = "n",
    xpd = NA
  )
  
  mtext(
    paste0(
      trajectory_title,
      ": top ",
      length(genes_use),
      " ",
      method_name,
      " functional curves"
    ),
    outer = TRUE,
    side = 3,
    line = 0.8,
    font = 2,
    cex = 1.2,
    col = method_cols_curve[method_name]
  )
  
  invisible(
    list(
      genes = genes_use,
      method = method_name,
      trajectory = trajectory_title
    )
  )
}

## ============================================================
## Build top-10 gene vectors
## ============================================================

top10_scFPCDE_traj1 <- get_top_method_genes(
  res = traj_results$traj1,
  method = "scFPC-DE",
  top_k = K_curve
)

top10_PseudotimeDE_traj1 <- get_top_method_genes(
  res = traj_results$traj1,
  method = "PseudotimeDE",
  top_k = K_curve
)

top10_scFPCDE_traj2 <- get_top_method_genes(
  res = traj_results$traj2,
  method = "scFPC-DE",
  top_k = K_curve
)

top10_PseudotimeDE_traj2 <- get_top_method_genes(
  res = traj_results$traj2,
  method = "PseudotimeDE",
  top_k = K_curve
)

print(top10_scFPCDE_traj1)
print(top10_PseudotimeDE_traj1)
print(top10_scFPCDE_traj2)
print(top10_PseudotimeDE_traj2)

## ============================================================
## Plot in RStudio Plots pane
## ============================================================

plot_top_method_curves(
  res = traj_results$traj1,
  gene_vec = top10_scFPCDE_traj1,
  method_name = "scFPC-DE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1
)

plot_top_method_curves(
  res = traj_results$traj1,
  gene_vec = top10_PseudotimeDE_traj1,
  method_name = "PseudotimeDE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1
)

plot_top_method_curves(
  res = traj_results$traj2,
  gene_vec = top10_scFPCDE_traj2,
  method_name = "scFPC-DE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2
)

plot_top_method_curves(
  res = traj_results$traj2,
  gene_vec = top10_PseudotimeDE_traj2,
  method_name = "PseudotimeDE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2
)

## ============================================================
## Save TIFF helper
## ============================================================

save_top_method_curves_tiff <- function(filename,
                                        res,
                                        gene_vec,
                                        method_name,
                                        trajectory_title,
                                        cluster_order,
                                        cluster_cols,
                                        width = 14,
                                        height = 7,
                                        dpi = 300) {
  
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
  
  plot_top_method_curves(
    res = res,
    gene_vec = gene_vec,
    method_name = method_name,
    trajectory_title = trajectory_title,
    cluster_order = cluster_order,
    cluster_cols = cluster_cols,
    nrow_panel = 2,
    ncol_panel = 5,
    point_cex = 0.45,
    point_alpha = 0.65,
    curve_lwd = 2.5
  )
  
  grDevices::dev.off()
}

## ============================================================
## Save four TIFF figures
## ============================================================

save_top_method_curves_tiff(
  filename = "traj1_scFPCDE_top10_functional_curves.tiff",
  res = traj_results$traj1,
  gene_vec = top10_scFPCDE_traj1,
  method_name = "scFPC-DE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1
)

save_top_method_curves_tiff(
  filename = "traj1_PseudotimeDE_top10_functional_curves.tiff",
  res = traj_results$traj1,
  gene_vec = top10_PseudotimeDE_traj1,
  method_name = "PseudotimeDE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1
)

save_top_method_curves_tiff(
  filename = "traj2_scFPCDE_top10_functional_curves.tiff",
  res = traj_results$traj2,
  gene_vec = top10_scFPCDE_traj2,
  method_name = "scFPC-DE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2
)

save_top_method_curves_tiff(
  filename = "traj2_PseudotimeDE_top10_functional_curves.tiff",
  res = traj_results$traj2,
  gene_vec = top10_PseudotimeDE_traj2,
  method_name = "PseudotimeDE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2
)
