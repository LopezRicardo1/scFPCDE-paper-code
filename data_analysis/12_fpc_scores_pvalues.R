
## ============================================================
## Individual TIFF plots for manual multipanel assembly
##
## FDR calculation:
##   - scFPC-DE: BH adjustment of D_test_result$p_value
##   - PseudotimeDE: BH adjustment of res.gauss$fix.pv
##   - Stored q_value and gauss.qvals objects are not used
##   - Significant TDEGs are defined as BH-adjusted p-value < 0.05
##
## Individual saved figures:
##
## Score-space plots:
##   traj1_TDEG_categories_FPC_gene_scores.tiff
##   traj2_TDEG_categories_FPC_gene_scores.tiff
##
## Individual FDR histograms:
##   traj1_scFPCDE_FDR_histogram.tiff
##   traj1_PseudotimeDE_FDR_histogram.tiff
##   traj2_scFPCDE_FDR_histogram.tiff
##   traj2_PseudotimeDE_FDR_histogram.tiff
##
## Curve panels:
##   traj1_scfpcde_top10_curves.tiff
##   traj2_scfpcde_top10_curves.tiff
##
## Required objects:
##   traj_results$traj1
##   traj_results$traj2
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

graphics.off()

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
## Global slide and output settings
## ============================================================

slide_width <- 25
slide_height <- 25

wide_panel_width <- slide_width / 2
wide_panel_height <- slide_height / 3

hist_panel_width <- wide_panel_width / 2
hist_panel_height <- wide_panel_height / 2

curve_panel_width <- wide_panel_width * 1.22
curve_panel_height <- wide_panel_height * 0.88

tiff_dpi <- 300

## ============================================================
## Score-space settings
## ============================================================

fdr_cut <- 0.05
ellipse_level <- 0.99

standardize_score_plot <- TRUE

score_xlim <- c(-3, 3)
score_ylim <- c(-3, 3)

score_point_cex <- 1.5

score_main_cex <- 1.35
score_axis_cex <- 1.25
score_label_cex <- 1.40
score_legend_cex <- 1.05
score_legend_pt_cex <- 1.30

score_axis_lwd <- 1.6
ellipse_lwd <- 3

scfpcde_only_col <- "green"
pseudotime_only_col <- "royalblue"
both_col <- "purple"
background_col <- "gray"

null_region_fill <- grDevices::adjustcolor(
  "red3",
  alpha.f = 0.10
)

null_region_border <- "red3"

deg_category_cols <- c(
  "scFPC-DE only"     = scfpcde_only_col,
  "PseudotimeDE only" = pseudotime_only_col,
  "Both"              = both_col,
  "Background"        = background_col
)

## ============================================================
## Histogram settings
## ============================================================

hist_breaks <- seq(
  0,
  1,
  length.out = 41
)

hist_col <- "gray"
hist_border_col <- "black"

hist_main_cex <- 1.05
hist_axis_cex <- 1.00
hist_label_cex <- 1.15
hist_legend_cex <- 0.85
hist_threshold_lwd <- 2.4

hist_ylim_traj1 <- c(0, 900)
hist_ylim_traj2 <- c(0, 650)

## ============================================================
## Curve settings
## ============================================================

K_top <- 10

curve_nrow_panel <- 2
curve_ncol_panel <- 5

curve_point_cex_console <- 0.72
curve_point_cex_tiff <- 0.78
curve_point_alpha <- 0.75

curve_lwd <- 3.3
zero_tick_lwd <- 1.1

gene_title_cex <- 1.65
curve_axis_cex <- 0.78
curve_axis_label_cex <- 0.82
curve_legend_cex <- 1.35
curve_main_title_cex <- 1.25

## ============================================================
## Cluster settings
## ============================================================

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

## ============================================================
## BH adjustment helper
##
## Only finite raw p-values are adjusted.
## Missing or invalid values remain NA.
## ============================================================

bh_adjust_finite <- function(pvalues) {
  
  pvalues <- suppressWarnings(
    as.numeric(pvalues)
  )
  
  qvalues <- rep(
    NA_real_,
    length(pvalues)
  )
  
  valid <- is.finite(pvalues) &
    pvalues >= 0 &
    pvalues <= 1
  
  if (any(valid)) {
    
    qvalues[valid] <- stats::p.adjust(
      pvalues[valid],
      method = "BH"
    )
  }
  
  qvalues
}

## ============================================================
## Standardize FPC score coordinates
## ============================================================

standardize_score_coordinates <- function(score_tbl,
                                          null_ellipse = NULL) {
  
  center <- c(
    FPC1 = mean(score_tbl$FPC1, na.rm = TRUE),
    FPC2 = mean(score_tbl$FPC2, na.rm = TRUE)
  )
  
  scale <- c(
    FPC1 = stats::sd(score_tbl$FPC1, na.rm = TRUE),
    FPC2 = stats::sd(score_tbl$FPC2, na.rm = TRUE)
  )
  
  scale[
    scale == 0 |
      !is.finite(scale)
  ] <- 1
  
  score_tbl_std <- score_tbl
  
  score_tbl_std$FPC1 <- (
    score_tbl$FPC1 - center["FPC1"]
  ) / scale["FPC1"]
  
  score_tbl_std$FPC2 <- (
    score_tbl$FPC2 - center["FPC2"]
  ) / scale["FPC2"]
  
  null_ellipse_std <- null_ellipse
  
  if (!is.null(null_ellipse)) {
    
    null_ellipse_std[, 1] <- (
      null_ellipse[, 1] - center["FPC1"]
    ) / scale["FPC1"]
    
    null_ellipse_std[, 2] <- (
      null_ellipse[, 2] - center["FPC2"]
    ) / scale["FPC2"]
  }
  
  list(
    score_tbl = score_tbl_std,
    null_ellipse = null_ellipse_std,
    center = center,
    scale = scale
  )
}

## ============================================================
## Extract FPCA gene scores
## ============================================================

get_fpca_scores <- function(fpca.res) {
  
  if (is.null(fpca.res$fpca_result)) {
    stop("fpca.res$fpca_result was not found")
  }
  
  scores <- NULL
  
  if (!is.null(fpca.res$fpca_result$scores)) {
    scores <- fpca.res$fpca_result$scores
  }
  
  if (
    is.null(scores) &&
    !is.null(fpca.res$fpca_result$fda_fpca$scores)
  ) {
    scores <- fpca.res$fpca_result$fda_fpca$scores
  }
  
  if (is.null(scores)) {
    stop("No FPCA score matrix was found")
  }
  
  scores <- as.matrix(scores)
  storage.mode(scores) <- "numeric"
  
  if (ncol(scores) < 2) {
    stop("The FPCA score matrix must contain at least two columns")
  }
  
  scores
}

## ============================================================
## Extract gene IDs for FPCA score rows
## ============================================================

get_fpca_gene_ids <- function(res,
                              scores) {
  
  gene_ids <- NULL
  
  if (
    !is.null(rownames(scores)) &&
    length(rownames(scores)) == nrow(scores) &&
    all(nzchar(rownames(scores)))
  ) {
    gene_ids <- rownames(scores)
  }
  
  if (
    is.null(gene_ids) &&
    !is.null(res$fpca.res$D_test_result$ID) &&
    length(res$fpca.res$D_test_result$ID) == nrow(scores)
  ) {
    gene_ids <- res$fpca.res$D_test_result$ID
  }
  
  if (is.null(gene_ids)) {
    stop("Could not identify genes corresponding to FPCA score rows")
  }
  
  toupper(
    trimws(
      as.character(gene_ids)
    )
  )
}

## ============================================================
## Extract scFPC-DE FDR table
##
## BH-adjusted p-values are always recomputed from:
##   D_test_result$p_value
##
## Stored D_test_result$q_value values are ignored.
## ============================================================

get_scfpcde_fdr_table <- function(res) {
  
  if (
    is.null(res$fpca.res) ||
    is.null(res$fpca.res$D_test_result)
  ) {
    stop("scFPC-DE D_test_result was not found")
  }
  
  sc_tbl <- as.data.frame(
    res$fpca.res$D_test_result
  )
  
  required_columns <- c(
    "ID",
    "p_value"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(sc_tbl)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      "D_test_result is missing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  pvalues <- suppressWarnings(
    as.numeric(sc_tbl$p_value)
  )
  
  qvalues <- bh_adjust_finite(
    pvalues
  )
  
  if ("D_obs" %in% names(sc_tbl)) {
    
    D_obs <- suppressWarnings(
      as.numeric(sc_tbl$D_obs)
    )
    
  } else {
    
    D_obs <- rep(
      NA_real_,
      nrow(sc_tbl)
    )
  }
  
  out <- data.frame(
    Gene_ID = toupper(
      trimws(
        as.character(sc_tbl$ID)
      )
    ),
    p_value = pvalues,
    q_value = qvalues,
    D_obs = D_obs,
    stringsAsFactors = FALSE
  )
  
  out <- out[
    !is.na(out$Gene_ID) &
      nzchar(out$Gene_ID),
    ,
    drop = FALSE
  ]
  
  out <- out[
    order(
      out$q_value,
      out$p_value,
      -out$D_obs,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  out <- out[
    !duplicated(out$Gene_ID),
    ,
    drop = FALSE
  ]
  
  rownames(out) <- NULL
  
  out
}

## ============================================================
## Extract PseudotimeDE FDR table
##
## BH-adjusted p-values are always recomputed from:
##   res.gauss$fix.pv
##
## Stored gauss.qvals values are ignored.
## ============================================================

get_pseudotimeDE_fdr_table <- function(res) {
  
  if (is.null(res$res.gauss)) {
    stop("PseudotimeDE res.gauss was not found")
  }
  
  pt_tbl <- as.data.frame(
    res$res.gauss
  )
  
  required_columns <- c(
    "gene",
    "fix.pv"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(pt_tbl)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      "PseudotimeDE result is missing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  pvalues <- suppressWarnings(
    as.numeric(pt_tbl$fix.pv)
  )
  
  qvalues <- bh_adjust_finite(
    pvalues
  )
  
  out <- data.frame(
    Gene_ID = toupper(
      trimws(
        as.character(pt_tbl$gene)
      )
    ),
    p_value = pvalues,
    q_value = qvalues,
    stringsAsFactors = FALSE
  )
  
  out <- out[
    !is.na(out$Gene_ID) &
      nzchar(out$Gene_ID),
    ,
    drop = FALSE
  ]
  
  out <- out[
    order(
      out$q_value,
      out$p_value,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  out <- out[
    !duplicated(out$Gene_ID),
    ,
    drop = FALSE
  ]
  
  rownames(out) <- NULL
  
  out
}

## ============================================================
## Construct consistent BH-based TDEG sets
##
## Use these sets for:
##   - Euler diagrams
##   - overlap counts
##   - enrichment analyses
##   - manuscript totals
## ============================================================

get_bh_tdeg_sets <- function(res,
                             q_cut = 0.05) {
  
  sc_tbl <- get_scfpcde_fdr_table(
    res
  )
  
  pt_tbl <- get_pseudotimeDE_fdr_table(
    res
  )
  
  sc_keep <- is.finite(sc_tbl$q_value) &
    sc_tbl$q_value < q_cut
  
  pt_keep <- is.finite(pt_tbl$q_value) &
    pt_tbl$q_value < q_cut
  
  list(
    scFPCDE = unique(
      sc_tbl$Gene_ID[sc_keep]
    ),
    PseudotimeDE = unique(
      pt_tbl$Gene_ID[pt_keep]
    )
  )
}

## ============================================================
## Count complete BH-based TDEG results
## ============================================================

count_bh_tdegs <- function(res,
                           q_cut = 0.05) {
  
  sc_tbl <- get_scfpcde_fdr_table(
    res
  )
  
  pt_tbl <- get_pseudotimeDE_fdr_table(
    res
  )
  
  c(
    scFPCDE_total_tested = sum(
      is.finite(sc_tbl$q_value)
    ),
    scFPCDE_significant = sum(
      is.finite(sc_tbl$q_value) &
        sc_tbl$q_value < q_cut
    ),
    PseudotimeDE_total_tested = sum(
      is.finite(pt_tbl$q_value)
    ),
    PseudotimeDE_significant = sum(
      is.finite(pt_tbl$q_value) &
        pt_tbl$q_value < q_cut
    )
  )
}

## ============================================================
## Summarize TDEG overlap
## ============================================================

summarize_tdeg_overlap <- function(res,
                                   q_cut = 0.05) {
  
  sets <- get_bh_tdeg_sets(
    res = res,
    q_cut = q_cut
  )
  
  shared <- intersect(
    sets$scFPCDE,
    sets$PseudotimeDE
  )
  
  sc_only <- setdiff(
    sets$scFPCDE,
    sets$PseudotimeDE
  )
  
  pt_only <- setdiff(
    sets$PseudotimeDE,
    sets$scFPCDE
  )
  
  c(
    scFPCDE_only = length(sc_only),
    shared = length(shared),
    PseudotimeDE_only = length(pt_only),
    scFPCDE_total = length(sets$scFPCDE),
    PseudotimeDE_total = length(sets$PseudotimeDE)
  )
}

## ============================================================
## Prepare score/FDR data
## ============================================================

prepare_score_fdr_data <- function(res) {
  
  scores <- get_fpca_scores(
    res$fpca.res
  )
  
  gene_ids <- get_fpca_gene_ids(
    res = res,
    scores = scores
  )
  
  sc_tbl <- get_scfpcde_fdr_table(
    res = res
  )
  
  pt_tbl <- get_pseudotimeDE_fdr_table(
    res = res
  )
  
  sc_match <- match(
    gene_ids,
    sc_tbl$Gene_ID
  )
  
  pt_match <- match(
    gene_ids,
    pt_tbl$Gene_ID
  )
  
  score_tbl <- data.frame(
    Gene_ID = gene_ids,
    FPC1 = scores[, 1],
    FPC2 = scores[, 2],
    scFPCDE_qvalue = sc_tbl$q_value[sc_match],
    PseudotimeDE_qvalue = pt_tbl$q_value[pt_match],
    stringsAsFactors = FALSE
  )
  
  score_tbl <- score_tbl[
    is.finite(score_tbl$FPC1) &
      is.finite(score_tbl$FPC2),
    ,
    drop = FALSE
  ]
  
  rownames(score_tbl) <- NULL
  
  list(
    score_table = score_tbl,
    scFPCDE_table = sc_tbl,
    PseudotimeDE_table = pt_tbl
  )
}

## ============================================================
## scFPC-DE non-DEG boundary from raw FPC scores
## ============================================================

get_scfpcde_null_ellipse <- function(score_tbl,
                                     q_cut = 0.05,
                                     level = 0.99,
                                     n_points = 400) {
  
  is_null <- is.finite(
    score_tbl$scFPCDE_qvalue
  ) &
    score_tbl$scFPCDE_qvalue >= q_cut
  
  X0 <- as.matrix(
    score_tbl[
      is_null,
      c(
        "FPC1",
        "FPC2"
      ),
      drop = FALSE
    ]
  )
  
  X0 <- X0[
    complete.cases(X0),
    ,
    drop = FALSE
  ]
  
  if (nrow(X0) < 5) {
    
    warning(
      "Too few nonsignificant genes to estimate the null ellipse"
    )
    
    return(NULL)
  }
  
  mu <- colMeans(
    X0
  )
  
  S <- stats::cov(
    X0
  )
  
  if (
    any(!is.finite(S)) ||
    nrow(S) != 2 ||
    ncol(S) != 2
  ) {
    
    warning(
      "The null covariance matrix is invalid"
    )
    
    return(NULL)
  }
  
  trace_S <- sum(
    diag(S)
  )
  
  if (
    !is.finite(trace_S) ||
    trace_S <= 0
  ) {
    trace_S <- 1
  }
  
  S <- S +
    diag(
      1e-4 * trace_S,
      nrow = 2
    )
  
  eig <- eigen(
    S,
    symmetric = TRUE
  )
  
  eig$values <- pmax(
    eig$values,
    0
  )
  
  A <- eig$vectors %*%
    diag(
      sqrt(eig$values),
      nrow = 2
    )
  
  theta <- seq(
    0,
    2 * pi,
    length.out = n_points
  )
  
  circle <- rbind(
    cos(theta),
    sin(theta)
  )
  
  radius <- sqrt(
    stats::qchisq(
      level,
      df = 2
    )
  )
  
  ellipse <- t(
    matrix(
      mu,
      nrow = 2,
      ncol = n_points
    ) +
      radius * (A %*% circle)
  )
  
  colnames(ellipse) <- c(
    "FPC1",
    "FPC2"
  )
  
  ellipse
}

## ============================================================
## Prepare full score plotting object
## ============================================================

prepare_tdeg_score_plot_object <- function(res,
                                           q_cut = 0.05,
                                           ellipse_level = 0.99) {
  
  plot_data <- prepare_score_fdr_data(
    res = res
  )
  
  score_tbl <- plot_data$score_table
  
  null_ellipse <- get_scfpcde_null_ellipse(
    score_tbl = score_tbl,
    q_cut = q_cut,
    level = ellipse_level
  )
  
  list(
    plot_data = plot_data,
    score_tbl = score_tbl,
    null_ellipse = null_ellipse
  )
}

## ============================================================
## TDEG category FPC score-space plot
## ============================================================

plot_tdeg_category_fpc_score <- function(plot_object,
                                         trajectory_title,
                                         q_cut = 0.05,
                                         point_cex = score_point_cex,
                                         show_null_fill = TRUE,
                                         standardize_plot = TRUE,
                                         xlim = score_xlim,
                                         ylim = score_ylim,
                                         legend_pos = "topright") {
  
  score_tbl_raw <- plot_object$score_tbl
  null_ellipse_raw <- plot_object$null_ellipse
  
  is_scfpcde <- is.finite(
    score_tbl_raw$scFPCDE_qvalue
  ) &
    score_tbl_raw$scFPCDE_qvalue < q_cut
  
  is_pseudotime <- is.finite(
    score_tbl_raw$PseudotimeDE_qvalue
  ) &
    score_tbl_raw$PseudotimeDE_qvalue < q_cut
  
  deg_category <- ifelse(
    is_scfpcde & is_pseudotime,
    "Both",
    ifelse(
      is_scfpcde & !is_pseudotime,
      "scFPC-DE only",
      ifelse(
        !is_scfpcde & is_pseudotime,
        "PseudotimeDE only",
        "Background"
      )
    )
  )
  
  deg_category <- factor(
    deg_category,
    levels = c(
      "scFPC-DE only",
      "PseudotimeDE only",
      "Both",
      "Background"
    )
  )
  
  point_cols <- deg_category_cols[
    as.character(deg_category)
  ]
  
  if (standardize_plot) {
    
    std_out <- standardize_score_coordinates(
      score_tbl = score_tbl_raw,
      null_ellipse = null_ellipse_raw
    )
    
    score_tbl <- std_out$score_tbl
    null_ellipse <- std_out$null_ellipse
    
    xlab_use <- "Standardized gene FPC1 score"
    ylab_use <- "Standardized gene FPC2 score"
    
  } else {
    
    score_tbl <- score_tbl_raw
    null_ellipse <- null_ellipse_raw
    
    xlim <- range(
      score_tbl$FPC1,
      finite = TRUE
    )
    
    ylim <- range(
      score_tbl$FPC2,
      finite = TRUE
    )
    
    xlab_use <- "Gene FPC1 score"
    ylab_use <- "Gene FPC2 score"
  }
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  par(
    mar = c(5.0, 5.8, 3.2, 1.6),
    mgp = c(3.0, 0.9, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    xaxs = "i",
    yaxs = "i",
    cex.axis = score_axis_cex,
    cex.lab = score_label_cex,
    cex.main = score_main_cex,
    lwd = score_axis_lwd,
    font.main = 2
  )
  
  plot(
    score_tbl$FPC1,
    score_tbl$FPC2,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = xlab_use,
    ylab = ylab_use,
    main = paste0(
      trajectory_title,
      ": TDEG categories in FPC score space"
    )
  )
  
  if (
    !is.null(null_ellipse) &&
    show_null_fill
  ) {
    
    polygon(
      null_ellipse[, 1],
      null_ellipse[, 2],
      col = null_region_fill,
      border = NA
    )
  }
  
  idx_bg <- deg_category == "Background"
  
  points(
    score_tbl$FPC1[idx_bg],
    score_tbl$FPC2[idx_bg],
    pch = 16,
    cex = point_cex * 0.85,
    col = point_cols[idx_bg]
  )
  
  idx_fg <- deg_category != "Background"
  
  points(
    score_tbl$FPC1[idx_fg],
    score_tbl$FPC2[idx_fg],
    pch = 16,
    cex = point_cex,
    col = point_cols[idx_fg]
  )
  
  abline(
    h = 0,
    v = 0,
    lty = 2,
    lwd = 1.4,
    col = "gray55"
  )
  
  if (!is.null(null_ellipse)) {
    
    lines(
      null_ellipse[, 1],
      null_ellipse[, 2],
      col = null_region_border,
      lty = 2,
      lwd = ellipse_lwd
    )
  }
  
  legend(
    legend_pos,
    legend = c(
      "scFPC-DE only",
      "PseudotimeDE only",
      "Both",
      "Background",
      "scFPC-DE non-DEG boundary"
    ),
    col = c(
      scfpcde_only_col,
      pseudotime_only_col,
      both_col,
      background_col,
      null_region_border
    ),
    pch = c(
      16,
      16,
      16,
      16,
      NA
    ),
    lty = c(
      NA,
      NA,
      NA,
      NA,
      2
    ),
    lwd = c(
      NA,
      NA,
      NA,
      NA,
      ellipse_lwd
    ),
    pt.cex = score_legend_pt_cex,
    bty = "n",
    cex = score_legend_cex
  )
  
  invisible(
    list(
      score_tbl_raw = score_tbl_raw,
      score_tbl_plot = score_tbl,
      null_ellipse_raw = null_ellipse_raw,
      null_ellipse_plot = null_ellipse,
      deg_category = deg_category,
      colors = point_cols,
      is_scfpcde = is_scfpcde,
      is_pseudotime = is_pseudotime
    )
  )
}

## ============================================================
## Individual FDR histogram
## ============================================================

plot_fdr_histogram_single <- function(qvalues,
                                      method_label,
                                      trajectory_label,
                                      q_cut = 0.05,
                                      breaks = hist_breaks,
                                      ylim = NULL,
                                      add_legend = TRUE) {
  
  qvalues <- suppressWarnings(
    as.numeric(qvalues)
  )
  
  qvalues <- qvalues[
    is.finite(qvalues) &
      qvalues >= 0 &
      qvalues <= 1
  ]
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  par(
    mar = c(5.0, 5.5, 3.2, 1.4),
    mgp = c(3.0, 0.9, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    xaxs = "i",
    yaxs = "i",
    cex.axis = hist_axis_cex,
    cex.lab = hist_label_cex,
    cex.main = hist_main_cex,
    font.main = 2
  )
  
  if (length(qvalues) == 0) {
    
    plot.new()
    
    title(
      main = paste0(
        trajectory_label,
        ": ",
        method_label,
        " BH-adjusted p-values"
      )
    )
    
    text(
      x = 0.5,
      y = 0.5,
      labels = "No valid adjusted p-values",
      cex = 1.0
    )
    
    return(
      invisible(NULL)
    )
  }
  
  hist_obj <- hist(
    qvalues,
    breaks = breaks,
    col = hist_col,
    border = hist_border_col,
    main = paste0(
      trajectory_label,
      ": ",
      method_label,
      " BH-adjusted p-values"
    ),
    xlab = "BH-adjusted p-value",
    ylab = "Count",
    xlim = c(0, 1),
    ylim = ylim,
    xaxt = "n"
  )
  
  axis(
    side = 1,
    at = seq(
      0,
      1,
      by = 0.25
    ),
    labels = sprintf(
      "%.2f",
      seq(
        0,
        1,
        by = 0.25
      )
    ),
    cex.axis = hist_axis_cex
  )
  
  abline(
    v = q_cut,
    col = "red3",
    lty = 2,
    lwd = hist_threshold_lwd
  )
  
  n_significant <- sum(
    qvalues < q_cut
  )
  
  if (add_legend) {
    
    legend(
      "topright",
      legend = c(
        paste0(
          "FDR = ",
          q_cut
        ),
        paste0(
          "Sig. = ",
          n_significant
        ),
        paste0(
          "Genes = ",
          length(qvalues)
        )
      ),
      col = c(
        "red3",
        NA,
        NA
      ),
      lty = c(
        2,
        NA,
        NA
      ),
      lwd = c(
        hist_threshold_lwd,
        NA,
        NA
      ),
      bty = "n",
      cex = hist_legend_cex
    )
  }
  
  invisible(
    list(
      histogram = hist_obj,
      qvalues = qvalues,
      significant_count = n_significant
    )
  )
}

## ============================================================
## Keep genes present in full FPCA object
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
## Top scFPC-DE genes by D_obs
##
## Significance is based on BH adjustment of raw p-values.
## Significant genes are then ranked by decreasing D_obs.
## ============================================================

get_top_scfpcde_by_Dobs <- function(res,
                                    top_k = K_top,
                                    q_cut = 0.05) {
  
  sc_tbl <- get_scfpcde_fdr_table(
    res = res
  )
  
  selected_genes <- sc_tbl %>%
    dplyr::filter(
      is.finite(q_value),
      q_value < q_cut,
      is.finite(D_obs)
    ) %>%
    dplyr::arrange(
      dplyr::desc(D_obs),
      p_value
    ) %>%
    dplyr::pull(
      Gene_ID
    ) %>%
    unique()
  
  selected_genes <- keep_genes_in_full_fpca(
    fpca.obj = res$fpca.res,
    gene_vec = selected_genes
  )
  
  head(
    selected_genes,
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
                                           point_cex = 0.78,
                                           point_alpha = 0.75,
                                           curve_lwd = 3.3,
                                           zero_tick_lwd = 1.1,
                                           cluster_cols = NULL,
                                           cluster_order = NULL,
                                           main_title = NULL,
                                           gene_title_cex = 1.65,
                                           axis_cex = 0.78,
                                           axis_label_cex = 0.82,
                                           legend_cex = 1.35,
                                           main_title_cex = 1.25) {
  
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
  y_ctr <- fpca.obj$fpca_result$fda_splines$y
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
  
  label_vec <- label_vec[keep]
  idx_xt <- idx_xt[keep]
  idx_y <- idx_y[keep]
  
  panel_n <- nrow_panel * ncol_panel
  
  keep_n <- seq_len(
    min(
      panel_n,
      length(label_vec)
    )
  )
  
  label_vec <- label_vec[keep_n]
  idx_xt <- idx_xt[keep_n]
  idx_y <- idx_y[keep_n]
  
  genes_use <- gene_ids_y[
    idx_y
  ]
  
  ## ----------------------------------------------------------
  ## Cluster information
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
  
  outer_margin <- if (!is.null(main_title)) {
    c(1, 0, 1.7, 0)
  } else {
    c(0, 0, 0, 0)
  }
  
  par(
    mar = c(3.9, 3.2, 2.35, 0.8),
    oma = outer_margin,
    mgp = c(2.45, 0.75, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    col.sub = "black",
    pty = "m"
  )
  
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
      
      if (
        !is.finite(y_pad) ||
        y_pad <= 0
      ) {
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
        main = "",
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
        cex.axis = axis_cex,
        padj = 0.45
      )
      
      axis(
        side = 2,
        cex.axis = axis_cex,
        las = 1
      )
      
      box()
      
      title(
        main = g_nm,
        xlab = "Pseudotime",
        ylab = "Expression",
        cex.main = gene_title_cex,
        cex.lab = axis_label_cex,
        line = 1.0
      )
      
      zero_idx <- which(
        raw_abs == 0
      )
      
      if (length(zero_idx) > 0) {
        
        rug(
          tt_raw[zero_idx],
          side = 1,
          ticksize = 0.032,
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
      rep(1.45, length(clu_levels)),
      NA,
      1.45
    ),
    horiz = TRUE,
    bty = "n",
    xpd = NA,
    text.col = "black",
    cex = legend_cex
  )
  
  if (!is.null(main_title)) {
    
    mtext(
      main_title,
      outer = TRUE,
      side = 3,
      line = 0.25,
      font = 2,
      cex = main_title_cex
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
## Save helpers
## ============================================================

save_score_tiff <- function(filename,
                            plot_object,
                            trajectory_title,
                            legend_pos = "topright") {
  
  grDevices::tiff(
    filename = filename,
    width = wide_panel_width,
    height = wide_panel_height,
    units = "in",
    res = tiff_dpi,
    compression = "lzw",
    type = "cairo",
    bg = "white"
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  plot_tdeg_category_fpc_score(
    plot_object = plot_object,
    trajectory_title = trajectory_title,
    q_cut = fdr_cut,
    point_cex = score_point_cex,
    show_null_fill = TRUE,
    standardize_plot = standardize_score_plot,
    xlim = score_xlim,
    ylim = score_ylim,
    legend_pos = legend_pos
  )
  
  invisible(filename)
}

save_hist_tiff <- function(filename,
                           qvalues,
                           method_label,
                           trajectory_label,
                           ylim) {
  
  grDevices::tiff(
    filename = filename,
    width = hist_panel_width,
    height = hist_panel_height,
    units = "in",
    res = tiff_dpi,
    compression = "lzw",
    type = "cairo",
    bg = "white"
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  plot_fdr_histogram_single(
    qvalues = qvalues,
    method_label = method_label,
    trajectory_label = trajectory_label,
    q_cut = fdr_cut,
    breaks = hist_breaks,
    ylim = ylim,
    add_legend = TRUE
  )
  
  invisible(filename)
}

save_gene_curve_tiff <- function(filename,
                                 fpca.obj,
                                 label_vec,
                                 cds_obj,
                                 cluster_cols,
                                 cluster_order,
                                 main_title) {
  
  grDevices::tiff(
    filename = filename,
    width = curve_panel_width,
    height = curve_panel_height,
    units = "in",
    res = tiff_dpi,
    compression = "lzw",
    type = "cairo",
    bg = "white"
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  plot_gene_panels_base_absolute(
    fpca.obj = fpca.obj,
    label_vec = label_vec,
    cds_obj = cds_obj,
    nrow_panel = curve_nrow_panel,
    ncol_panel = curve_ncol_panel,
    point_cex = curve_point_cex_tiff,
    point_alpha = curve_point_alpha,
    curve_lwd = curve_lwd,
    zero_tick_lwd = zero_tick_lwd,
    cluster_cols = cluster_cols,
    cluster_order = cluster_order,
    main_title = main_title,
    gene_title_cex = gene_title_cex,
    axis_cex = curve_axis_cex,
    axis_label_cex = curve_axis_label_cex,
    legend_cex = curve_legend_cex,
    main_title_cex = curve_main_title_cex
  )
  
  invisible(filename)
}

## ============================================================
## Prepare BH-adjusted result tables
## ============================================================

scfpcde_table_traj1 <- get_scfpcde_fdr_table(
  traj_results$traj1
)

pseudotimeDE_table_traj1 <- get_pseudotimeDE_fdr_table(
  traj_results$traj1
)

scfpcde_table_traj2 <- get_scfpcde_fdr_table(
  traj_results$traj2
)

pseudotimeDE_table_traj2 <- get_pseudotimeDE_fdr_table(
  traj_results$traj2
)

## ============================================================
## Prepare TDEG sets for Euler plots and enrichment
## ============================================================

tdeg_sets_traj1 <- get_bh_tdeg_sets(
  res = traj_results$traj1,
  q_cut = fdr_cut
)

tdeg_sets_traj2 <- get_bh_tdeg_sets(
  res = traj_results$traj2,
  q_cut = fdr_cut
)

## These are the sets that should be supplied to eulerr:
##
## tdeg_sets_traj1$scFPCDE
## tdeg_sets_traj1$PseudotimeDE
## tdeg_sets_traj2$scFPCDE
## tdeg_sets_traj2$PseudotimeDE

## ============================================================
## BH-adjusted TDEG counts and overlaps
## ============================================================

bh_count_table <- rbind(
  Trajectory_1 = count_bh_tdegs(
    traj_results$traj1,
    q_cut = fdr_cut
  ),
  Trajectory_2 = count_bh_tdegs(
    traj_results$traj2,
    q_cut = fdr_cut
  )
)

bh_overlap_table <- rbind(
  Trajectory_1 = summarize_tdeg_overlap(
    traj_results$traj1,
    q_cut = fdr_cut
  ),
  Trajectory_2 = summarize_tdeg_overlap(
    traj_results$traj2,
    q_cut = fdr_cut
  )
)

cat("\nBH-adjusted TDEG counts:\n")

print(
  bh_count_table
)

cat("\nBH-adjusted TDEG overlap counts:\n")

print(
  bh_overlap_table
)

## ============================================================
## Prepare score-space objects
## ============================================================

tdeg_score_traj1 <- prepare_tdeg_score_plot_object(
  res = traj_results$traj1,
  q_cut = fdr_cut,
  ellipse_level = ellipse_level
)

tdeg_score_traj2 <- prepare_tdeg_score_plot_object(
  res = traj_results$traj2,
  q_cut = fdr_cut,
  ellipse_level = ellipse_level
)

## ============================================================
## Select top curve genes
## ============================================================

topK_scfpcde_curves_traj1 <- get_top_scfpcde_by_Dobs(
  res = traj_results$traj1,
  top_k = K_top,
  q_cut = fdr_cut
)

topK_scfpcde_curves_traj2 <- get_top_scfpcde_by_Dobs(
  res = traj_results$traj2,
  top_k = K_top,
  q_cut = fdr_cut
)

cat("\nTop scFPC-DE curve genes, trajectory 1:\n")

print(
  topK_scfpcde_curves_traj1
)

cat("\nTop scFPC-DE curve genes, trajectory 2:\n")

print(
  topK_scfpcde_curves_traj2
)

## ============================================================
## Display individual figures in RStudio Plots pane
## ============================================================

plot_tdeg_category_fpc_score(
  plot_object = tdeg_score_traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  q_cut = fdr_cut,
  point_cex = score_point_cex,
  show_null_fill = TRUE,
  standardize_plot = standardize_score_plot,
  xlim = score_xlim,
  ylim = score_ylim,
  legend_pos = "topright"
)

plot_tdeg_category_fpc_score(
  plot_object = tdeg_score_traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  q_cut = fdr_cut,
  point_cex = score_point_cex,
  show_null_fill = TRUE,
  standardize_plot = standardize_score_plot,
  xlim = score_xlim,
  ylim = score_ylim,
  legend_pos = "topleft"
)

plot_fdr_histogram_single(
  qvalues = scfpcde_table_traj1$q_value,
  method_label = "scFPC-DE",
  trajectory_label = "Trajectory 1",
  q_cut = fdr_cut,
  ylim = hist_ylim_traj1
)

plot_fdr_histogram_single(
  qvalues = pseudotimeDE_table_traj1$q_value,
  method_label = "PseudotimeDE",
  trajectory_label = "Trajectory 1",
  q_cut = fdr_cut,
  ylim = hist_ylim_traj1
)

plot_fdr_histogram_single(
  qvalues = scfpcde_table_traj2$q_value,
  method_label = "scFPC-DE",
  trajectory_label = "Trajectory 2",
  q_cut = fdr_cut,
  ylim = hist_ylim_traj2
)

plot_fdr_histogram_single(
  qvalues = pseudotimeDE_table_traj2$q_value,
  method_label = "PseudotimeDE",
  trajectory_label = "Trajectory 2",
  q_cut = fdr_cut,
  ylim = hist_ylim_traj2
)

plot_gene_panels_base_absolute(
  fpca.obj = traj_results$traj1$fpca.res,
  label_vec = topK_scfpcde_curves_traj1,
  cds_obj = traj_results$traj1$cds_sub2,
  nrow_panel = curve_nrow_panel,
  ncol_panel = curve_ncol_panel,
  point_cex = curve_point_cex_console,
  point_alpha = curve_point_alpha,
  curve_lwd = curve_lwd,
  zero_tick_lwd = zero_tick_lwd,
  cluster_cols = cluster_cols_traj1,
  cluster_order = cluster_order_traj1,
  main_title = "Trajectory 1: scFPC-DE top 10 curves",
  gene_title_cex = gene_title_cex,
  axis_cex = curve_axis_cex,
  axis_label_cex = curve_axis_label_cex,
  legend_cex = curve_legend_cex,
  main_title_cex = curve_main_title_cex
)

plot_gene_panels_base_absolute(
  fpca.obj = traj_results$traj2$fpca.res,
  label_vec = topK_scfpcde_curves_traj2,
  cds_obj = traj_results$traj2$cds_sub2,
  nrow_panel = curve_nrow_panel,
  ncol_panel = curve_ncol_panel,
  point_cex = curve_point_cex_console,
  point_alpha = curve_point_alpha,
  curve_lwd = curve_lwd,
  zero_tick_lwd = zero_tick_lwd,
  cluster_cols = cluster_cols_traj2,
  cluster_order = cluster_order_traj2,
  main_title = "Trajectory 2: scFPC-DE top 10 curves",
  gene_title_cex = gene_title_cex,
  axis_cex = curve_axis_cex,
  axis_label_cex = curve_axis_label_cex,
  legend_cex = curve_legend_cex,
  main_title_cex = curve_main_title_cex
)

## ============================================================
## Save individual TIFF files
## ============================================================

save_score_tiff(
  filename = "traj1_TDEG_categories_FPC_gene_scores.tiff",
  plot_object = tdeg_score_traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  legend_pos = "topright"
)

save_score_tiff(
  filename = "traj2_TDEG_categories_FPC_gene_scores.tiff",
  plot_object = tdeg_score_traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  legend_pos = "topleft"
)

save_hist_tiff(
  filename = "traj1_scFPCDE_FDR_histogram.tiff",
  qvalues = scfpcde_table_traj1$q_value,
  method_label = "scFPC-DE",
  trajectory_label = "Trajectory 1",
  ylim = hist_ylim_traj1
)

save_hist_tiff(
  filename = "traj1_PseudotimeDE_FDR_histogram.tiff",
  qvalues = pseudotimeDE_table_traj1$q_value,
  method_label = "PseudotimeDE",
  trajectory_label = "Trajectory 1",
  ylim = hist_ylim_traj1
)

save_hist_tiff(
  filename = "traj2_scFPCDE_FDR_histogram.tiff",
  qvalues = scfpcde_table_traj2$q_value,
  method_label = "scFPC-DE",
  trajectory_label = "Trajectory 2",
  ylim = hist_ylim_traj2
)

save_hist_tiff(
  filename = "traj2_PseudotimeDE_FDR_histogram.tiff",
  qvalues = pseudotimeDE_table_traj2$q_value,
  method_label = "PseudotimeDE",
  trajectory_label = "Trajectory 2",
  ylim = hist_ylim_traj2
)

save_gene_curve_tiff(
  filename = "traj1_scfpcde_top10_curves.tiff",
  fpca.obj = traj_results$traj1$fpca.res,
  label_vec = topK_scfpcde_curves_traj1,
  cds_obj = traj_results$traj1$cds_sub2,
  cluster_cols = cluster_cols_traj1,
  cluster_order = cluster_order_traj1,
  main_title = "Trajectory 1: scFPC-DE top 10 curves"
)

save_gene_curve_tiff(
  filename = "traj2_scfpcde_top10_curves.tiff",
  fpca.obj = traj_results$traj2$fpca.res,
  label_vec = topK_scfpcde_curves_traj2,
  cds_obj = traj_results$traj2$cds_sub2,
  cluster_cols = cluster_cols_traj2,
  cluster_order = cluster_order_traj2,
  main_title = "Trajectory 2: scFPC-DE top 10 curves"
)

## ============================================================
## Confirm saved files
## ============================================================

saved_files <- c(
  "traj1_TDEG_categories_FPC_gene_scores.tiff",
  "traj2_TDEG_categories_FPC_gene_scores.tiff",
  "traj1_scFPCDE_FDR_histogram.tiff",
  "traj1_PseudotimeDE_FDR_histogram.tiff",
  "traj2_scFPCDE_FDR_histogram.tiff",
  "traj2_PseudotimeDE_FDR_histogram.tiff",
  "traj1_scfpcde_top10_curves.tiff",
  "traj2_scfpcde_top10_curves.tiff"
)

cat("\nCurrent working directory:\n")

print(
  getwd()
)

cat("\nSlide size:\n")

print(
  c(
    width_in = slide_width,
    height_in = slide_height
  )
)

cat("\nWide panel size for score plots:\n")

print(
  c(
    width_in = wide_panel_width,
    height_in = wide_panel_height,
    dpi = tiff_dpi
  )
)

cat("\nCurve panel size:\n")

print(
  c(
    width_in = curve_panel_width,
    height_in = curve_panel_height,
    dpi = tiff_dpi
  )
)

cat("\nIndividual histogram size:\n")

print(
  c(
    width_in = hist_panel_width,
    height_in = hist_panel_height,
    dpi = tiff_dpi
  )
)

cat("\nHistogram y-axis limits:\n")

print(
  list(
    trajectory_1 = hist_ylim_traj1,
    trajectory_2 = hist_ylim_traj2
  )
)

cat("\nFunctional curve panel style:\n")

print(
  list(
    pty = "m",
    gene_title_cex = gene_title_cex,
    curve_legend_cex = curve_legend_cex,
    curve_main_title_cex = curve_main_title_cex,
    curve_lwd = curve_lwd,
    curve_axis_cex = curve_axis_cex,
    curve_axis_label_cex = curve_axis_label_cex,
    curve_panel_width = curve_panel_width,
    curve_panel_height = curve_panel_height
  )
)

cat("\nBH-adjusted TDEG counts:\n")

print(
  bh_count_table
)

cat("\nBH-adjusted TDEG overlap counts:\n")

print(
  bh_overlap_table
)

cat("\nSaved TIFF files:\n")

print(
  normalizePath(
    saved_files,
    mustWork = FALSE
  )
)

cat("\nFile exists check:\n")

print(
  file.exists(saved_files)
)
