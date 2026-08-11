## ============================================================
## Independent base R FPC gene-score plots and
## FDR-corrected p-value histograms
##
## For each trajectory:
##
##   1) scFPC-DE FPC1-FPC2 gene-score plot
##   2) PseudotimeDE FPC1-FPC2 gene-score plot
##   3) scFPC-DE FDR-corrected p-value histogram
##   4) PseudotimeDE FDR-corrected p-value histogram
##
## All TIFF files are saved directly in the current working directory.
##
## Required objects:
##   traj_results$traj1
##   traj_results$traj2
##
## Saves:
##   traj1_scFPCDE_FPC_gene_scores.tiff
##   traj1_PseudotimeDE_FPC_gene_scores.tiff
##   traj1_scFPCDE_FDR_histogram.tiff
##   traj1_PseudotimeDE_FDR_histogram.tiff
##
##   traj2_scFPCDE_FPC_gene_scores.tiff
##   traj2_PseudotimeDE_FPC_gene_scores.tiff
##   traj2_scFPCDE_FDR_histogram.tiff
##   traj2_PseudotimeDE_FDR_histogram.tiff
## ============================================================

## ============================================================
## User settings
## ============================================================

fdr_cut <- 0.05

ellipse_level <- 0.99

hist_breaks <- seq(
  0,
  1,
  length.out = 41
)

score_point_cex <- 0.55

deg_point_col <- "black"
null_point_col <- "gray80"

null_region_fill <- grDevices::adjustcolor(
  "mistyrose",
  alpha.f = 0.35
)

null_region_border <- "red3"

## TIFF settings: save directly in current working directory
tiff_width <- 8
tiff_height <- 6
tiff_dpi <- 300

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
## Extract FPCA gene-score matrix
## ============================================================

get_fpca_scores <- function(fpca.res) {
  
  if (is.null(fpca.res$fpca_result)) {
    stop("fpca.res$fpca_result was not found")
  }
  
  scores <- NULL
  
  if (!is.null(fpca.res$fpca_result$scores)) {
    
    scores <- fpca.res$fpca_result$scores
    
  } else if (!is.null(fpca.res$fpca_result$fda_fpca$scores)) {
    
    scores <- fpca.res$fpca_result$fda_fpca$scores
  }
  
  if (is.null(scores)) {
    
    stop(
      paste0(
        "No FPCA score matrix was found in:\n",
        "fpca_result$scores or ",
        "fpca_result$fda_fpca$scores"
      )
    )
  }
  
  scores <- as.matrix(scores)
  storage.mode(scores) <- "numeric"
  
  if (ncol(scores) < 2) {
    stop("The FPCA score matrix must contain at least two columns")
  }
  
  scores
}

## ============================================================
## Extract gene IDs corresponding to FPCA score rows
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
    
  } else if (
    !is.null(res$fpca.res$D_test_result$ID) &&
    length(res$fpca.res$D_test_result$ID) == nrow(scores)
  ) {
    
    gene_ids <- res$fpca.res$D_test_result$ID
  }
  
  if (is.null(gene_ids)) {
    
    stop(
      paste0(
        "Could not identify the genes corresponding to the ",
        "rows of the FPCA score matrix"
      )
    )
  }
  
  toupper(
    trimws(
      as.character(gene_ids)
    )
  )
}

## ============================================================
## Extract scFPC-DE gene-level p-values and FDR values
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
  
  if (!"ID" %in% names(sc_tbl)) {
    stop("D_test_result does not contain ID")
  }
  
  if ("p_value" %in% names(sc_tbl)) {
    
    pvalues <- suppressWarnings(
      as.numeric(sc_tbl$p_value)
    )
    
  } else {
    
    pvalues <- rep(
      NA_real_,
      nrow(sc_tbl)
    )
  }
  
  if ("q_value" %in% names(sc_tbl)) {
    
    qvalues <- suppressWarnings(
      as.numeric(sc_tbl$q_value)
    )
    
  } else if ("p_value" %in% names(sc_tbl)) {
    
    qvalues <- stats::p.adjust(
      pvalues,
      method = "BH"
    )
    
  } else {
    
    stop("D_test_result contains neither q_value nor p_value")
  }
  
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
  
  ord <- order(
    out$q_value,
    out$p_value,
    -out$D_obs,
    na.last = TRUE
  )
  
  out <- out[
    ord,
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
## Extract PseudotimeDE gene-level p-values and FDR values
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
  
  if (
    !is.null(res$gauss.qvals) &&
    length(res$gauss.qvals) == nrow(pt_tbl)
  ) {
    
    qvalues <- suppressWarnings(
      as.numeric(res$gauss.qvals)
    )
    
  } else {
    
    qvalues <- stats::p.adjust(
      pvalues,
      method = "BH"
    )
  }
  
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
  
  ord <- order(
    out$q_value,
    out$p_value,
    na.last = TRUE
  )
  
  out <- out[
    ord,
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
## Prepare score and FDR data for one trajectory
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
  
  keep <- is.finite(score_tbl$FPC1) &
    is.finite(score_tbl$FPC2)
  
  score_tbl <- score_tbl[
    keep,
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
## Estimate the scFPC-DE null-region ellipse
##
## The ellipse is estimated from genes with:
##   scFPC-DE q-value >= fdr_cut
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
      "Too few scFPC-DE nonsignificant genes to estimate the null ellipse"
    )
    
    return(NULL)
  }
  
  mu <- colMeans(X0)
  S <- stats::cov(X0)
  
  if (
    any(!is.finite(S)) ||
    nrow(S) != 2 ||
    ncol(S) != 2
  ) {
    
    warning("The scFPC-DE null covariance matrix is invalid")
    
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
      nrow = 2,
      ncol = 2
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
      nrow = 2,
      ncol = 2
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
## Prepare trajectory-specific score plot object
## ============================================================

prepare_independent_fpc_fdr_plots <- function(res,
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
  
  x_values <- score_tbl$FPC1
  y_values <- score_tbl$FPC2
  
  if (!is.null(null_ellipse)) {
    
    x_values <- c(
      x_values,
      null_ellipse[, 1]
    )
    
    y_values <- c(
      y_values,
      null_ellipse[, 2]
    )
  }
  
  xlim <- range(
    x_values,
    finite = TRUE
  )
  
  ylim <- range(
    y_values,
    finite = TRUE
  )
  
  x_pad <- 0.06 * diff(xlim)
  y_pad <- 0.06 * diff(ylim)
  
  if (!is.finite(x_pad) || x_pad <= 0) {
    x_pad <- 0.5
  }
  
  if (!is.finite(y_pad) || y_pad <= 0) {
    y_pad <- 0.5
  }
  
  xlim <- xlim + c(
    -x_pad,
    x_pad
  )
  
  ylim <- ylim + c(
    -y_pad,
    y_pad
  )
  
  list(
    plot_data = plot_data,
    score_tbl = score_tbl,
    null_ellipse = null_ellipse,
    xlim = xlim,
    ylim = ylim
  )
}

## ============================================================
## Independent FPC gene-score plot
## ============================================================

plot_independent_fpc_score <- function(plot_object,
                                       method = c(
                                         "scFPC-DE",
                                         "PseudotimeDE"
                                       ),
                                       trajectory_title,
                                       q_cut = 0.05,
                                       point_cex = 0.55,
                                       show_null_fill = TRUE) {
  
  method <- match.arg(method)
  
  score_tbl <- plot_object$score_tbl
  null_ellipse <- plot_object$null_ellipse
  
  qvalues <- if (method == "scFPC-DE") {
    score_tbl$scFPCDE_qvalue
  } else {
    score_tbl$PseudotimeDE_qvalue
  }
  
  is_deg <- is.finite(qvalues) &
    qvalues < q_cut
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  par(
    mar = c(5.2, 5.2, 4.6, 2.0),
    mgp = c(2.5, 0.75, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    xaxs = "r",
    yaxs = "r"
  )
  
  plot(
    score_tbl$FPC1,
    score_tbl$FPC2,
    type = "n",
    xlim = plot_object$xlim,
    ylim = plot_object$ylim,
    xlab = "FPC1",
    ylab = "FPC2",
    main = paste0(
      trajectory_title,
      "\n",
      method,
      " gene scores colored by test decision"
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
  
  points(
    score_tbl$FPC1[!is_deg],
    score_tbl$FPC2[!is_deg],
    pch = 16,
    cex = point_cex,
    col = null_point_col
  )
  
  points(
    score_tbl$FPC1[is_deg],
    score_tbl$FPC2[is_deg],
    pch = 16,
    cex = point_cex,
    col = deg_point_col
  )
  
  abline(
    h = 0,
    v = 0,
    lty = 2,
    lwd = 1,
    col = "gray55"
  )
  
  if (!is.null(null_ellipse)) {
    
    lines(
      null_ellipse[, 1],
      null_ellipse[, 2],
      col = null_region_border,
      lty = 2,
      lwd = 2
    )
  }
  
  legend(
    "topright",
    legend = c(
      paste0(
        method,
        " significant gene"
      ),
      "Nonsignificant gene",
      "scFPC-DE null region"
    ),
    col = c(
      deg_point_col,
      null_point_col,
      null_region_border
    ),
    pch = c(
      16,
      16,
      NA
    ),
    lty = c(
      NA,
      NA,
      2
    ),
    lwd = c(
      NA,
      NA,
      2
    ),
    pt.cex = 1.1,
    bty = "n",
    cex = 0.85
  )
  
  invisible(
    list(
      method = method,
      qvalues = qvalues,
      is_deg = is_deg,
      significant_genes = score_tbl$Gene_ID[is_deg]
    )
  )
}

## ============================================================
## Independent FDR-corrected p-value histogram
## ============================================================

plot_independent_fdr_histogram <- function(qvalues,
                                           method = c(
                                             "scFPC-DE",
                                             "PseudotimeDE"
                                           ),
                                           trajectory_title,
                                           q_cut = 0.05,
                                           breaks = hist_breaks,
                                           hist_col = "white",
                                           border_col = "black") {
  
  method <- match.arg(method)
  
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
    mar = c(5.2, 5.2, 4.6, 2.0),
    mgp = c(2.5, 0.75, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    xaxs = "i"
  )
  
  if (length(qvalues) == 0) {
    
    plot.new()
    
    title(
      main = paste0(
        trajectory_title,
        "\n",
        method,
        " FDR-corrected p-values"
      )
    )
    
    text(
      x = 0.5,
      y = 0.5,
      labels = "No valid FDR values",
      cex = 1.1
    )
    
    return(
      invisible(NULL)
    )
  }
  
  hist_obj <- hist(
    qvalues,
    breaks = breaks,
    col = hist_col,
    border = border_col,
    main = paste0(
      trajectory_title,
      "\n",
      method,
      " FDR-corrected p-values"
    ),
    xlab = "FDR-corrected p-value",
    ylab = "Count",
    xlim = c(0, 1),
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
    )
  )
  
  abline(
    v = q_cut,
    col = "red3",
    lty = 2,
    lwd = 2
  )
  
  n_significant <- sum(
    qvalues < q_cut
  )
  
  legend(
    "topright",
    legend = c(
      paste0(
        "FDR threshold = ",
        q_cut
      ),
      paste0(
        "Significant genes = ",
        n_significant
      ),
      paste0(
        "Total genes = ",
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
      2,
      NA,
      NA
    ),
    bty = "n",
    cex = 0.90
  )
  
  invisible(
    list(
      histogram = hist_obj,
      qvalues = qvalues,
      significant_count = n_significant
    )
  )
}

## ============================================================
## Prepare both trajectories
## ============================================================

independent_fpc_fdr_traj1 <- prepare_independent_fpc_fdr_plots(
  res = traj_results$traj1,
  q_cut = fdr_cut,
  ellipse_level = ellipse_level
)

independent_fpc_fdr_traj2 <- prepare_independent_fpc_fdr_plots(
  res = traj_results$traj2,
  q_cut = fdr_cut,
  ellipse_level = ellipse_level
)

## ============================================================
## Display all eight figures independently in RStudio
## ============================================================

traj1_scfpcde_score_plot <- plot_independent_fpc_score(
  plot_object = independent_fpc_fdr_traj1,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  q_cut = fdr_cut,
  point_cex = score_point_cex
)

traj1_ptime_score_plot <- plot_independent_fpc_score(
  plot_object = independent_fpc_fdr_traj1,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  q_cut = fdr_cut,
  point_cex = score_point_cex
)

traj1_scfpcde_fdr_plot <- plot_independent_fdr_histogram(
  qvalues = independent_fpc_fdr_traj1$plot_data$scFPCDE_table$q_value,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  q_cut = fdr_cut,
  breaks = hist_breaks
)

traj1_ptime_fdr_plot <- plot_independent_fdr_histogram(
  qvalues = independent_fpc_fdr_traj1$plot_data$PseudotimeDE_table$q_value,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  q_cut = fdr_cut,
  breaks = hist_breaks
)

traj2_scfpcde_score_plot <- plot_independent_fpc_score(
  plot_object = independent_fpc_fdr_traj2,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  q_cut = fdr_cut,
  point_cex = score_point_cex
)

traj2_ptime_score_plot <- plot_independent_fpc_score(
  plot_object = independent_fpc_fdr_traj2,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  q_cut = fdr_cut,
  point_cex = score_point_cex
)

traj2_scfpcde_fdr_plot <- plot_independent_fdr_histogram(
  qvalues = independent_fpc_fdr_traj2$plot_data$scFPCDE_table$q_value,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  q_cut = fdr_cut,
  breaks = hist_breaks
)

traj2_ptime_fdr_plot <- plot_independent_fdr_histogram(
  qvalues = independent_fpc_fdr_traj2$plot_data$PseudotimeDE_table$q_value,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  q_cut = fdr_cut,
  breaks = hist_breaks
)

## ============================================================
## Save independent FPC score plot as TIFF
## ============================================================

save_independent_fpc_score_tiff <- function(filename,
                                            plot_object,
                                            method,
                                            trajectory_title,
                                            q_cut = fdr_cut,
                                            point_cex = score_point_cex,
                                            width = tiff_width,
                                            height = tiff_height,
                                            dpi = tiff_dpi) {
  
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
  
  plot_independent_fpc_score(
    plot_object = plot_object,
    method = method,
    trajectory_title = trajectory_title,
    q_cut = q_cut,
    point_cex = point_cex
  )
  
  invisible(filename)
}

## ============================================================
## Save independent FDR histogram as TIFF
## ============================================================

save_independent_fdr_histogram_tiff <- function(filename,
                                                qvalues,
                                                method,
                                                trajectory_title,
                                                q_cut = fdr_cut,
                                                breaks = hist_breaks,
                                                width = tiff_width,
                                                height = tiff_height,
                                                dpi = tiff_dpi) {
  
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
  
  plot_independent_fdr_histogram(
    qvalues = qvalues,
    method = method,
    trajectory_title = trajectory_title,
    q_cut = q_cut,
    breaks = breaks
  )
  
  invisible(filename)
}

## ============================================================
## Save all eight TIFF files directly in current working directory
## ============================================================

saved_tiff_files <- c(
  "traj1_scFPCDE_FPC_gene_scores.tiff",
  "traj1_PseudotimeDE_FPC_gene_scores.tiff",
  "traj1_scFPCDE_FDR_histogram.tiff",
  "traj1_PseudotimeDE_FDR_histogram.tiff",
  "traj2_scFPCDE_FPC_gene_scores.tiff",
  "traj2_PseudotimeDE_FPC_gene_scores.tiff",
  "traj2_scFPCDE_FDR_histogram.tiff",
  "traj2_PseudotimeDE_FDR_histogram.tiff"
)

save_independent_fpc_score_tiff(
  filename = saved_tiff_files[1],
  plot_object = independent_fpc_fdr_traj1,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  q_cut = fdr_cut,
  point_cex = score_point_cex
)

save_independent_fpc_score_tiff(
  filename = saved_tiff_files[2],
  plot_object = independent_fpc_fdr_traj1,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  q_cut = fdr_cut,
  point_cex = score_point_cex
)

save_independent_fdr_histogram_tiff(
  filename = saved_tiff_files[3],
  qvalues = independent_fpc_fdr_traj1$plot_data$scFPCDE_table$q_value,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  q_cut = fdr_cut,
  breaks = hist_breaks
)

save_independent_fdr_histogram_tiff(
  filename = saved_tiff_files[4],
  qvalues = independent_fpc_fdr_traj1$plot_data$PseudotimeDE_table$q_value,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  q_cut = fdr_cut,
  breaks = hist_breaks
)

save_independent_fpc_score_tiff(
  filename = saved_tiff_files[5],
  plot_object = independent_fpc_fdr_traj2,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  q_cut = fdr_cut,
  point_cex = score_point_cex
)

save_independent_fpc_score_tiff(
  filename = saved_tiff_files[6],
  plot_object = independent_fpc_fdr_traj2,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  q_cut = fdr_cut,
  point_cex = score_point_cex
)

save_independent_fdr_histogram_tiff(
  filename = saved_tiff_files[7],
  qvalues = independent_fpc_fdr_traj2$plot_data$scFPCDE_table$q_value,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  q_cut = fdr_cut,
  breaks = hist_breaks
)

save_independent_fdr_histogram_tiff(
  filename = saved_tiff_files[8],
  qvalues = independent_fpc_fdr_traj2$plot_data$PseudotimeDE_table$q_value,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  q_cut = fdr_cut,
  breaks = hist_breaks
)

## ============================================================
## Confirm saved files
## ============================================================

cat("\nCurrent working directory:\n")
print(getwd())

cat("\nSaved TIFF files:\n")
print(
  normalizePath(
    saved_tiff_files,
    mustWork = FALSE
  )
)

cat("\nFile exists check:\n")
print(
  file.exists(saved_tiff_files)
)
