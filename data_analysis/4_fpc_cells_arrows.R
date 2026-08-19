## ============================================================
## Observed-cell FPC projection with fitted trajectory curve
##   - Cells:        Y %*% gene FPC scores
##   - Smooth curve: xt_hat %*% gene FPC scores
##   - Arrows:       finite differences along fitted curve
##
## FPC dimensions:
##   - FPC1 vs FPC2
##   - FPC1 vs FPC3
##   - FPC2 vs FPC3
##
## Saves:
##   traj1_observed_cell_fpc1_fpc2_fitted_path.tiff
##   traj1_observed_cell_fpc1_fpc3_fitted_path.tiff
##   traj1_observed_cell_fpc2_fpc3_fitted_path.tiff
##   traj2_observed_cell_fpc1_fpc2_fitted_path.tiff
##   traj2_observed_cell_fpc1_fpc3_fitted_path.tiff
##   traj2_observed_cell_fpc2_fpc3_fitted_path.tiff
## ============================================================

get_analysis_clusters <- function(res) {
  if (!is.null(res$tt.cluster)) {
    return(res$tt.cluster)
  }

  cds_obj <- res$cds_sub2
  if (is.null(cds_obj)) {
    stop("No cluster vector or cds_sub2 object was found")
  }

  cell_data <- SummarizedExperiment::colData(cds_obj)
  for (field in c("cluster", "cell_type")) {
    if (field %in% colnames(cell_data)) {
      return(cell_data[[field]])
    }
  }

  if (inherits(cds_obj, "cell_data_set") &&
      requireNamespace("monocle3", quietly = TRUE)) {
    return(monocle3::clusters(cds_obj))
  }

  stop("Could not identify cell clusters in the analysis result")
}
## ------------------------------------------------------------
## Cluster settings
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## FPC pairs
## ------------------------------------------------------------

fpc_pairs <- list(
  fpc1_fpc2 = c(1, 2),
  fpc1_fpc3 = c(1, 3),
  fpc2_fpc3 = c(2, 3)
)

## ------------------------------------------------------------
## Panel size for a 29 x 23 inch 3 x 3 multipanel figure
## ------------------------------------------------------------

panel_width <- 29 / 3
panel_height <- 23 / 3

## ============================================================
## Finite-difference derivative
## ============================================================

gradient_fd <- function(tt, y) {
  
  n <- length(tt)
  
  if (n != length(y)) {
    stop("tt and y must have the same length")
  }
  
  if (n < 3) {
    stop("At least three time points are required")
  }
  
  dy <- numeric(n)
  
  dy[1] <- (y[2] - y[1]) /
    (tt[2] - tt[1])
  
  dy[n] <- (y[n] - y[n - 1]) /
    (tt[n] - tt[n - 1])
  
  for (i in 2:(n - 1)) {
    
    dy[i] <- (y[i + 1] - y[i - 1]) /
      (tt[i + 1] - tt[i - 1])
  }
  
  dy[!is.finite(dy)] <- 0
  
  dy
}

## ============================================================
## Extract observed and fitted projections for selected FPC pair
## ============================================================

get_cell_and_fitted_fpc_projection <- function(res,
                                               pc_pair = c(1, 2)) {
  
  if (length(pc_pair) != 2) {
    stop("pc_pair must contain exactly two FPC indices")
  }
  
  if (any(pc_pair < 1)) {
    stop("pc_pair must contain positive FPC indices")
  }
  
  ## observed centered data: cells x genes
  Y <- as.matrix(
    res$y
  )
  
  ## fitted centered/smoothed expression: cells x genes
  Xhat <- as.matrix(
    res$fpca.res$fpca_result$xt_hat
  )
  
  ## gene-level FPC scores: genes x PCs
  scores <- res$fpca.res$fpca_result$scores
  
  if (is.null(scores)) {
    scores <- res$fpca.res$fpca_result$fda_fpca$scores
  }
  
  if (is.null(scores)) {
    stop("Could not find FPC scores in res$fpca.res")
  }
  
  scores <- as.matrix(
    scores
  )
  
  if (ncol(scores) < max(pc_pair)) {
    stop(
      "FPC score matrix has fewer columns than requested pc_pair. Available PCs = ",
      ncol(scores)
    )
  }
  
  genes_y <- colnames(Y)
  genes_x <- colnames(Xhat)
  genes_s <- rownames(scores)
  
  if (is.null(genes_y)) {
    stop("res$y must have gene names as column names")
  }
  
  if (is.null(genes_x)) {
    stop("xt_hat must have gene names as column names")
  }
  
  if (is.null(genes_s)) {
    stop("FPC score matrix must have gene names as row names")
  }
  
  common_genes <- Reduce(
    intersect,
    list(
      genes_y,
      genes_x,
      genes_s
    )
  )
  
  if (length(common_genes) < 5) {
    stop("Too few genes overlap between res$y, xt_hat, and FPC scores")
  }
  
  Y_use <- Y[
    ,
    common_genes,
    drop = FALSE
  ]
  
  Xhat_use <- Xhat[
    ,
    common_genes,
    drop = FALSE
  ]
  
  S_use <- scores[
    common_genes,
    pc_pair,
    drop = FALSE
  ]
  
  ## normalize selected FPC score directions
  score_norm <- sqrt(
    colSums(S_use^2)
  )
  
  score_norm[
    score_norm == 0 |
      !is.finite(score_norm)
  ] <- 1
  
  S_use <- sweep(
    S_use,
    MARGIN = 2,
    STATS = score_norm,
    FUN = "/"
  )
  
  S_use[!is.finite(S_use)] <- 0
  
  ## observed cell projection
  Z_obs <- Y_use %*% S_use
  
  ## fitted trajectory projection
  Z_fit <- Xhat_use %*% S_use
  
  colnames(Z_obs) <- paste0(
    "FPC",
    pc_pair
  )
  
  colnames(Z_fit) <- paste0(
    "FPC",
    pc_pair
  )
  
  tt <- as.numeric(
    res$tt
  )
  
  if (length(tt) != nrow(Y)) {
    stop("Length of res$tt does not match number of rows in res$y")
  }
  
  clu <- get_analysis_clusters(res)
  
  clu <- as.character(
    clu
  )
  
  if (length(clu) != nrow(Y)) {
    stop("Cluster vector length does not match number of rows in res$y")
  }
  
  ord <- order(
    tt
  )
  
  list(
    Z_obs = Z_obs[
      ord,
      ,
      drop = FALSE
    ],
    Z_fit = Z_fit[
      ord,
      ,
      drop = FALSE
    ],
    tt = tt[ord],
    cluster = clu[ord],
    genes_used = common_genes,
    pc_pair = pc_pair
  )
}

## ============================================================
## Cluster legend outside plot, below x-axis title
## ============================================================

draw_bottom_cluster_legend <- function(clu_levels,
                                       cluster_cols_use,
                                       cex = 1.20,
                                       point_cex = 1.70) {
  
  usr <- par("usr")
  
  x_span <- usr[2] - usr[1]
  y_span <- usr[4] - usr[3]
  
  legend_y <- usr[3] - 0.260 * y_span
  
  legend_x <- seq(
    usr[1] + 0.25 * x_span,
    usr[2] - 0.25 * x_span,
    length.out = length(clu_levels)
  )
  
  par(xpd = NA)
  
  for (j in seq_along(clu_levels)) {
    
    points(
      x = legend_x[j],
      y = legend_y,
      pch = 16,
      col = unname(
        cluster_cols_use[clu_levels[j]]
      ),
      cex = point_cex
    )
    
    text(
      x = legend_x[j] + 0.022 * x_span,
      y = legend_y,
      labels = clu_levels[j],
      adj = c(0, 0.5),
      cex = cex,
      col = "black"
    )
  }
  
  par(xpd = FALSE)
}

## ============================================================
## Model legend inside plot box, bottom-left
## ============================================================

draw_inside_model_legend <- function(path_lwd = 4.0,
                                     cex = 1.15) {
  
  usr <- par("usr")
  
  x_span <- usr[2] - usr[1]
  y_span <- usr[4] - usr[3]
  
  line_y <- usr[3] + 0.115 * y_span
  arrow_y <- usr[3] + 0.060 * y_span
  
  x0 <- usr[1] + 0.055 * x_span
  x1 <- usr[1] + 0.135 * x_span
  text_x <- usr[1] + 0.155 * x_span
  
  par(xpd = FALSE)
  
  ## Fitted trajectory line
  segments(
    x0 = x0,
    y0 = line_y,
    x1 = x1,
    y1 = line_y,
    col = "black",
    lwd = path_lwd
  )
  
  text(
    x = text_x,
    y = line_y,
    labels = "Fitted trajectory",
    adj = c(0, 0.5),
    cex = cex,
    col = "black"
  )
  
  ## Gene expression direction arrow
  arrows(
    x0 = x0,
    y0 = arrow_y,
    x1 = x1,
    y1 = arrow_y,
    length = 0.14,
    angle = 30,
    col = "black",
    lwd = 3.5
  )
  
  text(
    x = text_x,
    y = arrow_y,
    labels = "Gene expression direction",
    adj = c(0, 0.5),
    cex = cex,
    col = "black"
  )
}

## ============================================================
## Plot observed cells + fitted FPC trajectory + arrows
## ============================================================

plot_cell_fpc_fitted_path <- function(res,
                                      pc_pair = c(1, 2),
                                      traj_title = "Trajectory",
                                      cluster_order = NULL,
                                      cluster_cols = NULL,
                                      arrow_n = 18,
                                      arrow_scale = 0.45,
                                      point_cex = 0.95,
                                      point_alpha = 0.75,
                                      path_lwd = 4.0,
                                      arrow_lwd = 4.0,
                                      legend_cex = 1.20,
                                      axis_cex = 1.25,
                                      label_cex = 1.40,
                                      title_cex = 1.45) {
  
  proj <- get_cell_and_fitted_fpc_projection(
    res = res,
    pc_pair = pc_pair
  )
  
  Z_obs <- proj$Z_obs
  Z_fit <- proj$Z_fit
  tt <- proj$tt
  clu_raw <- proj$cluster
  
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
  
  if (anyNA(clu)) {
    
    missing_clu <- sort(
      unique(
        clu_raw[is.na(clu)]
      )
    )
    
    stop(
      "Some clusters were not included in cluster_order: ",
      paste(missing_clu, collapse = ", ")
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
    
    cluster_cols <- setNames(
      cluster_cols[seq_along(clu_levels)],
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
      paste(missing_cols, collapse = ", ")
    )
  }
  
  cluster_cols_use <- cluster_cols[
    clu_levels
  ]
  
  point_cols <- grDevices::adjustcolor(
    unname(cluster_cols_use),
    alpha.f = point_alpha
  )
  
  names(point_cols) <- clu_levels
  
  cell_cols <- point_cols[
    as.character(clu)
  ]
  
  if (anyNA(cell_cols)) {
    stop("Cell colors contain NA values")
  }
  
  ## ----------------------------------------------------------
  ## Arrow directions from fitted trajectory
  ## ----------------------------------------------------------
  
  dZ1 <- gradient_fd(
    tt,
    Z_fit[, 1]
  )
  
  dZ2 <- gradient_fd(
    tt,
    Z_fit[, 2]
  )
  
  arrow_idx <- unique(
    round(
      seq(
        1,
        length(tt),
        length.out = arrow_n
      )
    )
  )
  
  x0 <- Z_fit[
    arrow_idx,
    1
  ]
  
  y0 <- Z_fit[
    arrow_idx,
    2
  ]
  
  dx <- dZ1[
    arrow_idx
  ]
  
  dy <- dZ2[
    arrow_idx
  ]
  
  norm_arrow <- sqrt(
    dx^2 + dy^2
  )
  
  norm_arrow[
    norm_arrow == 0 |
      !is.finite(norm_arrow)
  ] <- 1
  
  dx <- arrow_scale * dx /
    norm_arrow
  
  dy <- arrow_scale * dy /
    norm_arrow
  
  ## ----------------------------------------------------------
  ## Plot limits
  ## ----------------------------------------------------------
  
  x_rng <- range(
    c(
      Z_obs[, 1],
      Z_fit[, 1]
    ),
    finite = TRUE
  )
  
  y_rng <- range(
    c(
      Z_obs[, 2],
      Z_fit[, 2]
    ),
    finite = TRUE
  )
  
  x_pad <- 0.12 * diff(x_rng)
  y_pad <- 0.12 * diff(y_rng)
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  par(
    mar = c(8.2, 6.0, 4.4, 1.6),
    mgp = c(3.2, 1.0, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    col.sub = "black",
    cex.axis = axis_cex,
    cex.lab = label_cex,
    cex.main = title_cex,
    font.main = 2,
    lwd = 1.5,
    xpd = FALSE
  )
  
  plot(
    Z_obs[, 1],
    Z_obs[, 2],
    type = "n",
    xlab = "",
    ylab = "",
    main = "",
    xlim = c(
      x_rng[1] - x_pad,
      x_rng[2] + x_pad
    ),
    ylim = c(
      y_rng[1] - y_pad,
      y_rng[2] + y_pad
    ),
    axes = FALSE,
    ann = FALSE
  )
  
  usr <- par("usr")
  
  rect(
    usr[1],
    usr[3],
    usr[2],
    usr[4],
    col = "white",
    border = NA
  )
  
  points(
    Z_obs[, 1],
    Z_obs[, 2],
    col = cell_cols,
    pch = 16,
    cex = point_cex
  )
  
  lines(
    Z_fit[, 1],
    Z_fit[, 2],
    col = "black",
    lwd = path_lwd
  )
  
  arrows(
    x0 = x0,
    y0 = y0,
    x1 = x0 + dx,
    y1 = y0 + dy,
    length = 0.18,
    angle = 30,
    col = "black",
    lwd = arrow_lwd
  )
  
  axis(
    side = 1,
    col = "black",
    col.axis = "black",
    cex.axis = axis_cex,
    lwd = 1.5
  )
  
  axis(
    side = 2,
    col = "black",
    col.axis = "black",
    las = 1,
    cex.axis = axis_cex,
    lwd = 1.5
  )
  
  box(
    col = "black",
    lwd = 1.5
  )
  
  title(
    main = paste0(
      traj_title,
      ": FPC",
      pc_pair[1],
      "-FPC",
      pc_pair[2],
      " Cell Trajectory"
    ),
    ylab = paste0(
      "Cell projection on FPC",
      pc_pair[2]
    ),
    col.main = "black",
    col.lab = "black",
    cex.main = title_cex,
    cex.lab = label_cex,
    line = 2.0
  )
  
  mtext(
    text = paste0(
      "Cell projection on FPC",
      pc_pair[1]
    ),
    side = 1,
    line = 2.85,
    cex = label_cex
  )
  
  draw_inside_model_legend(
    path_lwd = path_lwd,
    cex = legend_cex
  )
  
  draw_bottom_cluster_legend(
    clu_levels = clu_levels,
    cluster_cols_use = cluster_cols_use,
    cex = legend_cex,
    point_cex = 1.70
  )
  
  invisible(
    proj
  )
}

## ============================================================
## Plot all combinations in RStudio Plots pane
## ============================================================

cell_path_traj1_fpc1_fpc2 <- plot_cell_fpc_fitted_path(
  res = traj_results$traj1,
  pc_pair = c(1, 2),
  traj_title = "Trajectory 1: Trans / Naive / C-mem1",
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1
)

cell_path_traj1_fpc1_fpc3 <- plot_cell_fpc_fitted_path(
  res = traj_results$traj1,
  pc_pair = c(1, 3),
  traj_title = "Trajectory 1: Trans / Naive / C-mem1",
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1
)

cell_path_traj1_fpc2_fpc3 <- plot_cell_fpc_fitted_path(
  res = traj_results$traj1,
  pc_pair = c(2, 3),
  traj_title = "Trajectory 1: Trans / Naive / C-mem1",
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1
)

cell_path_traj2_fpc1_fpc2 <- plot_cell_fpc_fitted_path(
  res = traj_results$traj2,
  pc_pair = c(1, 2),
  traj_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2
)

cell_path_traj2_fpc1_fpc3 <- plot_cell_fpc_fitted_path(
  res = traj_results$traj2,
  pc_pair = c(1, 3),
  traj_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2
)

cell_path_traj2_fpc2_fpc3 <- plot_cell_fpc_fitted_path(
  res = traj_results$traj2,
  pc_pair = c(2, 3),
  traj_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2
)

## ============================================================
## Save TIFF helper
## ============================================================

save_tiff_plot <- function(filename,
                           plot_fun,
                           width = panel_width,
                           height = panel_height,
                           dpi = 300,
                           ...) {
  
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
  
  plot_fun(...)
}

## ============================================================
## Save TIFF files
## ============================================================

saved_files <- character()

for (traj_name in c("traj1", "traj2")) {
  
  if (traj_name == "traj1") {
    
    res_use <- traj_results$traj1
    traj_title_use <- "Trajectory 1: Trans / Naive / C-mem1"
    cluster_order_use <- cluster_order_traj1
    cluster_cols_use <- cluster_cols_traj1
    
  } else {
    
    res_use <- traj_results$traj2
    traj_title_use <- "Trajectory 2: M-mem1 / DN3 / DN2"
    cluster_order_use <- cluster_order_traj2
    cluster_cols_use <- cluster_cols_traj2
  }
  
  for (pair_name in names(fpc_pairs)) {
    
    pc_pair_use <- fpc_pairs[[pair_name]]
    
    filename_use <- paste0(
      traj_name,
      "_observed_cell_",
      pair_name,
      "_fitted_path.tiff"
    )
    
    save_tiff_plot(
      filename = filename_use,
      plot_fun = plot_cell_fpc_fitted_path,
      width = panel_width,
      height = panel_height,
      dpi = 300,
      res = res_use,
      pc_pair = pc_pair_use,
      traj_title = traj_title_use,
      cluster_order = cluster_order_use,
      cluster_cols = cluster_cols_use,
      arrow_n = 18,
      arrow_scale = 0.45,
      point_cex = 0.95,
      point_alpha = 0.75,
      path_lwd = 4.0,
      arrow_lwd = 4.0,
      legend_cex = 1.20,
      axis_cex = 1.25,
      label_cex = 1.40,
      title_cex = 1.45
    )
    
    saved_files <- c(
      saved_files,
      filename_use
    )
  }
}

## ============================================================
## Confirm saved files
## ============================================================

print(
  normalizePath(
    saved_files,
    mustWork = FALSE
  )
)

print(
  file.exists(saved_files)
)
