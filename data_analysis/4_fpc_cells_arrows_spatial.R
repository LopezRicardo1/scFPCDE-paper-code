## ============================================================
## 3D observed-cell FPC projection with fitted trajectory curve
##   - Cells:        Y %*% gene FPC scores
##   - Smooth curve: xt_hat %*% gene FPC scores
##   - 3D space:     FPC1, FPC2, FPC3
##   - Arrows:       finite differences along fitted 3D curve
##
## Standardization:
##   - Coordinates are standardized to z-scores
##   - Axes are truncated to [-2, 2] by default
##
## PowerPoint:
##   - Final figure size: 30 x 24 inches
##   - Layout: 3 x 3
##   - Individual panel size: 10 x 8 inches
##
## Saves:
##   traj1_observed_cell_fpc123_fitted_path_3d.tiff
##   traj2_observed_cell_fpc123_fitted_path_3d.tiff
## ============================================================

suppressPackageStartupMessages({
  library(scatterplot3d)
})

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
## Global figure settings
## ------------------------------------------------------------

view_angle <- 300

view_angle_traj1 <- view_angle - 50
view_angle_traj2 <- view_angle

standard_axis_limit <- 2

ppt_width <- 30
ppt_height <- 24

n_col <- 3
n_row <- 3

panel_width <- ppt_width / n_col
panel_height <- ppt_height / n_row

## ------------------------------------------------------------
## Legend positions
##
## These are normalized plot-region coordinates:
##   x = 0 means left side of plot region
##   x = 1 means right side of plot region
##   y = 0 means bottom of plot region
##   y = 1 means top of plot region
##
## Change only these four vectors to move legends.
## ------------------------------------------------------------

cluster_legend_npc_traj1 <- c(0.73, 0.18)
model_legend_npc_traj1   <- c(0.10, 0.18)

cluster_legend_npc_traj2 <- c(0.73, 0.18)
model_legend_npc_traj2   <- c(0.10, 0.18)

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

## ============================================================
## Convert normalized plot-region coordinates to user coordinates
## ============================================================

npc_to_usr <- function(xy) {
  
  usr <- par("usr")
  
  x <- usr[1] + xy[1] * (usr[2] - usr[1])
  y <- usr[3] + xy[2] * (usr[4] - usr[3])
  
  c(
    x = x,
    y = y
  )
}

## ============================================================
## Custom model legend with arrow
## ============================================================

draw_model_legend_arrow <- function(model_legend_npc = c(0.10, 0.18),
                                    path_lwd = 5,
                                    arrow_lwd = 4,
                                    legend_cex = 1.5,
                                    line_len_npc = 0.080,
                                    row_gap_npc = 0.042,
                                    text_gap_npc = 0.030,
                                    arrow_head = 0.15) {
  
  usr <- par("usr")
  
  x_span <- usr[2] - usr[1]
  y_span <- usr[4] - usr[3]
  
  base_xy <- npc_to_usr(
    model_legend_npc
  )
  
  x0 <- base_xy["x"]
  y0 <- base_xy["y"]
  
  x1 <- x0 + line_len_npc * x_span
  text_x <- x1 + text_gap_npc * x_span
  
  y_fit <- y0
  y_arrow <- y0 - row_gap_npc * y_span
  
  par(
    xpd = NA
  )
  
  segments(
    x0 = x0,
    y0 = y_fit,
    x1 = x1,
    y1 = y_fit,
    col = "black",
    lwd = path_lwd,
    lend = 2
  )
  
  text(
    x = text_x,
    y = y_fit,
    labels = "Fitted trajectory",
    adj = c(0, 0.5),
    cex = legend_cex,
    col = "black"
  )
  
  arrows(
    x0 = x0,
    y0 = y_arrow,
    x1 = x1,
    y1 = y_arrow,
    length = arrow_head,
    angle = 30,
    col = "black",
    lwd = arrow_lwd
  )
  
  text(
    x = text_x,
    y = y_arrow,
    labels = "Expression direction",
    adj = c(0, 0.5),
    cex = legend_cex,
    col = "black"
  )
  
  par(
    xpd = FALSE
  )
}

## ============================================================
## Custom cluster legend
## ============================================================

draw_cluster_legend_npc <- function(clu_levels,
                                    cluster_cols_use,
                                    cluster_legend_npc = c(0.73, 0.18),
                                    legend_cex = 1.5,
                                    pt_cex = 1.4,
                                    row_gap_npc = 0.042,
                                    text_gap_npc = 0.030) {
  
  usr <- par("usr")
  
  x_span <- usr[2] - usr[1]
  y_span <- usr[4] - usr[3]
  
  base_xy <- npc_to_usr(
    cluster_legend_npc
  )
  
  x0 <- base_xy["x"]
  y0 <- base_xy["y"]
  
  text_x <- x0 + text_gap_npc * x_span
  
  par(
    xpd = NA
  )
  
  for (j in seq_along(clu_levels)) {
    
    yj <- y0 - (j - 1) * row_gap_npc * y_span
    
    points(
      x = x0,
      y = yj,
      pch = 16,
      col = unname(
        cluster_cols_use[clu_levels[j]]
      ),
      cex = pt_cex
    )
    
    text(
      x = text_x,
      y = yj,
      labels = clu_levels[j],
      adj = c(0, 0.5),
      cex = legend_cex,
      col = "black"
    )
  }
  
  par(
    xpd = FALSE
  )
}

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
## Standardize 3D FPC coordinates to z-scores
## ============================================================

standardize_fpc_coordinates_3d <- function(Z_obs,
                                           Z_fit) {
  
  Z_all <- rbind(
    Z_obs,
    Z_fit
  )
  
  center <- colMeans(
    Z_all,
    na.rm = TRUE
  )
  
  scale <- apply(
    Z_all,
    MARGIN = 2,
    FUN = sd,
    na.rm = TRUE
  )
  
  scale[
    scale == 0 |
      !is.finite(scale)
  ] <- 1
  
  Z_obs_std <- sweep(
    Z_obs,
    MARGIN = 2,
    STATS = center,
    FUN = "-"
  )
  
  Z_obs_std <- sweep(
    Z_obs_std,
    MARGIN = 2,
    STATS = scale,
    FUN = "/"
  )
  
  Z_fit_std <- sweep(
    Z_fit,
    MARGIN = 2,
    STATS = center,
    FUN = "-"
  )
  
  Z_fit_std <- sweep(
    Z_fit_std,
    MARGIN = 2,
    STATS = scale,
    FUN = "/"
  )
  
  list(
    Z_obs = Z_obs_std,
    Z_fit = Z_fit_std,
    center = center,
    scale = scale
  )
}

## ============================================================
## Extract observed and fitted 3D FPC projection
## ============================================================

get_cell_and_fitted_fpc_projection_3d <- function(res,
                                                  pc_use = c(1, 2, 3)) {
  
  if (length(pc_use) != 3) {
    stop("pc_use must contain exactly three FPC indices")
  }
  
  if (any(pc_use < 1)) {
    stop("pc_use must contain positive FPC indices")
  }
  
  Y <- as.matrix(
    res$y
  )
  
  Xhat <- as.matrix(
    res$fpca.res$fpca_result$xt_hat
  )
  
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
  
  if (ncol(scores) < max(pc_use)) {
    stop(
      "FPC score matrix has fewer columns than requested pc_use. Available PCs = ",
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
    pc_use,
    drop = FALSE
  ]
  
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
  
  Z_obs <- Y_use %*% S_use
  Z_fit <- Xhat_use %*% S_use
  
  colnames(Z_obs) <- paste0(
    "FPC",
    pc_use
  )
  
  colnames(Z_fit) <- paste0(
    "FPC",
    pc_use
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
    pc_use = pc_use
  )
}

## ============================================================
## Add fitted trajectory arrows in projected 3D space
## ============================================================

add_3d_arrows <- function(s3d,
                          Z_fit,
                          tt,
                          arrow_n = 30,
                          arrow_scale = 0.40,
                          arrow_lwd = 4,
                          arrow_col = "black",
                          arrow_head = 0.25) {
  
  dZ1 <- gradient_fd(
    tt,
    Z_fit[, 1]
  )
  
  dZ2 <- gradient_fd(
    tt,
    Z_fit[, 2]
  )
  
  dZ3 <- gradient_fd(
    tt,
    Z_fit[, 3]
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
  
  z0 <- Z_fit[
    arrow_idx,
    3
  ]
  
  dx <- dZ1[
    arrow_idx
  ]
  
  dy <- dZ2[
    arrow_idx
  ]
  
  dz <- dZ3[
    arrow_idx
  ]
  
  norm_arrow <- sqrt(
    dx^2 + dy^2 + dz^2
  )
  
  norm_arrow[
    norm_arrow == 0 |
      !is.finite(norm_arrow)
  ] <- 1
  
  dx <- arrow_scale * dx /
    norm_arrow
  
  dy <- arrow_scale * dy /
    norm_arrow
  
  dz <- arrow_scale * dz /
    norm_arrow
  
  p0 <- s3d$xyz.convert(
    x0,
    y0,
    z0
  )
  
  p1 <- s3d$xyz.convert(
    x0 + dx,
    y0 + dy,
    z0 + dz
  )
  
  arrows(
    x0 = p0$x,
    y0 = p0$y,
    x1 = p1$x,
    y1 = p1$y,
    length = arrow_head,
    angle = 30,
    col = arrow_col,
    lwd = arrow_lwd
  )
  
  invisible(
    data.frame(
      x0 = x0,
      y0 = y0,
      z0 = z0,
      x1 = x0 + dx,
      y1 = y0 + dy,
      z1 = z0 + dz
    )
  )
}

## ============================================================
## 3D plot
## ============================================================

plot_cell_fpc_fitted_path_3d <- function(res,
                                         traj_title = "Trajectory",
                                         cluster_order = NULL,
                                         cluster_cols = NULL,
                                         pc_use = c(1, 2, 3),
                                         point_cex = 1.5,
                                         point_alpha = 0.60,
                                         path_lwd = 5,
                                         arrow_n = 30,
                                         arrow_scale = 0.40,
                                         arrow_lwd = 4,
                                         arrow_head = 0.25,
                                         angle = view_angle,
                                         color_pseudotime = FALSE,
                                         grid = FALSE,
                                         box = FALSE,
                                         standardize_axes = TRUE,
                                         standard_axis_limit = 2,
                                         legend_cex = 1.5,
                                         axis_cex = 1,
                                         label_cex = 1.5,
                                         title_cex = 2,
                                         title_line = -1.0,
                                         cluster_legend_npc = c(0.73, 0.18),
                                         model_legend_npc = c(0.10, 0.18)) {
  
  proj <- get_cell_and_fitted_fpc_projection_3d(
    res = res,
    pc_use = pc_use
  )
  
  Z_obs <- proj$Z_obs
  Z_fit <- proj$Z_fit
  tt <- proj$tt
  clu_raw <- proj$cluster
  
  if (standardize_axes) {
    
    std_out <- standardize_fpc_coordinates_3d(
      Z_obs = Z_obs,
      Z_fit = Z_fit
    )
    
    Z_obs <- std_out$Z_obs
    Z_fit <- std_out$Z_fit
  }
  
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
  
  if (color_pseudotime) {
    
    tt_scaled <- (tt - min(tt, na.rm = TRUE)) /
      diff(range(tt, na.rm = TRUE))
    
    tt_scaled[!is.finite(tt_scaled)] <- 0
    
    pal <- colorRampPalette(
      c(
        "navy",
        "deepskyblue",
        "cyan",
        "yellow",
        "orange",
        "red3"
      )
    )(100)
    
    cell_cols <- pal[
      pmax(
        1,
        pmin(
          100,
          round(1 + 99 * tt_scaled)
        )
      )
    ]
    
  } else {
    
    point_cols <- grDevices::adjustcolor(
      unname(cluster_cols_use),
      alpha.f = point_alpha
    )
    
    names(point_cols) <- clu_levels
    
    cell_cols <- point_cols[
      as.character(clu)
    ]
  }
  
  if (anyNA(cell_cols)) {
    stop("Cell colors contain NA values")
  }
  
  if (standardize_axes) {
    
    axis_lims <- list(
      xlim = c(
        -standard_axis_limit,
        standard_axis_limit
      ),
      ylim = c(
        -standard_axis_limit,
        standard_axis_limit
      ),
      zlim = c(
        -standard_axis_limit,
        standard_axis_limit
      )
    )
    
    xlab_use <- paste0(
      "FPC",
      pc_use[1]
    )
    
    ylab_use <- paste0(
      "FPC",
      pc_use[2]
    )
    
    zlab_use <- paste0(
      "FPC",
      pc_use[3]
    )
    
  } else {
    
    axis_lims <- list(
      xlim = range(
        c(
          Z_obs[, 1],
          Z_fit[, 1]
        ),
        finite = TRUE
      ),
      ylim = range(
        c(
          Z_obs[, 2],
          Z_fit[, 2]
        ),
        finite = TRUE
      ),
      zlim = range(
        c(
          Z_obs[, 3],
          Z_fit[, 3]
        ),
        finite = TRUE
      )
    )
    
    xlab_use <- paste0(
      "Cell projection on FPC",
      pc_use[1]
    )
    
    ylab_use <- paste0(
      "Cell projection on FPC",
      pc_use[2]
    )
    
    zlab_use <- paste0(
      "Cell projection on FPC",
      pc_use[3]
    )
  }
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  par(
    mar = c(4.2, 4.4, 3.2, 1.4),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    col.sub = "black",
    cex.axis = axis_cex,
    cex.lab = label_cex,
    cex.main = title_cex,
    font.main = 2
  )
  
  s3d <- scatterplot3d(
    x = Z_obs[, 1],
    y = Z_obs[, 2],
    z = Z_obs[, 3],
    color = cell_cols,
    pch = 16,
    cex.symbols = point_cex,
    angle = angle,
    grid = grid,
    box = box,
    xlim = axis_lims$xlim,
    ylim = axis_lims$ylim,
    zlim = axis_lims$zlim,
    xlab = xlab_use,
    ylab = ylab_use,
    zlab = zlab_use,
    main = ""
  )
  
  title(
    main = paste0(
      traj_title,
      ": in FPC Space"
    ),
    line = title_line,
    cex.main = title_cex,
    font.main = 2
  )
  
  s3d$points3d(
    x = Z_fit[, 1],
    y = Z_fit[, 2],
    z = Z_fit[, 3],
    type = "l",
    col = "black",
    lwd = path_lwd
  )
  
  arrows_out <- add_3d_arrows(
    s3d = s3d,
    Z_fit = Z_fit,
    tt = tt,
    arrow_n = arrow_n,
    arrow_scale = arrow_scale,
    arrow_lwd = arrow_lwd,
    arrow_col = "black",
    arrow_head = arrow_head
  )
  
  if (!color_pseudotime) {
    
    draw_cluster_legend_npc(
      clu_levels = clu_levels,
      cluster_cols_use = cluster_cols_use,
      cluster_legend_npc = cluster_legend_npc,
      legend_cex = legend_cex,
      pt_cex = 1.4
    )
    
  } else {
    
    pal_cols <- c(
      "navy",
      "red3"
    )
    
    cluster_legend_xy <- npc_to_usr(
      cluster_legend_npc
    )
    
    legend(
      x = cluster_legend_xy["x"],
      y = cluster_legend_xy["y"],
      legend = c(
        "Low pseudotime",
        "High pseudotime"
      ),
      col = pal_cols,
      pch = 16,
      bty = "n",
      cex = legend_cex,
      pt.cex = 1.4,
      xpd = NA
    )
  }
  
  draw_model_legend_arrow(
    model_legend_npc = model_legend_npc,
    path_lwd = path_lwd,
    arrow_lwd = arrow_lwd,
    legend_cex = legend_cex,
    arrow_head = arrow_head
  )
  
  invisible(
    list(
      projection = proj,
      scatterplot3d = s3d,
      arrows = arrows_out,
      standardize_axes = standardize_axes,
      standard_axis_limit = standard_axis_limit
    )
  )
}

## ============================================================
## Plot in RStudio Plots pane
## ============================================================

cell_3d_traj1 <- plot_cell_fpc_fitted_path_3d(
  res = traj_results$traj1,
  traj_title = "Trajectory 1: Trans / Naive / C-mem1",
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1,
  pc_use = c(1, 2, 3),
  point_cex = 1,
  point_alpha = 0.65,
  path_lwd = 2.5,
  arrow_n = 30,
  arrow_scale = 0.30,
  arrow_lwd = 3,
  arrow_head = 0.15,
  angle = view_angle_traj1,
  color_pseudotime = FALSE,
  grid = F,
  box = FALSE,
  standardize_axes = TRUE,
  standard_axis_limit = standard_axis_limit,
  legend_cex = 1.1,
  axis_cex = 1,
  label_cex = 1.5,
  title_cex = 2,
  title_line = -1.0,
  cluster_legend_npc = c(0.65, 0.18),
  model_legend_npc = c(.050,.80)
)

cell_3d_traj2 <- plot_cell_fpc_fitted_path_3d(
  res = traj_results$traj2,
  traj_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2,
  pc_use = c(1, 2, 3),
  point_cex = 1,
  point_alpha = 0.65,
  path_lwd = 4.5,
  arrow_n = 30,
  arrow_scale = 0.30,
  arrow_lwd = 3,
  arrow_head = 0.15,
  angle = view_angle_traj2,
  color_pseudotime = FALSE,
  grid = FALSE,
  box = FALSE,
  standardize_axes = TRUE,
  standard_axis_limit = standard_axis_limit,
  legend_cex = 1.5,
  axis_cex = 1,
  label_cex = 1.5,
  title_cex = 2,
  title_line = -1.0,
  cluster_legend_npc =c(.75, .85),
  model_legend_npc = c(0.3, 0.13)
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
  
  plot_fun(...)
}

## ============================================================
## Save TIFF files
## ============================================================

save_tiff_plot(
  filename = "traj1_observed_cell_fpc123_fitted_path_3d.tiff",
  plot_fun = plot_cell_fpc_fitted_path_3d,
  width = panel_width,
  height = panel_height,
  dpi = 300,
  res = traj_results$traj1,
  traj_title = "Trajectory 1: Trans / Naive / C-mem1",
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1,
  pc_use = c(1, 2, 3),
  point_cex = 1.8,
  point_alpha = 0.6,
  path_lwd = 5,
  arrow_n = 30,
  arrow_scale = 0.40,
  arrow_lwd = 4,
  arrow_head = 0.25,
  angle = view_angle_traj1,
  color_pseudotime = FALSE,
  grid = FALSE,
  box = FALSE,
  standardize_axes = TRUE,
  standard_axis_limit = standard_axis_limit,
  legend_cex = 1.5,
  axis_cex = 2,
  label_cex = 2,
  title_cex = 1.5,
  title_line = -1.0,
  cluster_legend_npc = c(0.65, 0.18),
  model_legend_npc = c(.050,.80)
)

save_tiff_plot(
  filename = "traj2_observed_cell_fpc123_fitted_path_3d.tiff",
  plot_fun = plot_cell_fpc_fitted_path_3d,
  width = panel_width,
  height = panel_height,
  dpi = 300,
  res = traj_results$traj2,
  traj_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2,
  pc_use = c(1, 2, 3),
  point_cex = 1.8,
  point_alpha = 0.6,
  path_lwd = 5,
  arrow_n = 30,
  arrow_scale = 0.40,
  arrow_lwd = 4,
  arrow_head = 0.25,
  angle = view_angle_traj2,
  color_pseudotime = FALSE,
  grid = FALSE,
  box = FALSE,
  standardize_axes = TRUE,
  standard_axis_limit = standard_axis_limit,
  legend_cex = 1.5,
  axis_cex = 2,
  label_cex = 2,
  title_cex = 1.5,
  title_line = -1.0,
  cluster_legend_npc =c(.75, .85),
  model_legend_npc = c(0.3, 0.13)
)

## ============================================================
## Confirm saved files
## ============================================================

saved_files <- c(
  "traj1_observed_cell_fpc123_fitted_path_3d.tiff",
  "traj2_observed_cell_fpc123_fitted_path_3d.tiff"
)

print(
  normalizePath(
    saved_files,
    mustWork = FALSE
  )
)

print(
  file.exists(saved_files)
)
