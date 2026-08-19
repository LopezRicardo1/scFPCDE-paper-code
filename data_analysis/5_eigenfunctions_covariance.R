## ============================================================
## FPC eigenfunctions, PEV, and covariance contours
##
## For each trajectory:
##   1) Plot estimated FPC eigenfunctions
##   2) Include PEV in the eigenfunction legend
##   3) Add cell-type pseudotime tick tracks below the x-axis
##   4) Plot the PEV barplot
##   5) Reconstruct the covariance surface:
##
##        C(t, s) = sum_k lambda_k phi_k(t) phi_k(s)
##
##   6) Plot the covariance using filled.contour()
##
## PowerPoint:
##   - Final figure size: 30 x 24 inches
##   - Layout: 3 x 3
##   - Individual panel size: 10 x 8 inches
##
## Fixes:
##   - Different eigenfunction colors for trajectory 1 and trajectory 2
##   - Eigenfunction colors are not matched to cell-type colors
##   - No duplicated y-axis labels in eigenfunction plots
##   - Larger bottom margin so pseudotime tracks are not chopped
##
## Saves:
##   traj1_FPC_eigenfunctions.tiff
##   traj1_FPC_PEV.tiff
##   traj1_FPC_covariance_filled_contour.tiff
##   traj2_FPC_eigenfunctions.tiff
##   traj2_FPC_PEV.tiff
##   traj2_FPC_covariance_filled_contour.tiff
## ============================================================

suppressPackageStartupMessages({
  library(fda)
  library(plot3D)
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
## PowerPoint multipanel settings
## ------------------------------------------------------------

ppt_width <- 30
ppt_height <- 24

n_col <- 3
n_row <- 3

panel_width <- ppt_width / n_col
panel_height <- ppt_height / n_row

## ------------------------------------------------------------
## User-defined settings
## ------------------------------------------------------------

L_plot <- 3
cov_grid_n <- 150
cov_nlevels <- 20

## ------------------------------------------------------------
## FPC eigenfunction colors
## Different FPC color sets by trajectory.
## These are intentionally not matched to cell-type colors.
## ------------------------------------------------------------

eig_cols_traj1 <- c(
  "black",
  "gray35",
  "darkorange3"
)

eig_cols_traj2 <- c(
  "brown4",
  "navy",
  "darkolivegreen4"
)

## ------------------------------------------------------------
## Cell-type / cluster settings
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

cov_palette <- grDevices::colorRampPalette(
  plot3D::jet2.col()
)

## ============================================================
## Extract FPCA components
## ============================================================

extract_fpca_components <- function(res,
                                    L = 3,
                                    regular_grid_n = 150) {
  
  fpca_result <- res$fpca.res$fpca_result
  
  if (is.null(fpca_result)) {
    stop("res$fpca.res$fpca_result was not found")
  }
  
  ## ----------------------------------------------------------
  ## Original pseudotime grid
  ## ----------------------------------------------------------
  
  tt_original <- fpca_result$fda_splines$argvals
  
  if (is.null(tt_original)) {
    tt_original <- res$tt
  }
  
  if (is.null(tt_original)) {
    stop("The pseudotime grid was not found")
  }
  
  if (is.matrix(tt_original) || is.data.frame(tt_original)) {
    tt_original <- tt_original[, 1]
  }
  
  tt_original <- as.numeric(
    tt_original
  )
  
  if (any(!is.finite(tt_original))) {
    stop("The pseudotime grid contains non-finite values")
  }
  
  ord <- order(
    tt_original
  )
  
  tt_original_sorted <- tt_original[
    ord
  ]
  
  tt_regular <- seq(
    min(tt_original_sorted),
    max(tt_original_sorted),
    length.out = regular_grid_n
  )
  
  ## ----------------------------------------------------------
  ## Locate internal FPCA object
  ## ----------------------------------------------------------
  
  fpca_fit <- NULL
  
  if (!is.null(fpca_result$fda_fpca)) {
    
    fpca_fit <- fpca_result$fda_fpca
    
  } else if (!is.null(fpca_result$fpca)) {
    
    fpca_fit <- fpca_result$fpca
    
  } else if (!is.null(fpca_result$pca)) {
    
    fpca_fit <- fpca_result$pca
  }
  
  ## ----------------------------------------------------------
  ## Extract eigenfunctions
  ## ----------------------------------------------------------
  
  harmonics <- NULL
  eigenfunctions_original <- NULL
  
  if (!is.null(fpca_result$harmonics)) {
    
    harmonics <- fpca_result$harmonics
    
  } else if (!is.null(fpca_fit$harmonics)) {
    
    harmonics <- fpca_fit$harmonics
    
  } else if (!is.null(fpca_result$eigenfunctions)) {
    
    eigenfunctions_original <- fpca_result$eigenfunctions
    
  } else if (!is.null(fpca_fit$eigenfunctions)) {
    
    eigenfunctions_original <- fpca_fit$eigenfunctions
    
  } else {
    
    stop("The FPC eigenfunctions were not found")
  }
  
  if (!is.null(harmonics) && inherits(harmonics, "fd")) {
    
    eigenfunctions_regular <- fda::eval.fd(
      evalarg = tt_regular,
      fdobj = harmonics
    )
    
  } else {
    
    if (!is.null(harmonics)) {
      eigenfunctions_original <- harmonics
    }
    
    eigenfunctions_original <- as.matrix(
      eigenfunctions_original
    )
    
    if (
      nrow(eigenfunctions_original) != length(tt_original) &&
      ncol(eigenfunctions_original) == length(tt_original)
    ) {
      
      eigenfunctions_original <- t(
        eigenfunctions_original
      )
    }
    
    if (nrow(eigenfunctions_original) != length(tt_original)) {
      stop("Eigenfunction rows do not match the original pseudotime-grid length")
    }
    
    eigenfunctions_original <- eigenfunctions_original[
      ord,
      ,
      drop = FALSE
    ]
    
    tt_unique <- sort(
      unique(tt_original_sorted)
    )
    
    eigenfunctions_unique <- vapply(
      seq_len(ncol(eigenfunctions_original)),
      function(k) {
        
        vapply(
          tt_unique,
          function(tt_i) {
            
            mean(
              eigenfunctions_original[
                tt_original_sorted == tt_i,
                k
              ],
              na.rm = TRUE
            )
          },
          FUN.VALUE = numeric(1)
        )
      },
      FUN.VALUE = numeric(length(tt_unique))
    )
    
    eigenfunctions_unique <- as.matrix(
      eigenfunctions_unique
    )
    
    eigenfunctions_regular <- vapply(
      seq_len(ncol(eigenfunctions_unique)),
      function(k) {
        
        stats::approx(
          x = tt_unique,
          y = eigenfunctions_unique[, k],
          xout = tt_regular,
          rule = 2,
          ties = mean
        )$y
      },
      FUN.VALUE = numeric(length(tt_regular))
    )
    
    eigenfunctions_regular <- as.matrix(
      eigenfunctions_regular
    )
  }
  
  ## ----------------------------------------------------------
  ## Extract eigenvalues
  ## ----------------------------------------------------------
  
  eigenvalues <- NULL
  
  if (!is.null(fpca_result$eigenvalues)) {
    
    eigenvalues <- fpca_result$eigenvalues
    
  } else if (!is.null(fpca_result$values)) {
    
    eigenvalues <- fpca_result$values
    
  } else if (!is.null(fpca_fit$eigenvalues)) {
    
    eigenvalues <- fpca_fit$eigenvalues
    
  } else if (!is.null(fpca_fit$values)) {
    
    eigenvalues <- fpca_fit$values
  }
  
  if (!is.null(eigenvalues)) {
    
    eigenvalues <- as.numeric(
      eigenvalues
    )
  }
  
  ## ----------------------------------------------------------
  ## Extract stored PEV
  ## ----------------------------------------------------------
  
  varprop <- NULL
  
  if (!is.null(fpca_result$varprop)) {
    
    varprop <- as.numeric(
      fpca_result$varprop
    )
    
  } else if (!is.null(fpca_fit$varprop)) {
    
    varprop <- as.numeric(
      fpca_fit$varprop
    )
  }
  
  if (!is.null(varprop)) {
    
    varprop[!is.finite(varprop)] <- 0
    
    if (sum(varprop, na.rm = TRUE) > 1.5) {
      varprop <- varprop / 100
    }
  }
  
  ## ----------------------------------------------------------
  ## Estimate eigenvalues from gene FPC scores if necessary
  ## ----------------------------------------------------------
  
  if (is.null(eigenvalues)) {
    
    scores <- fpca_result$scores
    
    if (is.null(scores) && !is.null(fpca_fit$scores)) {
      scores <- fpca_fit$scores
    }
    
    if (!is.null(scores)) {
      
      scores <- as.matrix(
        scores
      )
      
      eigenvalues <- apply(
        scores,
        MARGIN = 2,
        FUN = stats::var,
        na.rm = TRUE
      )
    }
  }
  
  ## ----------------------------------------------------------
  ## Final covariance-weight fallback
  ## ----------------------------------------------------------
  
  normalized_covariance <- FALSE
  
  if (is.null(eigenvalues)) {
    
    if (is.null(varprop)) {
      stop("Neither FPCA eigenvalues, FPC scores, nor PEV were found")
    }
    
    eigenvalues <- varprop
    normalized_covariance <- TRUE
  }
  
  eigenvalues <- as.numeric(
    eigenvalues
  )
  
  eigenvalues[!is.finite(eigenvalues)] <- 0
  eigenvalues[eigenvalues < 0] <- 0
  
  if (sum(eigenvalues) <= 0) {
    stop("The extracted FPCA eigenvalues are invalid")
  }
  
  if (!is.null(varprop) && sum(varprop) > 0) {
    
    pev_all <- varprop / sum(varprop)
    
  } else {
    
    pev_all <- eigenvalues / sum(eigenvalues)
  }
  
  L_use <- min(
    L,
    ncol(eigenfunctions_regular),
    length(eigenvalues),
    length(pev_all)
  )
  
  if (L_use < 1) {
    stop("No valid FPC components were found")
  }
  
  eigenfunctions_use <- eigenfunctions_regular[
    ,
    seq_len(L_use),
    drop = FALSE
  ]
  
  eigenvalues_use <- eigenvalues[
    seq_len(L_use)
  ]
  
  pev_use <- pev_all[
    seq_len(L_use)
  ]
  
  colnames(eigenfunctions_use) <- paste0(
    "FPC",
    seq_len(L_use)
  )
  
  ## ----------------------------------------------------------
  ## Covariance reconstruction
  ## ----------------------------------------------------------
  
  covariance_surface <- eigenfunctions_use %*%
    diag(
      eigenvalues_use,
      nrow = L_use,
      ncol = L_use
    ) %*%
    t(eigenfunctions_use)
  
  covariance_surface <- 0.5 * (
    covariance_surface + t(covariance_surface)
  )
  
  list(
    tt = tt_regular,
    eigenfunctions = eigenfunctions_use,
    eigenvalues = eigenvalues_use,
    pev = pev_use,
    pev_all = pev_all,
    covariance = covariance_surface,
    L = L_use,
    normalized_covariance = normalized_covariance
  )
}

## ============================================================
## Extract cell-level pseudotime and cell-type information
## ============================================================

extract_celltype_pseudotime <- function(res,
                                        cluster_order,
                                        cluster_cols) {
  
  tt_cells <- as.numeric(
    res$tt
  )
  
  cluster_cells <- get_analysis_clusters(res)
  
  cluster_cells <- as.character(
    cluster_cells
  )
  
  if (length(tt_cells) != length(cluster_cells)) {
    stop("Cell-level pseudotime and cluster vectors have different lengths")
  }
  
  keep <- is.finite(tt_cells) &
    !is.na(cluster_cells)
  
  tt_cells <- tt_cells[
    keep
  ]
  
  cluster_cells <- cluster_cells[
    keep
  ]
  
  cluster_order_use <- cluster_order[
    cluster_order %in% unique(cluster_cells)
  ]
  
  if (length(cluster_order_use) == 0) {
    stop("None of the requested clusters were found")
  }
  
  cluster_factor <- factor(
    cluster_cells,
    levels = cluster_order_use
  )
  
  if (anyNA(cluster_factor)) {
    
    missing_clusters <- unique(
      cluster_cells[is.na(cluster_factor)]
    )
    
    stop(
      "Clusters missing from cluster_order: ",
      paste(missing_clusters, collapse = ", ")
    )
  }
  
  if (is.null(names(cluster_cols))) {
    
    cluster_cols <- setNames(
      cluster_cols[
        seq_along(cluster_order_use)
      ],
      cluster_order_use
    )
  }
  
  missing_cols <- setdiff(
    cluster_order_use,
    names(cluster_cols)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      "Missing colors for clusters: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  list(
    tt = tt_cells,
    cluster = cluster_factor,
    cluster_order = cluster_order_use,
    cluster_cols = cluster_cols[
      cluster_order_use
    ]
  )
}

## ============================================================
## Eigenfunction plot with PEV and cell-type pseudotime tracks
## ============================================================

plot_fpca_eigenfunctions <- function(res,
                                     trajectory_title,
                                     L = L_plot,
                                     line_cols,
                                     cluster_order,
                                     cluster_cols,
                                     line_lwd = 4.0,
                                     regular_grid_n = cov_grid_n,
                                     tick_lwd = 1.10,
                                     tick_height = 0.018,
                                     track_gap = 0.052,
                                     upper_y_expansion = 0.26,
                                     lower_y_expansion = 0.04,
                                     rug_offset = 0.145,
                                     title_offset = 0.085,
                                     title_cex = 1.60,
                                     axis_cex = 1.15,
                                     label_cex = 1.35,
                                     legend_cex = 1.15,
                                     track_label_cex = 1.05,
                                     x_label_cex = 1.25) {
  
  fpca_comp <- extract_fpca_components(
    res = res,
    L = L,
    regular_grid_n = regular_grid_n
  )
  
  tt_grid <- fpca_comp$tt
  phi <- fpca_comp$eigenfunctions
  L_use <- fpca_comp$L
  
  pev <- res$fpca.res.deg$fpca_result$fda_fpca$varprop
  
  if (is.null(pev)) {
    stop("PEV/varprop not found at res$fpca.res.deg$fpca_result$fda_fpca$varprop")
  }
  
  pev <- as.numeric(
    pev
  )
  
  pev <- pev[
    is.finite(pev)
  ]
  
  if (length(pev) == 0) {
    stop("PEV/varprop vector is empty")
  }
  
  if (max(pev, na.rm = TRUE) > 1) {
    pev <- pev / 100
  }
  
  L_use <- min(
    L_use,
    length(pev),
    ncol(phi)
  )
  
  phi <- phi[
    ,
    seq_len(L_use),
    drop = FALSE
  ]
  
  pev <- pev[
    seq_len(L_use)
  ]
  
  line_cols <- rep(
    line_cols,
    length.out = L_use
  )
  
  cell_info <- extract_celltype_pseudotime(
    res = res,
    cluster_order = cluster_order,
    cluster_cols = cluster_cols
  )
  
  tt_cells <- cell_info$tt
  cluster_factor <- cell_info$cluster
  cluster_order_use <- cell_info$cluster_order
  cluster_cols_use <- cell_info$cluster_cols
  
  phi_range <- range(
    phi,
    finite = TRUE
  )
  
  phi_span <- diff(
    phi_range
  )
  
  if (!is.finite(phi_span) || phi_span <= 0) {
    phi_span <- 1
  }
  
  ylim_use <- c(
    phi_range[1] - lower_y_expansion * phi_span,
    phi_range[2] + upper_y_expansion * phi_span
  )
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  par(
    mar = c(10.2, 5.8, 4.8, 2.0),
    mgp = c(3.2, 0.95, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    xpd = FALSE,
    lwd = 1.5
  )
  
  matplot(
    x = tt_grid,
    y = phi,
    type = "l",
    lty = 1,
    lwd = line_lwd,
    col = line_cols,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = expression(paste("Eigenfunctions ", phi(t))),
    main = paste0(
      trajectory_title,
      "\nFPC eigenfunctions"
    ),
    ylim = ylim_use,
    cex.main = title_cex,
    cex.lab = label_cex,
    cex.axis = axis_cex
  )
  
  abline(
    h = 0,
    lty = 2,
    lwd = 1.5,
    col = "gray60"
  )
  
  axis(
    side = 1,
    line = 0,
    tck = -0.012,
    col = "black",
    col.axis = "black",
    las = 1,
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
    lwd = 1.5
  )
  
  usr <- par("usr")
  
  x_range <- usr[2] - usr[1]
  y_range <- usr[4] - usr[3]
  
  track_height_y <- tick_height * y_range
  track_gap_y <- track_gap * y_range
  
  first_track_y <- usr[3] -
    rug_offset * y_range
  
  track_y <- first_track_y -
    track_gap_y *
    (
      seq_along(cluster_order_use) - 1
    )
  
  par(
    xpd = NA
  )
  
  for (j in seq_along(cluster_order_use)) {
    
    cluster_name <- cluster_order_use[
      j
    ]
    
    tt_cluster <- tt_cells[
      cluster_factor == cluster_name
    ]
    
    segments(
      x0 = tt_cluster,
      y0 = track_y[j] - track_height_y / 2,
      x1 = tt_cluster,
      y1 = track_y[j] + track_height_y / 2,
      col = unname(
        cluster_cols_use[cluster_name]
      ),
      lwd = tick_lwd
    )
    
    text(
      x = usr[1] - 0.035 * x_range,
      y = track_y[j],
      labels = cluster_name,
      adj = c(1, 0.5),
      col = unname(
        cluster_cols_use[cluster_name]
      ),
      cex = track_label_cex,
      font = 2
    )
  }
  
  title_y <- min(track_y) -
    title_offset * y_range
  
  text(
    x = mean(usr[1:2]),
    y = title_y,
    labels = "Pseudotime (t)",
    cex = x_label_cex
  )
  
  par(
    xpd = FALSE
  )
  
  legend_labels <- paste0(
    "FPC",
    seq_len(L_use),
    " (PEV = ",
    sprintf(
      "%.1f",
      100 * pev
    ),
    "%)"
  )
  
  legend(
    "topleft",
    inset = c(0.012, 0.018),
    legend = legend_labels,
    col = line_cols,
    lty = 1,
    lwd = line_lwd,
    bty = "n",
    cex = legend_cex,
    seg.len = 2.8,
    y.intersp = 1.10
  )
  
  invisible(
    list(
      fpca_components = fpca_comp,
      pev = pev,
      pev_percent = 100 * pev,
      cumulative_pev_percent = cumsum(100 * pev),
      top3_pev_percent = sum(
        100 * pev[
          seq_len(
            min(3, length(pev))
          )
        ]
      ),
      track_y = track_y,
      cluster_order = cluster_order_use,
      cluster_cols = cluster_cols_use,
      ylim = ylim_use
    )
  )
}

## ============================================================
## PEV barplot
## ============================================================

plot_fpca_pev <- function(res,
                          trajectory_title,
                          bar_cols,
                          L = L_plot,
                          regular_grid_n = cov_grid_n,
                          title_cex = 1.60,
                          axis_cex = 1.20,
                          label_cex = 1.35,
                          value_cex = 1.05) {
  
  fpca_comp <- extract_fpca_components(
    res = res,
    L = L,
    regular_grid_n = regular_grid_n
  )
  
  pev <- fpca_comp$pev
  L_use <- fpca_comp$L
  
  bar_cols <- rep(
    bar_cols,
    length.out = L_use
  )
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  par(
    mar = c(5.8, 5.8, 4.8, 2.0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    lwd = 1.5
  )
  
  bp <- barplot(
    height = 100 * pev,
    names.arg = paste0(
      "FPC",
      seq_len(L_use)
    ),
    col = bar_cols,
    border = "black",
    ylim = c(
      0,
      max(100 * pev, na.rm = TRUE) * 1.22
    ),
    ylab = "Percent explained variance",
    main = paste0(
      trajectory_title,
      "\nFPC percent explained variance"
    ),
    cex.main = title_cex,
    cex.lab = label_cex,
    cex.axis = axis_cex,
    cex.names = axis_cex,
    las = 1
  )
  
  text(
    x = bp,
    y = 100 * pev,
    labels = paste0(
      sprintf(
        "%.1f",
        100 * pev
      ),
      "%"
    ),
    pos = 3,
    cex = value_cex
  )
  
  box(
    lwd = 1.5
  )
  
  invisible(
    fpca_comp
  )
}

## ============================================================
## Filled covariance contour
## ============================================================

plot_fpca_covariance_filled <- function(res,
                                        trajectory_title,
                                        L = L_plot,
                                        regular_grid_n = cov_grid_n,
                                        nlevels = cov_nlevels,
                                        color_palette = cov_palette,
                                        title_cex = 1.60,
                                        axis_cex = 1.15,
                                        label_cex = 1.35,
                                        key_cex = 0.95) {
  
  fpca_comp <- extract_fpca_components(
    res = res,
    L = L,
    regular_grid_n = regular_grid_n
  )
  
  grid_t1 <- fpca_comp$tt
  grid_t2 <- fpca_comp$tt
  cov_surface <- fpca_comp$covariance
  
  covariance_label <- if (
    fpca_comp$normalized_covariance
  ) {
    
    "Normalized covariance"
    
  } else {
    
    "Covariance"
  }
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  par(
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black"
  )
  
  filled.contour(
    x = grid_t1,
    y = grid_t2,
    z = cov_surface,
    color.palette = color_palette,
    nlevels = nlevels,
    xlab = "Pseudotime (t)",
    ylab = "Pseudotime (s)",
    main = paste0(
      trajectory_title,
      "\nCovariance surface"
    ),
    cex.main = title_cex,
    cex.lab = label_cex,
    cex.axis = axis_cex,
    key.title = title(
      main = covariance_label,
      cex.main = key_cex
    ),
    key.axes = axis(
      side = 4,
      las = 1,
      cex.axis = axis_cex
    ),
    plot.axes = {
      
      axis(
        side = 1,
        cex.axis = axis_cex
      )
      
      axis(
        side = 2,
        las = 1,
        cex.axis = axis_cex
      )
      
      contour(
        x = grid_t1,
        y = grid_t2,
        z = cov_surface,
        nlevels = 10,
        add = TRUE,
        drawlabels = FALSE,
        col = grDevices::adjustcolor(
          "white",
          alpha.f = 0.70
        ),
        lwd = 0.9
      )
      
      box(
        lwd = 1.5
      )
    }
  )
  
  invisible(
    fpca_comp
  )
}

## ============================================================
## Extract FPCA results and PEV tables
## ============================================================

fpca_components_traj1 <- extract_fpca_components(
  res = traj_results$traj1,
  L = L_plot,
  regular_grid_n = cov_grid_n
)

fpca_components_traj2 <- extract_fpca_components(
  res = traj_results$traj2,
  L = L_plot,
  regular_grid_n = cov_grid_n
)

PEV_table_traj1 <- data.frame(
  FPC = paste0(
    "FPC",
    seq_len(fpca_components_traj1$L)
  ),
  Eigenvalue = fpca_components_traj1$eigenvalues,
  PEV = fpca_components_traj1$pev,
  PEV_percent = 100 * fpca_components_traj1$pev,
  Cumulative_PEV = cumsum(
    fpca_components_traj1$pev
  ),
  Cumulative_PEV_percent = 100 * cumsum(
    fpca_components_traj1$pev
  )
)

PEV_table_traj2 <- data.frame(
  FPC = paste0(
    "FPC",
    seq_len(fpca_components_traj2$L)
  ),
  Eigenvalue = fpca_components_traj2$eigenvalues,
  PEV = fpca_components_traj2$pev,
  PEV_percent = 100 * fpca_components_traj2$pev,
  Cumulative_PEV = cumsum(
    fpca_components_traj2$pev
  ),
  Cumulative_PEV_percent = 100 * cumsum(
    fpca_components_traj2$pev
  )
)

print(
  PEV_table_traj1
)

print(
  PEV_table_traj2
)

## ============================================================
## Display plots in RStudio Plots pane
## ============================================================

plot_fpca_eigenfunctions(
  res = traj_results$traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  line_cols = eig_cols_traj1,
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1,
  L = L_plot,
  regular_grid_n = cov_grid_n
)

plot_fpca_pev(
  res = traj_results$traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  bar_cols = eig_cols_traj1,
  L = L_plot,
  regular_grid_n = cov_grid_n
)

plot_fpca_covariance_filled(
  res = traj_results$traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  L = L_plot,
  regular_grid_n = cov_grid_n,
  nlevels = cov_nlevels,
  color_palette = cov_palette
)

plot_fpca_eigenfunctions(
  res = traj_results$traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  line_cols = eig_cols_traj2,
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2,
  L = L_plot,
  regular_grid_n = cov_grid_n
)

plot_fpca_pev(
  res = traj_results$traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  bar_cols = eig_cols_traj2,
  L = L_plot,
  regular_grid_n = cov_grid_n
)

plot_fpca_covariance_filled(
  res = traj_results$traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  L = L_plot,
  regular_grid_n = cov_grid_n,
  nlevels = cov_nlevels,
  color_palette = cov_palette
)

## ============================================================
## TIFF-saving functions
## ============================================================

save_eigenfunction_tiff <- function(filename,
                                    res,
                                    trajectory_title,
                                    line_cols,
                                    cluster_order,
                                    cluster_cols,
                                    width = panel_width,
                                    height = panel_height,
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
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  plot_fpca_eigenfunctions(
    res = res,
    trajectory_title = trajectory_title,
    line_cols = line_cols,
    cluster_order = cluster_order,
    cluster_cols = cluster_cols,
    L = L_plot,
    regular_grid_n = cov_grid_n
  )
}

save_pev_tiff <- function(filename,
                          res,
                          trajectory_title,
                          bar_cols,
                          width = panel_width,
                          height = panel_height,
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
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  plot_fpca_pev(
    res = res,
    trajectory_title = trajectory_title,
    bar_cols = bar_cols,
    L = L_plot,
    regular_grid_n = cov_grid_n
  )
}

save_covariance_tiff <- function(filename,
                                 res,
                                 trajectory_title,
                                 width = panel_width,
                                 height = panel_height,
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
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  plot_fpca_covariance_filled(
    res = res,
    trajectory_title = trajectory_title,
    L = L_plot,
    regular_grid_n = cov_grid_n,
    nlevels = cov_nlevels,
    color_palette = cov_palette
  )
}

## ============================================================
## Save all TIFF figures
## ============================================================

saved_files <- c(
  "traj1_FPC_eigenfunctions.tiff",
  "traj1_FPC_PEV.tiff",
  "traj1_FPC_covariance_filled_contour.tiff",
  "traj2_FPC_eigenfunctions.tiff",
  "traj2_FPC_PEV.tiff",
  "traj2_FPC_covariance_filled_contour.tiff"
)

save_eigenfunction_tiff(
  filename = saved_files[1],
  res = traj_results$traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  line_cols = eig_cols_traj1,
  cluster_order = cluster_order_traj1,
  cluster_cols = cluster_cols_traj1
)

save_pev_tiff(
  filename = saved_files[2],
  res = traj_results$traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  bar_cols = eig_cols_traj1
)

save_covariance_tiff(
  filename = saved_files[3],
  res = traj_results$traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1"
)

save_eigenfunction_tiff(
  filename = saved_files[4],
  res = traj_results$traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  line_cols = eig_cols_traj2,
  cluster_order = cluster_order_traj2,
  cluster_cols = cluster_cols_traj2
)

save_pev_tiff(
  filename = saved_files[5],
  res = traj_results$traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  bar_cols = eig_cols_traj2
)

save_covariance_tiff(
  filename = saved_files[6],
  res = traj_results$traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2"
)

## ============================================================
## Confirm saved files
## ============================================================

cat("Saved TIFF files:\n")

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
