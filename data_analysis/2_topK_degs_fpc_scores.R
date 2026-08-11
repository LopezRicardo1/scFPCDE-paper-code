## ============================================================
## Independent FPC score-space plots
##   scFPC-DE     : top K genes per quadrant
##   PseudotimeDE : top K_total genes by p-value
##
## Null boundary:
##   Always based on scFPC-DE non-TDEGs
##
## Visual transform:
##   Uses signed log10 transform:
##
##     sign(x) * log10(1 + abs(x))
##
##   This preserves quadrant signs while compressing extreme scores.
##
## PowerPoint:
##   - Final figure size: 30 x 24 inches
##   - Layout: 3 x 3
##   - Individual panel size: 10 x 8 inches
##
## Saves TIFF images directly in current working directory:
##   traj1_scfpcde_fpc_score_topK.tiff
##   traj1_pseudotimeDE_fpc_score_topK.tiff
##   traj2_scfpcde_fpc_score_topK.tiff
##   traj2_pseudotimeDE_fpc_score_topK.tiff
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

## ============================================================
## PowerPoint multipanel settings
## ============================================================

ppt_width <- 30
ppt_height <- 24

n_col <- 3
n_row <- 3

panel_width <- ppt_width / n_col
panel_height <- ppt_height / n_row

## ============================================================
## User-defined settings
## ============================================================

K_quad  <- 4
K_total <- 4 * K_quad

q_cut <- 0.05

ellipse_level <- 0.99

tiff_width <- panel_width
tiff_height <- panel_height
tiff_dpi <- 300

cex_pt <- 0.85
cex_lab <- 1.00
cex_corner <- 1.00

deg_col <- "black"
bg_col <- "gray75"

null_boundary_col <- "red3"
null_boundary_lwd <- 2.5
null_boundary_lty <- 2

pseudotime_inside_col <- "black"
pseudotime_inside_pch <- 1
pseudotime_inside_cex <- 1.7
pseudotime_inside_lwd <- 1.8

use_signed_log_scale <- TRUE

## ============================================================
## Signed log transform
## ============================================================

signed_log10 <- function(x) {
  
  sign(x) * log10(
    1 + abs(x)
  )
}

transform_score_matrix <- function(scores,
                                   use_signed_log = TRUE) {
  
  scores_out <- scores
  
  if (use_signed_log) {
    
    scores_out[, 1] <- signed_log10(
      scores_out[, 1]
    )
    
    scores_out[, 2] <- signed_log10(
      scores_out[, 2]
    )
  }
  
  scores_out
}

transform_null_object_for_plot <- function(null_object,
                                           use_signed_log = TRUE) {
  
  if (is.null(null_object)) {
    return(NULL)
  }
  
  null_plot <- null_object
  
  if (use_signed_log) {
    
    null_plot$ellipse[, 1] <- signed_log10(
      null_plot$ellipse[, 1]
    )
    
    null_plot$ellipse[, 2] <- signed_log10(
      null_plot$ellipse[, 2]
    )
  }
  
  null_plot
}

## ============================================================
## Safe score extractor
## ============================================================

get_fpca_scores <- function(fpca.res) {
  
  if (!is.null(fpca.res$fpca_result$scores)) {
    
    scores <- fpca.res$fpca_result$scores
    
  } else if (!is.null(fpca.res$fpca_result$fda_fpca$scores)) {
    
    scores <- fpca.res$fpca_result$fda_fpca$scores
    
  } else {
    
    stop(
      paste0(
        "No FPCA scores found in ",
        "fpca.res$fpca_result$scores or ",
        "fpca.res$fpca_result$fda_fpca$scores"
      )
    )
  }
  
  scores <- as.matrix(
    scores
  )
  
  storage.mode(scores) <- "numeric"
  
  if (ncol(scores) < 2) {
    stop("FPCA score matrix must contain at least two columns")
  }
  
  scores
}

## ============================================================
## Safe gene ID extractor
## ============================================================

get_fpca_gene_ids <- function(res,
                              scores) {
  
  gene_ids <- NULL
  
  if (
    !is.null(res$fpca.res$D_test_result$ID) &&
    length(res$fpca.res$D_test_result$ID) == nrow(scores)
  ) {
    
    gene_ids <- res$fpca.res$D_test_result$ID
    
  } else if (
    !is.null(rownames(scores)) &&
    length(rownames(scores)) == nrow(scores)
  ) {
    
    gene_ids <- rownames(scores)
  }
  
  if (is.null(gene_ids)) {
    stop("Gene IDs not found in D_test_result$ID or score rownames")
  }
  
  as.character(
    gene_ids
  )
}

## ============================================================
## Get scFPC-DE q-values
## ============================================================

get_scfpcde_qvalues <- function(res,
                                gene_ids) {
  
  dres <- as.data.frame(
    res$fpca.res$D_test_result
  )
  
  if (!"ID" %in% names(dres)) {
    stop("res$fpca.res$D_test_result must contain ID")
  }
  
  if ("q_value" %in% names(dres)) {
    
    qvals_all <- suppressWarnings(
      as.numeric(dres$q_value)
    )
    
  } else if ("p_value" %in% names(dres)) {
    
    qvals_all <- stats::p.adjust(
      suppressWarnings(
        as.numeric(dres$p_value)
      ),
      method = "BH"
    )
    
  } else {
    
    stop("D_test_result must contain q_value or p_value")
  }
  
  gene_map <- toupper(
    trimws(
      as.character(dres$ID)
    )
  )
  
  m <- match(
    toupper(trimws(gene_ids)),
    gene_map
  )
  
  qvals_all[m]
}

## ============================================================
## Get PseudotimeDE q-values
## ============================================================

get_pseudotimeDE_qvalues <- function(res,
                                     gene_ids) {
  
  if (is.null(res$res.gauss)) {
    stop("res$res.gauss was not found")
  }
  
  pt_df <- as.data.frame(
    res$res.gauss
  )
  
  if (!"gene" %in% names(pt_df)) {
    stop("res$res.gauss must contain gene")
  }
  
  if (!"fix.pv" %in% names(pt_df)) {
    stop("res$res.gauss must contain fix.pv")
  }
  
  pvals <- suppressWarnings(
    as.numeric(pt_df$fix.pv)
  )
  
  if (
    !is.null(res$gauss.qvals) &&
    length(res$gauss.qvals) == nrow(pt_df)
  ) {
    
    qvals_all <- suppressWarnings(
      as.numeric(res$gauss.qvals)
    )
    
  } else {
    
    qvals_all <- stats::p.adjust(
      pvals,
      method = "BH"
    )
  }
  
  gene_map <- toupper(
    trimws(
      as.character(pt_df$gene)
    )
  )
  
  m <- match(
    toupper(trimws(gene_ids)),
    gene_map
  )
  
  qvals_all[m]
}

## ============================================================
## Top scFPC-DE genes by quadrant
## ============================================================

get_top_scfpcde_by_quadrant <- function(res,
                                        top_k = K_quad,
                                        use_pcs = 1:3) {
  
  scores_all <- get_fpca_scores(
    res$fpca.res
  )
  
  gene_ids <- get_fpca_gene_ids(
    res = res,
    scores = scores_all
  )
  
  use_pcs <- use_pcs[
    use_pcs <= ncol(scores_all)
  ]
  
  scores_use <- scores_all[
    ,
    use_pcs,
    drop = FALSE
  ]
  
  quad_df <- data.frame(
    gene = as.character(gene_ids),
    FPC1 = scores_all[, 1],
    FPC2 = scores_all[, 2],
    dist_sq = rowSums(scores_use^2),
    stringsAsFactors = FALSE
  )
  
  quad_df$quadrant <- ifelse(
    quad_df$FPC1 > 0 & quad_df$FPC2 > 0,
    "Q1",
    ifelse(
      quad_df$FPC1 < 0 & quad_df$FPC2 > 0,
      "Q2",
      ifelse(
        quad_df$FPC1 < 0 & quad_df$FPC2 < 0,
        "Q3",
        "Q4"
      )
    )
  )
  
  quad_df %>%
    dplyr::filter(
      is.finite(FPC1),
      is.finite(FPC2),
      is.finite(dist_sq)
    ) %>%
    dplyr::group_by(
      quadrant
    ) %>%
    dplyr::arrange(
      dplyr::desc(dist_sq),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(
      n = top_k
    ) %>%
    dplyr::ungroup()
}

## ============================================================
## Top PseudotimeDE genes by p-value
## ============================================================

get_top_pseudotimeDE_by_pvalue <- function(res,
                                           top_k = K_total) {
  
  if (is.null(res$res.gauss)) {
    stop("res$res.gauss was not found")
  }
  
  if (!"gene" %in% names(res$res.gauss)) {
    stop("res$res.gauss must contain gene")
  }
  
  if (!"fix.pv" %in% names(res$res.gauss)) {
    stop("res$res.gauss must contain fix.pv")
  }
  
  res$res.gauss %>%
    dplyr::arrange(
      fix.pv
    ) %>%
    dplyr::pull(
      gene
    ) %>%
    unique() %>%
    head(
      top_k
    )
}

## ============================================================
## Estimate scFPC-DE null boundary on RAW FPC scores
## ============================================================

get_scfpcde_null_boundary <- function(res,
                                      q_cut = 0.05,
                                      ellipse_level = 0.99,
                                      n_points = 400) {
  
  scores <- get_fpca_scores(
    res$fpca.res
  )[
    ,
    1:2,
    drop = FALSE
  ]
  
  gene_ids <- get_fpca_gene_ids(
    res = res,
    scores = scores
  )
  
  qvals <- get_scfpcde_qvalues(
    res = res,
    gene_ids = gene_ids
  )
  
  is_scfpcde_deg <- is.finite(qvals) &
    qvals < q_cut
  
  X0 <- scores[
    !is_scfpcde_deg,
    ,
    drop = FALSE
  ]
  
  X0 <- X0[
    complete.cases(X0),
    ,
    drop = FALSE
  ]
  
  if (nrow(X0) < 5) {
    
    warning(
      "Too few scFPC-DE non-TDEGs to estimate the null boundary"
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
      "Invalid covariance matrix for scFPC-DE null boundary"
    )
    
    return(NULL)
  }
  
  trS <- sum(
    diag(S)
  )
  
  if (!is.finite(trS) || trS <= 0) {
    trS <- 1
  }
  
  S <- S +
    diag(
      1e-3 * trS,
      2
    )
  
  r2 <- stats::qchisq(
    ellipse_level,
    df = 2
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
  
  eig <- eigen(
    S,
    symmetric = TRUE
  )
  
  A <- eig$vectors %*%
    diag(
      sqrt(
        pmax(
          eig$values,
          0
        )
      ),
      2
    )
  
  ell <- t(
    matrix(
      mu,
      nrow = 2,
      ncol = n_points
    ) +
      sqrt(r2) * (A %*% circle)
  )
  
  colnames(ell) <- c(
    "FPC1",
    "FPC2"
  )
  
  list(
    ellipse = ell,
    center = mu,
    covariance = S,
    radius_sq = r2
  )
}

## ============================================================
## Flag points inside scFPC-DE null boundary using RAW scores
## ============================================================

flag_inside_scfpcde_null <- function(scores,
                                     null_object) {
  
  if (is.null(null_object)) {
    return(
      rep(FALSE, nrow(scores))
    )
  }
  
  X <- as.matrix(
    scores[
      ,
      1:2,
      drop = FALSE
    ]
  )
  
  d2 <- rep(
    NA_real_,
    nrow(X)
  )
  
  ok <- complete.cases(X)
  
  d2[ok] <- stats::mahalanobis(
    x = X[ok, , drop = FALSE],
    center = null_object$center,
    cov = null_object$covariance
  )
  
  is.finite(d2) &
    d2 <= null_object$radius_sq
}

## ============================================================
## Add corner labels
## ============================================================

add_corner_labels <- function(corner_labels = NULL,
                              corner_cols = NULL,
                              cex_corner = 1.00) {
  
  if (is.null(corner_labels)) {
    return(
      invisible(NULL)
    )
  }
  
  if (is.null(corner_cols)) {
    
    corner_cols <- list(
      tl = "darkgreen",
      tr = "blue3",
      bl = "orange3",
      br = "purple4"
    )
  }
  
  usr <- par("usr")
  
  xL <- usr[1]
  xR <- usr[2]
  yB <- usr[3]
  yT <- usr[4]
  
  if (!is.null(corner_labels$tl)) {
    
    text(
      x = xL + 0.025 * (xR - xL),
      y = yT - 0.030 * (yT - yB),
      labels = corner_labels$tl,
      col = corner_cols$tl,
      adj = c(0, 1),
      font = 2,
      cex = cex_corner
    )
  }
  
  if (!is.null(corner_labels$tr)) {
    
    text(
      x = xR - 0.025 * (xR - xL),
      y = yT - 0.030 * (yT - yB),
      labels = corner_labels$tr,
      col = corner_cols$tr,
      adj = c(1, 1),
      font = 2,
      cex = cex_corner
    )
  }
  
  if (!is.null(corner_labels$bl)) {
    
    text(
      x = xL + 0.025 * (xR - xL),
      y = yB + 0.035 * (yT - yB),
      labels = corner_labels$bl,
      col = corner_cols$bl,
      adj = c(0, 0),
      font = 2,
      cex = cex_corner
    )
  }
  
  if (!is.null(corner_labels$br)) {
    
    text(
      x = xR - 0.025 * (xR - xL),
      y = yB + 0.035 * (yT - yB),
      labels = corner_labels$br,
      col = corner_cols$br,
      adj = c(1, 0),
      font = 2,
      cex = cex_corner
    )
  }
  
  invisible(NULL)
}

## ============================================================
## Single independent FPC-score plot
## ============================================================

plot_fpc_gene_score_labelset <- function(res,
                                         label_genes,
                                         point_status = c(
                                           "scFPC-DE",
                                           "PseudotimeDE"
                                         ),
                                         q_cut = 0.05,
                                         xlim = NULL,
                                         ylim = NULL,
                                         main = NULL,
                                         cex_pt = 0.85,
                                         cex_lab = 1.00,
                                         corner_labels = NULL,
                                         corner_cols = NULL,
                                         cex_corner = 1.00,
                                         mark_pseudotime_inside_null = TRUE,
                                         use_signed_log = use_signed_log_scale,...) {
  
  point_status <- match.arg(
    point_status
  )
  
  scores_raw <- get_fpca_scores(
    res$fpca.res
  )[
    ,
    1:2,
    drop = FALSE
  ]
  
  scores_plot <- transform_score_matrix(
    scores = scores_raw,
    use_signed_log = use_signed_log
  )
  
  gene_ids <- get_fpca_gene_ids(
    res = res,
    scores = scores_raw
  )
  
  qvals_scfpcde <- get_scfpcde_qvalues(
    res = res,
    gene_ids = gene_ids
  )
  
  qvals_ptime <- get_pseudotimeDE_qvalues(
    res = res,
    gene_ids = gene_ids
  )
  
  is_scfpcde_sig <- is.finite(qvals_scfpcde) &
    qvals_scfpcde < q_cut
  
  is_ptime_sig <- is.finite(qvals_ptime) &
    qvals_ptime < q_cut
  
  is_sig <- if (point_status == "scFPC-DE") {
    is_scfpcde_sig
  } else {
    is_ptime_sig
  }
  
  point_col <- ifelse(
    is_sig,
    deg_col,
    bg_col
  )
  
  ## Null object stays raw for inside-null calculations
  null_object_raw <- get_scfpcde_null_boundary(
    res = res,
    q_cut = q_cut,
    ellipse_level = ellipse_level
  )
  
  inside_null <- flag_inside_scfpcde_null(
    scores = scores_raw,
    null_object = null_object_raw
  )
  
  ptime_sig_inside_null <- is_ptime_sig &
    inside_null
  
  ## Null ellipse transformed only for plotting
  null_object_plot <- transform_null_object_for_plot(
    null_object = null_object_raw,
    use_signed_log = use_signed_log
  )
  
  x_values <- scores_plot[, 1]
  y_values <- scores_plot[, 2]
  
  if (!is.null(null_object_plot)) {
    
    x_values <- c(
      x_values,
      null_object_plot$ellipse[, 1]
    )
    
    y_values <- c(
      y_values,
      null_object_plot$ellipse[, 2]
    )
  }
  
  if (is.null(xlim)) {
    xlim <- range(
      x_values,
      finite = TRUE
    )
  }
  
  if (is.null(ylim)) {
    ylim <- range(
      y_values,
      finite = TRUE
    )
  }
  
  xpad <- 0.08 * diff(
    xlim
  )
  
  ypad <- 0.08 * diff(
    ylim
  )
  
  if (!is.finite(xpad) || xpad <= 0) {
    xpad <- 0.5
  }
  
  if (!is.finite(ypad) || ypad <= 0) {
    ypad <- 0.5
  }
  
  xlim <- xlim + c(
    -xpad,
    xpad
  )
  
  ylim <- ylim + c(
    -ypad,
    ypad
  )
  
  xlab_use <- if (use_signed_log) {
    "Signed log10 FPC1 score"
  } else {
    "FPC1"
  }
  
  ylab_use <- if (use_signed_log) {
    "Signed log10 FPC2 score"
  } else {
    "FPC2"
  }
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par),
    add = TRUE
  )
  
  par(
    mar = c(6.8, 5.6, 4.8, 2.0),
    mgp = c(2.8, 0.85, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    col.sub = "black",
    xaxs = "r",
    yaxs = "r",
    xpd = FALSE,
    lwd = 1.5
  )
  
  plot(
    scores_plot[, 1],
    scores_plot[, 2],
    col = point_col,
    pch = 20,
    cex = cex_pt,
    xlab = xlab_use,
    ylab = ylab_use,
    main = main,
    xlim = xlim,
    ylim = ylim,
    cex.main = 1.35,
    cex.lab = 1.35,
    cex.axis = 1.10,...
  )
  
  abline(
    v = 0,
    h = 0,
    lty = 2,
    col = "gray45",
    lwd = 1.3
  )
  
  if (!is.null(null_object_plot)) {
    
    lines(
      null_object_plot$ellipse[, 1],
      null_object_plot$ellipse[, 2],
      col = null_boundary_col,
      lwd = null_boundary_lwd,
      lty = null_boundary_lty
    )
  }
  
  add_corner_labels(
    corner_labels = corner_labels,
    corner_cols = corner_cols,
    cex_corner = cex_corner
  )
  
  if (
    point_status == "PseudotimeDE" &&
    mark_pseudotime_inside_null &&
    any(ptime_sig_inside_null, na.rm = TRUE)
  ) {
    
    points(
      scores_plot[ptime_sig_inside_null, 1],
      scores_plot[ptime_sig_inside_null, 2],
      pch = pseudotime_inside_pch,
      cex = pseudotime_inside_cex,
      lwd = pseudotime_inside_lwd,
      col = pseudotime_inside_col
    )
  }
  
  label_genes <- as.character(
    label_genes
  )
  
  idx_lab <- match(
    toupper(label_genes),
    toupper(gene_ids)
  )
  
  idx_lab <- idx_lab[
    !is.na(idx_lab)
  ]
  
  if (length(idx_lab) > 0) {
    
    x <- scores_plot[idx_lab, 1]
    y <- scores_plot[idx_lab, 2]
    
    x_raw <- scores_raw[idx_lab, 1]
    y_raw <- scores_raw[idx_lab, 2]
    
    lab_col <- ifelse(
      x_raw > 0 & y_raw > 0,
      "blue3",
      ifelse(
        x_raw < 0 & y_raw > 0,
        "darkgreen",
        ifelse(
          x_raw < 0 & y_raw < 0,
          "orange3",
          "purple4"
        )
      )
    )
    
    y_offset <- 0.035 * diff(
      ylim
    )
    
    text(
      x,
      y + y_offset,
      labels = gene_ids[idx_lab],
      col = "white",
      font = 2,
      cex = cex_lab * 1.25
    )
    
    text(
      x,
      y + y_offset,
      labels = gene_ids[idx_lab],
      col = lab_col,
      font = 2,
      cex = cex_lab
    )
  }
  
  legend_items <- c(
    paste0(point_status, " TDEGs"),
    "Background genes",
    "scFPC-DE null boundary"
  )
  
  legend_cols <- c(
    deg_col,
    bg_col,
    null_boundary_col
  )
  
  legend_pch <- c(
    20,
    20,
    NA
  )
  
  legend_lty <- c(
    NA,
    NA,
    null_boundary_lty
  )
  
  legend_lwd <- c(
    NA,
    NA,
    null_boundary_lwd
  )
  
  if (point_status == "PseudotimeDE") {
    
    legend_items <- c(
      legend_items,
      "PseudotimeDE TDEG inside null"
    )
    
    legend_cols <- c(
      legend_cols,
      pseudotime_inside_col
    )
    
    legend_pch <- c(
      legend_pch,
      pseudotime_inside_pch
    )
    
    legend_lty <- c(
      legend_lty,
      NA
    )
    
    legend_lwd <- c(
      legend_lwd,
      pseudotime_inside_lwd
    )
  }
  
  legend(
    "bottom",
    inset = c(0, -0.32),
    xpd = NA,
    horiz = TRUE,
    bty = "n",
    legend = legend_items,
    col = legend_cols,
    pch = legend_pch,
    lty = legend_lty,
    lwd = legend_lwd,
    pt.cex = c(
      1.15,
      1.15,
      NA,
      1.25
    )[seq_along(legend_items)],
    cex = 0.95
  )
  
  invisible(
    list(
      gene_ids = gene_ids,
      scores_raw = scores_raw,
      scores_plot = scores_plot,
      is_scfpcde_sig = is_scfpcde_sig,
      is_ptime_sig = is_ptime_sig,
      inside_null = inside_null,
      ptime_sig_inside_null = ptime_sig_inside_null,
      null_object_raw = null_object_raw,
      null_object_plot = null_object_plot
    )
  )
}

## ============================================================
## Prepare top genes for one trajectory
## ============================================================

prepare_fpc_method_gene_sets <- function(res,
                                         top_k_quad = K_quad,
                                         top_k_ptime = K_total) {
  
  top_scfpcde_df <- get_top_scfpcde_by_quadrant(
    res = res,
    top_k = top_k_quad,
    use_pcs = 1:3
  )
  
  top_scfpcde <- top_scfpcde_df$gene
  
  top_ptime <- get_top_pseudotimeDE_by_pvalue(
    res = res,
    top_k = top_k_ptime
  )
  
  list(
    top_scfpcde = top_scfpcde,
    top_scfpcde_df = top_scfpcde_df,
    top_pseudotimeDE = top_ptime
  )
}

## ============================================================
## Quadrant profile labels
## ============================================================

corner_cols_default <- list(
  tl = "darkgreen",
  tr = "blue3",
  bl = "orange3",
  br = "purple4"
)

corner_traj1 <- list(
  tr = "Q1: FPC1+ / FPC2+\nlate-memory / terminal shift",
  tl = "Q2: FPC1- / FPC2+\nearly transitional",
  bl = "Q3: FPC1- / FPC2-\nintermediate / naive-like",
  br = "Q4: FPC1+ / FPC2-\nlate intermediate /\nactivation-shifted"
)

corner_traj2 <- list(
  tr = "Q1: FPC1+ / FPC2+\nearly-to-mid transition",
  tl = "Q2: FPC1- / FPC2+\nmid / DN3-shifted",
  bl = "Q3: FPC1- / FPC2-\nlate DN-like shift",
  br = "Q4: FPC1+ / FPC2-\nearly M-mem1-like"
)

## ============================================================
## Create gene sets
## ============================================================

compare_traj1 <- prepare_fpc_method_gene_sets(
  res = traj_results$traj1,
  top_k_quad = K_quad,
  top_k_ptime = K_total
)

compare_traj2 <- prepare_fpc_method_gene_sets(
  res = traj_results$traj2,
  top_k_quad = K_quad,
  top_k_ptime = K_total
)

topK_scfpcde_traj1 <- compare_traj1$top_scfpcde
topK_ptimeDE_traj1 <- compare_traj1$top_pseudotimeDE

topK_scfpcde_traj2 <- compare_traj2$top_scfpcde
topK_ptimeDE_traj2 <- compare_traj2$top_pseudotimeDE

## ============================================================
## Inspect selected genes
## ============================================================

print(topK_scfpcde_traj1)
print(topK_ptimeDE_traj1)

print(topK_scfpcde_traj2)
print(topK_ptimeDE_traj2)

print(compare_traj1$top_scfpcde_df)
print(compare_traj2$top_scfpcde_df)

## ============================================================
## Plot all four figures in RStudio Plots pane
## ============================================================

plot_fpc_gene_score_labelset(
  res = traj_results$traj1,
  label_genes = topK_scfpcde_traj1,
  point_status = "scFPC-DE",
  q_cut = q_cut,
  main = paste0(
    "Trajectory 1: Trans / Naive / C-mem1\n",
    "scFPC-DE top K per quadrant"
  ),
  corner_labels = corner_traj1,
  corner_cols = corner_cols_default,
  cex_pt = cex_pt,
  cex_lab = cex_lab,
  cex_corner = cex_corner,
  mark_pseudotime_inside_null = FALSE,
  use_signed_log = use_signed_log_scale, 
  ylim = c(-.35,.7),
  xlim = c(-.7,.3)
)

plot_fpc_gene_score_labelset(
  res = traj_results$traj1,
  label_genes = topK_ptimeDE_traj1,
  point_status = "PseudotimeDE",
  q_cut = q_cut,
  main = paste0(
    "Trajectory 1: Trans / Naive / C-mem1\n",
    "PseudotimeDE top genes by p-value"
  ),
  corner_labels = corner_traj1,
  corner_cols = corner_cols_default,
  cex_pt = cex_pt,
  cex_lab = cex_lab,
  cex_corner = cex_corner,
  mark_pseudotime_inside_null = TRUE,
  use_signed_log = use_signed_log_scale,
  ylim = c(-.35,.7),
  xlim = c(-.7,.3)
)

plot_fpc_gene_score_labelset(
  res = traj_results$traj2,
  label_genes = c(topK_scfpcde_traj2, "RHOB"),
  point_status = "scFPC-DE",
  q_cut = q_cut,
  main = paste0(
    "Trajectory 2: M-mem1 / DN3 / DN2\n",
    "scFPC-DE top K per quadrant"
  ),
  corner_labels = corner_traj2,
  corner_cols = corner_cols_default,
  cex_pt = cex_pt,
  cex_lab = cex_lab,
  cex_corner = cex_corner,
  mark_pseudotime_inside_null = FALSE,
  use_signed_log = use_signed_log_scale, 
  ylim = c(-.35,.55),
  xlim = c(-.4,.5)
)

plot_fpc_gene_score_labelset(
  res = traj_results$traj2,
  label_genes = topK_ptimeDE_traj2,
  point_status = "PseudotimeDE",
  q_cut = q_cut,
  main = paste0(
    "Trajectory 2: M-mem1 / DN3 / DN2\n",
    "PseudotimeDE top genes by p-value"
  ),
  corner_labels = corner_traj2,
  corner_cols = corner_cols_default,
  cex_pt = cex_pt,
  cex_lab = cex_lab,
  cex_corner = cex_corner,
  mark_pseudotime_inside_null = TRUE,
  use_signed_log = use_signed_log_scale
)

## ============================================================
## Save TIFF helper
## ============================================================

save_fpc_score_tiff <- function(filename,
                                res,
                                label_genes,
                                point_status,
                                main,
                                corner_labels,
                                corner_cols = corner_cols_default,
                                q_cut = 0.05,
                                xlim = NULL,
                                ylim = NULL,
                                width = tiff_width,
                                height = tiff_height,
                                dpi = tiff_dpi,
                                mark_pseudotime_inside_null = TRUE,
                                use_signed_log = use_signed_log_scale) {
  
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
  
  plot_fpc_gene_score_labelset(
    res = res,
    label_genes = label_genes,
    point_status = point_status,
    q_cut = q_cut,
    xlim = xlim,
    ylim = ylim,
    main = main,
    corner_labels = corner_labels,
    corner_cols = corner_cols,
    cex_pt = cex_pt,
    cex_lab = cex_lab,
    cex_corner = cex_corner,
    mark_pseudotime_inside_null = mark_pseudotime_inside_null,
    use_signed_log = use_signed_log
  )
  
  invisible(
    filename
  )
}

## ============================================================
## Save all four TIFF figures directly in current working dir
## ============================================================

save_fpc_score_tiff(
  filename = "traj1_scfpcde_fpc_score_topK.tiff",
  res = traj_results$traj1,
  label_genes = topK_scfpcde_traj1,
  point_status = "scFPC-DE",
  q_cut = q_cut,
  width = tiff_width,
  height = tiff_height,
  dpi = tiff_dpi,
  main = paste0(
    "Trajectory 1: Trans / Naive / C-mem1\n",
    "scFPC-DE top K per quadrant"
  ),
  corner_labels = corner_traj1,
  mark_pseudotime_inside_null = FALSE,
  use_signed_log = use_signed_log_scale,
  ylim = c(-.35,.7),
  xlim = c(-.7,.3)
)

save_fpc_score_tiff(
  filename = "traj1_pseudotimeDE_fpc_score_topK.tiff",
  res = traj_results$traj1,
  label_genes = topK_ptimeDE_traj1,
  point_status = "PseudotimeDE",
  q_cut = q_cut,
  width = tiff_width,
  height = tiff_height,
  dpi = tiff_dpi,
  main = paste0(
    "Trajectory 1: Trans / Naive / C-mem1\n",
    "PseudotimeDE top genes by p-value"
  ),
  corner_labels = corner_traj1,
  mark_pseudotime_inside_null = TRUE,
  use_signed_log = use_signed_log_scale
)

save_fpc_score_tiff(
  filename = "traj2_scfpcde_fpc_score_topK.tiff",
  res = traj_results$traj2,
  label_genes = topK_scfpcde_traj2,
  point_status = "scFPC-DE",
  q_cut = q_cut,
  width = tiff_width,
  height = tiff_height,
  dpi = tiff_dpi,
  main = paste0(
    "Trajectory 2: M-mem1 / DN3 / DN2\n",
    "scFPC-DE top K per quadrant"
  ),
  corner_labels = corner_traj2,
  mark_pseudotime_inside_null = FALSE,
  use_signed_log = use_signed_log_scale,
  ylim = c(-.35,.55),
  xlim = c(-.4,.5)
)

save_fpc_score_tiff(
  filename = "traj2_pseudotimeDE_fpc_score_topK.tiff",
  res = traj_results$traj2,
  label_genes = topK_ptimeDE_traj2,
  point_status = "PseudotimeDE",
  q_cut = q_cut,
  width = tiff_width,
  height = tiff_height,
  dpi = tiff_dpi,
  main = paste0(
    "Trajectory 2: M-mem1 / DN3 / DN2\n",
    "PseudotimeDE top genes by p-value"
  ),
  corner_labels = corner_traj2,
  mark_pseudotime_inside_null = TRUE,
  use_signed_log = use_signed_log_scale
)

## ============================================================
## Confirm saved files
## ============================================================

saved_fpc_topK_files <- c(
  "traj1_scfpcde_fpc_score_topK.tiff",
  "traj1_pseudotimeDE_fpc_score_topK.tiff",
  "traj2_scfpcde_fpc_score_topK.tiff",
  "traj2_pseudotimeDE_fpc_score_topK.tiff"
)

cat("\nCurrent working directory:\n")

print(
  getwd()
)

cat("\nSaved TIFF files:\n")

print(
  normalizePath(
    saved_fpc_topK_files,
    mustWork = FALSE
  )
)

cat("\nFile exists check:\n")

print(
  file.exists(saved_fpc_topK_files)
)
