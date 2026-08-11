## ============================================================
## Full code: Euler TDEG and Reactome plots saved with grid.arrange()
##
## TDEG calculation:
##   - scFPC-DE: BH adjustment of D_test_result$p_value
##   - PseudotimeDE: BH adjustment of res.gauss$fix.pv
##   - Stored q_value and gauss.qvals objects are not used
##   - Significant TDEGs: BH-adjusted p-value < 0.05
##
## Saves directly in current working directory:
##   traj1_DEG_overlap_euler.tiff
##   traj2_DEG_overlap_euler.tiff
##   traj1_Reactome_pathway_overlap_euler.tiff
##   traj2_Reactome_pathway_overlap_euler.tiff
##   DEG_and_Reactome_overlap_euler_2x2.tiff
##   DEG_overlap_summary.csv
##   Reactome_pathway_overlap_summary.csv
##
## Required objects:
##   traj_results$traj1
##   traj_results$traj2
##   reactome_ora_results$traj1
##   reactome_ora_results$traj2
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(eulerr)
  library(grid)
  library(gridExtra)
})

graphics.off()

## ============================================================
## Settings
## ============================================================

fdr_cut <- 0.05
ora_padj_cut <- 0.05

euler_fill_cols <- c(
  "#8ECAE6",
  "#F4A6B8"
)

euler_alpha <- 0.75

single_width <- 7
single_height <- 6

combined_width <- 11
combined_height <- 12

tiff_dpi <- 300

euler_label_cex <- 0
euler_count_cex <- 3
euler_edge_lwd <- 2

cat("\nCurrent working directory:\n")
print(getwd())

if (file.access(getwd(), 2) != 0) {
  stop("Current working directory is not writable: ", getwd())
}

## ============================================================
## BH adjustment helper
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
## Extract complete scFPC-DE result table
##
## BH-adjusted p-values are recomputed from:
##   D_test_result$p_value
## ============================================================

get_scfpcde_fdr_table <- function(res) {
  
  if (
    is.null(res$fpca.res) ||
    is.null(res$fpca.res$D_test_result)
  ) {
    stop("res$fpca.res$D_test_result was not found")
  }
  
  dres <- as.data.frame(
    res$fpca.res$D_test_result
  )
  
  required_columns <- c(
    "ID",
    "p_value"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(dres)
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
    as.numeric(dres$p_value)
  )
  
  qvalues <- bh_adjust_finite(
    pvalues
  )
  
  if ("D_obs" %in% names(dres)) {
    
    D_obs <- suppressWarnings(
      as.numeric(dres$D_obs)
    )
    
  } else {
    
    D_obs <- rep(
      NA_real_,
      nrow(dres)
    )
  }
  
  out <- data.frame(
    Gene_ID = toupper(
      trimws(
        as.character(dres$ID)
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
## Extract complete PseudotimeDE result table
##
## BH-adjusted p-values are recomputed from:
##   res.gauss$fix.pv
## ============================================================

get_pseudotimeDE_fdr_table <- function(res) {
  
  if (is.null(res$res.gauss)) {
    stop("res$res.gauss was not found")
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
## Extract scFPC-DE TDEG set
## ============================================================

get_scfpcde_deg_set <- function(res,
                                q_cut = 0.05) {
  
  sc_tbl <- get_scfpcde_fdr_table(
    res = res
  )
  
  sc_tbl %>%
    dplyr::filter(
      is.finite(q_value),
      q_value < q_cut
    ) %>%
    dplyr::pull(
      Gene_ID
    ) %>%
    unique()
}

## ============================================================
## Extract PseudotimeDE TDEG set
## ============================================================

get_pseudotimeDE_deg_set <- function(res,
                                     q_cut = 0.05) {
  
  pt_tbl <- get_pseudotimeDE_fdr_table(
    res = res
  )
  
  pt_tbl %>%
    dplyr::filter(
      is.finite(q_value),
      q_value < q_cut
    ) %>%
    dplyr::pull(
      Gene_ID
    ) %>%
    unique()
}

## ============================================================
## Extract Reactome pathway set
## ============================================================

get_reactome_pathway_set <- function(trajectory_result,
                                     method = c(
                                       "scFPC-DE",
                                       "PseudotimeDE"
                                     ),
                                     padj_cut = 0.05) {
  
  method <- match.arg(
    method
  )
  
  pathway_tbl <- if (method == "scFPC-DE") {
    trajectory_result$scFPCDE$table
  } else {
    trajectory_result$PseudotimeDE$table
  }
  
  if (
    is.null(pathway_tbl) ||
    nrow(pathway_tbl) == 0
  ) {
    return(
      character(0)
    )
  }
  
  pathway_tbl <- as.data.frame(
    pathway_tbl
  )
  
  if (!"ID" %in% names(pathway_tbl)) {
    stop(method, " Reactome table must contain ID")
  }
  
  if (!"p.adjust" %in% names(pathway_tbl)) {
    stop(method, " Reactome table must contain p.adjust")
  }
  
  pathway_tbl %>%
    dplyr::mutate(
      ID = as.character(ID),
      p.adjust = suppressWarnings(
        as.numeric(p.adjust)
      )
    ) %>%
    dplyr::filter(
      !is.na(ID),
      nzchar(ID),
      is.finite(p.adjust),
      p.adjust < padj_cut
    ) %>%
    dplyr::pull(
      ID
    ) %>%
    unique()
}

## ============================================================
## Build overlap object
## ============================================================

make_overlap_object <- function(sc_set,
                                pt_set) {
  
  sc_set <- unique(
    as.character(sc_set)
  )
  
  pt_set <- unique(
    as.character(pt_set)
  )
  
  sc_set <- sc_set[
    !is.na(sc_set) &
      nzchar(sc_set)
  ]
  
  pt_set <- pt_set[
    !is.na(pt_set) &
      nzchar(pt_set)
  ]
  
  common_set <- intersect(
    sc_set,
    pt_set
  )
  
  sc_only_set <- setdiff(
    sc_set,
    pt_set
  )
  
  pt_only_set <- setdiff(
    pt_set,
    sc_set
  )
  
  counts <- c(
    "scFPC-DE" = length(
      sc_only_set
    ),
    "PseudotimeDE" = length(
      pt_only_set
    ),
    "scFPC-DE&PseudotimeDE" = length(
      common_set
    )
  )
  
  list(
    sc_set = sc_set,
    pt_set = pt_set,
    common_set = common_set,
    sc_only_set = sc_only_set,
    pt_only_set = pt_only_set,
    counts = counts
  )
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
## Build TDEG overlap objects
## ============================================================

deg_overlap_traj1 <- make_overlap_object(
  sc_set = get_scfpcde_deg_set(
    res = traj_results$traj1,
    q_cut = fdr_cut
  ),
  pt_set = get_pseudotimeDE_deg_set(
    res = traj_results$traj1,
    q_cut = fdr_cut
  )
)

deg_overlap_traj2 <- make_overlap_object(
  sc_set = get_scfpcde_deg_set(
    res = traj_results$traj2,
    q_cut = fdr_cut
  ),
  pt_set = get_pseudotimeDE_deg_set(
    res = traj_results$traj2,
    q_cut = fdr_cut
  )
)

## ============================================================
## Build Reactome pathway overlap objects
## ============================================================

reactome_overlap_traj1 <- make_overlap_object(
  sc_set = get_reactome_pathway_set(
    trajectory_result = reactome_ora_results$traj1,
    method = "scFPC-DE",
    padj_cut = ora_padj_cut
  ),
  pt_set = get_reactome_pathway_set(
    trajectory_result = reactome_ora_results$traj1,
    method = "PseudotimeDE",
    padj_cut = ora_padj_cut
  )
)

reactome_overlap_traj2 <- make_overlap_object(
  sc_set = get_reactome_pathway_set(
    trajectory_result = reactome_ora_results$traj2,
    method = "scFPC-DE",
    padj_cut = ora_padj_cut
  ),
  pt_set = get_reactome_pathway_set(
    trajectory_result = reactome_ora_results$traj2,
    method = "PseudotimeDE",
    padj_cut = ora_padj_cut
  )
)

## ============================================================
## Summary tables
## ============================================================

deg_overlap_summary <- data.frame(
  Trajectory = c(
    "Trajectory 1",
    "Trajectory 2"
  ),
  scFPC_DE_only = c(
    length(deg_overlap_traj1$sc_only_set),
    length(deg_overlap_traj2$sc_only_set)
  ),
  Common = c(
    length(deg_overlap_traj1$common_set),
    length(deg_overlap_traj2$common_set)
  ),
  PseudotimeDE_only = c(
    length(deg_overlap_traj1$pt_only_set),
    length(deg_overlap_traj2$pt_only_set)
  ),
  scFPC_DE_total = c(
    length(deg_overlap_traj1$sc_set),
    length(deg_overlap_traj2$sc_set)
  ),
  PseudotimeDE_total = c(
    length(deg_overlap_traj1$pt_set),
    length(deg_overlap_traj2$pt_set)
  ),
  stringsAsFactors = FALSE
)

reactome_pathway_overlap_summary <- data.frame(
  Trajectory = c(
    "Trajectory 1",
    "Trajectory 2"
  ),
  scFPC_DE_only = c(
    length(reactome_overlap_traj1$sc_only_set),
    length(reactome_overlap_traj2$sc_only_set)
  ),
  Common = c(
    length(reactome_overlap_traj1$common_set),
    length(reactome_overlap_traj2$common_set)
  ),
  PseudotimeDE_only = c(
    length(reactome_overlap_traj1$pt_only_set),
    length(reactome_overlap_traj2$pt_only_set)
  ),
  scFPC_DE_total = c(
    length(reactome_overlap_traj1$sc_set),
    length(reactome_overlap_traj2$sc_set)
  ),
  PseudotimeDE_total = c(
    length(reactome_overlap_traj1$pt_set),
    length(reactome_overlap_traj2$pt_set)
  ),
  stringsAsFactors = FALSE
)

## ============================================================
## Consistency checks
## ============================================================

bh_tdeg_count_summary <- data.frame(
  Trajectory = c(
    "Trajectory 1",
    "Trajectory 2"
  ),
  scFPC_DE_tested = c(
    sum(is.finite(scfpcde_table_traj1$q_value)),
    sum(is.finite(scfpcde_table_traj2$q_value))
  ),
  scFPC_DE_significant = c(
    sum(
      is.finite(scfpcde_table_traj1$q_value) &
        scfpcde_table_traj1$q_value < fdr_cut
    ),
    sum(
      is.finite(scfpcde_table_traj2$q_value) &
        scfpcde_table_traj2$q_value < fdr_cut
    )
  ),
  PseudotimeDE_tested = c(
    sum(is.finite(pseudotimeDE_table_traj1$q_value)),
    sum(is.finite(pseudotimeDE_table_traj2$q_value))
  ),
  PseudotimeDE_significant = c(
    sum(
      is.finite(pseudotimeDE_table_traj1$q_value) &
        pseudotimeDE_table_traj1$q_value < fdr_cut
    ),
    sum(
      is.finite(pseudotimeDE_table_traj2$q_value) &
        pseudotimeDE_table_traj2$q_value < fdr_cut
    )
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  deg_overlap_summary$scFPC_DE_total ==
    bh_tdeg_count_summary$scFPC_DE_significant
)

stopifnot(
  deg_overlap_summary$PseudotimeDE_total ==
    bh_tdeg_count_summary$PseudotimeDE_significant
)

## ============================================================
## Fit Euler diagrams
## ============================================================

fit_deg_traj1 <- eulerr::euler(
  deg_overlap_traj1$counts,
  shape = "ellipse"
)

fit_deg_traj2 <- eulerr::euler(
  deg_overlap_traj2$counts,
  shape = "ellipse"
)

fit_reactome_traj1 <- eulerr::euler(
  reactome_overlap_traj1$counts,
  shape = "ellipse"
)

fit_reactome_traj2 <- eulerr::euler(
  reactome_overlap_traj2$counts,
  shape = "ellipse"
)

## ============================================================
## Euler plot object helper
## ============================================================

make_euler_plot_object <- function(fit) {
  
  plot(
    fit,
    fills = list(
      fill = euler_fill_cols,
      alpha = euler_alpha
    ),
    edges = list(
      col = "black",
      lwd = euler_edge_lwd
    ),
    labels = list(
      col = "black",
      cex = euler_label_cex
    ),
    quantities = list(
      col = "black",
      cex = euler_count_cex
    ),
    main = FALSE
  )
}

## ============================================================
## Create eulerr plot objects
## ============================================================

euler_degs_traj1 <- make_euler_plot_object(
  fit = fit_deg_traj1
)

euler_degs_traj2 <- make_euler_plot_object(
  fit = fit_deg_traj2
)

euler_reactome_traj1 <- make_euler_plot_object(
  fit = fit_reactome_traj1
)

euler_reactome_traj2 <- make_euler_plot_object(
  fit = fit_reactome_traj2
)

## ============================================================
## Arrange grobs with titles
## ============================================================

g_degs_traj1 <- gridExtra::arrangeGrob(
  grobs = list(
    euler_degs_traj1
  ),
  top = grid::textGrob(
    paste0(
      "Trajectory 1: Trans / Naive / C-mem1\n",
      "TDEG overlap between scFPC-DE and PseudotimeDE"
    ),
    gp = grid::gpar(
      fontsize = 15,
      fontface = "bold"
    )
  )
)

g_degs_traj2 <- gridExtra::arrangeGrob(
  grobs = list(
    euler_degs_traj2
  ),
  top = grid::textGrob(
    paste0(
      "Trajectory 2: M-mem1 / DN3 / DN2\n",
      "TDEG overlap between scFPC-DE and PseudotimeDE"
    ),
    gp = grid::gpar(
      fontsize = 15,
      fontface = "bold"
    )
  )
)

g_reactome_traj1 <- gridExtra::arrangeGrob(
  grobs = list(
    euler_reactome_traj1
  ),
  top = grid::textGrob(
    paste0(
      "Trajectory 1: Trans / Naive / C-mem1\n",
      "Significant Reactome pathway overlap"
    ),
    gp = grid::gpar(
      fontsize = 15,
      fontface = "bold"
    )
  )
)

g_reactome_traj2 <- gridExtra::arrangeGrob(
  grobs = list(
    euler_reactome_traj2
  ),
  top = grid::textGrob(
    paste0(
      "Trajectory 2: M-mem1 / DN3 / DN2\n",
      "Significant Reactome pathway overlap"
    ),
    gp = grid::gpar(
      fontsize = 15,
      fontface = "bold"
    )
  )
)

g_combined <- gridExtra::arrangeGrob(
  grobs = list(
    g_degs_traj1,
    g_degs_traj2,
    g_reactome_traj1,
    g_reactome_traj2
  ),
  ncol = 2,
  top = grid::textGrob(
    "Set-size overlap of TDEGs and Reactome pathways",
    gp = grid::gpar(
      fontsize = 19,
      fontface = "bold"
    )
  )
)

## ============================================================
## Display plots
## ============================================================

gridExtra::grid.arrange(
  g_degs_traj1
)

gridExtra::grid.arrange(
  g_degs_traj2
)

gridExtra::grid.arrange(
  g_reactome_traj1
)

gridExtra::grid.arrange(
  g_reactome_traj2
)

gridExtra::grid.arrange(
  g_combined
)

## ============================================================
## Save grob helper
## ============================================================

save_grob_tiff <- function(grob_obj,
                           filename,
                           width,
                           height,
                           dpi = 300) {
  
  filename <- file.path(
    getwd(),
    basename(filename)
  )
  
  if (file.exists(filename)) {
    unlink(filename)
  }
  
  grDevices::tiff(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  
  grid::grid.newpage()
  
  grid::grid.draw(
    grob_obj
  )
  
  grDevices::dev.off()
  
  if (!file.exists(filename)) {
    stop("File was not created: ", filename)
  }
  
  if (
    is.na(file.info(filename)$size) ||
    file.info(filename)$size == 0
  ) {
    stop("File is empty: ", filename)
  }
  
  invisible(
    filename
  )
}

## ============================================================
## Save TIFF files directly in current working directory
## ============================================================

file_traj1_deg <- save_grob_tiff(
  grob_obj = g_degs_traj1,
  filename = "traj1_DEG_overlap_euler.tiff",
  width = single_width,
  height = single_height,
  dpi = tiff_dpi
)

file_traj2_deg <- save_grob_tiff(
  grob_obj = g_degs_traj2,
  filename = "traj2_DEG_overlap_euler.tiff",
  width = single_width,
  height = single_height,
  dpi = tiff_dpi
)

file_traj1_reactome <- save_grob_tiff(
  grob_obj = g_reactome_traj1,
  filename = "traj1_Reactome_pathway_overlap_euler.tiff",
  width = single_width,
  height = single_height,
  dpi = tiff_dpi
)

file_traj2_reactome <- save_grob_tiff(
  grob_obj = g_reactome_traj2,
  filename = "traj2_Reactome_pathway_overlap_euler.tiff",
  width = single_width,
  height = single_height,
  dpi = tiff_dpi
)

file_combined <- save_grob_tiff(
  grob_obj = g_combined,
  filename = "DEG_and_Reactome_overlap_euler_2x2.tiff",
  width = combined_width,
  height = combined_height,
  dpi = tiff_dpi
)

## ============================================================
## Save CSV summaries directly in current working directory
## ============================================================

file_deg_summary <- file.path(
  getwd(),
  "DEG_overlap_summary.csv"
)

file_reactome_summary <- file.path(
  getwd(),
  "Reactome_pathway_overlap_summary.csv"
)

write.csv(
  deg_overlap_summary,
  file = file_deg_summary,
  row.names = FALSE
)

write.csv(
  reactome_pathway_overlap_summary,
  file = file_reactome_summary,
  row.names = FALSE
)

## ============================================================
## Confirm saved files
## ============================================================

saved_files <- c(
  file_traj1_deg,
  file_traj2_deg,
  file_traj1_reactome,
  file_traj2_reactome,
  file_combined,
  file_deg_summary,
  file_reactome_summary
)

cat("\nCurrent working directory:\n")
print(getwd())

cat("\nBH-adjusted TDEG counts:\n")
print(bh_tdeg_count_summary)

cat("\nTDEG overlap summary:\n")
print(deg_overlap_summary)

cat("\nReactome pathway overlap summary:\n")
print(reactome_pathway_overlap_summary)

cat("\nSaved files:\n")
print(saved_files)

cat("\nFile exists check:\n")
print(file.exists(saved_files))

cat("\nFile sizes:\n")
print(file.info(saved_files)$size)
