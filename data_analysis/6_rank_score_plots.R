## ============================================================
## Method colors
## ============================================================

method_cols_rank <- c(
  "scFPC-DE only"      = "red3",
  "PseudotimeDE only" = "blue3",
  "Both"              = "purple4",
  "Background"        = "gray75"
)

## ============================================================
## FPC score-space plot
## ============================================================

plot_ranked_genes_fpc_space <- function(res,
                                        rank_obj,
                                        trajectory_title,
                                        method_cols_use,
                                        label_top = TRUE,
                                        label_cex = 0.75,
                                        point_cex_bg = 0.55,
                                        point_cex_top = 1.05) {
  
  score_df <- get_fpca_gene_scores(res)
  gene_annot <- make_rank_gene_annotation(rank_obj)
  
  plot_df <- score_df %>%
    left_join(
      gene_annot %>%
        dplyr::select(Gene_ID, Method_Class, scFPC_DE_Rank, PseudotimeDE_Rank),
      by = "Gene_ID"
    ) %>%
    mutate(
      Method_Class = ifelse(is.na(Method_Class), "Background", Method_Class),
      Method_Class = factor(
        Method_Class,
        levels = c("Background", "scFPC-DE only", "PseudotimeDE only", "Both")
      )
    )
  
  xlim <- range(plot_df$FPC1, finite = TRUE)
  ylim <- range(plot_df$FPC2, finite = TRUE)
  
  xpad <- 0.08 * diff(xlim)
  ypad <- 0.08 * diff(ylim)
  
  old_par <- par(c("mar", "bg", "fg", "col.axis", "col.lab", "col.main"))
  on.exit(par(old_par), add = TRUE)
  
  par(
    mar = c(5, 5, 4, 2),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black"
  )
  
  plot(
    plot_df$FPC1,
    plot_df$FPC2,
    type = "n",
    xlab = "Gene FPC1 score",
    ylab = "Gene FPC2 score",
    main = paste0(trajectory_title, "\nTop-ranked genes in FPC score space"),
    xlim = c(xlim[1] - xpad, xlim[2] + xpad),
    ylim = c(ylim[1] - ypad, ylim[2] + ypad)
  )
  
  abline(v = 0, h = 0, lty = 2, col = "gray60")
  
  bg_df <- plot_df %>%
    filter(Method_Class == "Background")
  
  points(
    bg_df$FPC1,
    bg_df$FPC2,
    pch = 16,
    col = method_cols_use["Background"],
    cex = point_cex_bg
  )
  
  for (mc in c("scFPC-DE only", "PseudotimeDE only", "Both")) {
    
    tmp <- plot_df %>%
      filter(Method_Class == mc)
    
    if (nrow(tmp) > 0) {
      points(
        tmp$FPC1,
        tmp$FPC2,
        pch = 16,
        col = method_cols_use[mc],
        cex = point_cex_top
      )
    }
  }
  
  if (label_top) {
    
    lab_df <- plot_df %>%
      filter(Method_Class != "Background")
    
    if (nrow(lab_df) > 0) {
      text(
        lab_df$FPC1,
        lab_df$FPC2,
        labels = lab_df$Gene_ID,
        col = method_cols_use[as.character(lab_df$Method_Class)],
        pos = 3,
        cex = label_cex,
        font = 2
      )
    }
  }
  
  legend(
    "topleft",
    legend = c("scFPC-DE only", "PseudotimeDE only", "Both", "Background"),
    col = method_cols_use[c("scFPC-DE only", "PseudotimeDE only", "Both", "Background")],
    pch = 16,
    pt.cex = c(1.2, 1.2, 1.2, 0.9),
    bty = "n"
  )
  
  invisible(plot_df)
}

## ============================================================
## Run plots
## ============================================================

traj1_fpc_rank_plot <- plot_ranked_genes_fpc_space(
  res = traj_results$traj1,
  rank_obj = rank_tables_minimal$traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  method_cols_use = method_cols_rank
)

traj2_fpc_rank_plot <- plot_ranked_genes_fpc_space(
  res = traj_results$traj2,
  rank_obj = rank_tables_minimal$traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  method_cols_use = method_cols_rank
)





## ============================================================
## Save FPC rank-space plots as TIFF
## ============================================================

save_ranked_genes_fpc_space_tiff <- function(filename,
                                             res,
                                             rank_obj,
                                             trajectory_title,
                                             method_cols_use,
                                             width = 8,
                                             height = 6,
                                             dpi = 300,
                                             label_top = TRUE,
                                             label_cex = 0.75,
                                             point_cex_bg = 0.55,
                                             point_cex_top = 1.05) {
  
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
  
  plot_ranked_genes_fpc_space(
    res = res,
    rank_obj = rank_obj,
    trajectory_title = trajectory_title,
    method_cols_use = method_cols_use,
    label_top = label_top,
    label_cex = label_cex,
    point_cex_bg = point_cex_bg,
    point_cex_top = point_cex_top
  )
}

## ============================================================
## Display plots in RStudio Plots pane
## ============================================================

traj1_fpc_rank_plot <- plot_ranked_genes_fpc_space(
  res = traj_results$traj1,
  rank_obj = rank_tables_minimal$traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  method_cols_use = method_cols_rank
)

traj2_fpc_rank_plot <- plot_ranked_genes_fpc_space(
  res = traj_results$traj2,
  rank_obj = rank_tables_minimal$traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  method_cols_use = method_cols_rank
)

## ============================================================
## Save plots as TIFF
## ============================================================

save_ranked_genes_fpc_space_tiff(
  filename = "traj1_ranked_genes_FPC_score_space.tiff",
  res = traj_results$traj1,
  rank_obj = rank_tables_minimal$traj1,
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  method_cols_use = method_cols_rank,
  width = 8,
  height = 6,
  dpi = 300
)

save_ranked_genes_fpc_space_tiff(
  filename = "traj2_ranked_genes_FPC_score_space.tiff",
  res = traj_results$traj2,
  rank_obj = rank_tables_minimal$traj2,
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  method_cols_use = method_cols_rank,
  width = 8,
  height = 6,
  dpi = 300
)

## ============================================================
## Confirm saved files
## ============================================================

saved_rank_fpc_files <- c(
  "traj1_ranked_genes_FPC_score_space.tiff",
  "traj2_ranked_genes_FPC_score_space.tiff"
)

print(
  normalizePath(
    saved_rank_fpc_files,
    mustWork = FALSE
  )
)

print(
  file.exists(saved_rank_fpc_files)
)
