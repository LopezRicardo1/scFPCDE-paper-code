## ============================================================
## Base R GeneRatio-versus-BgRatio Reactome plots
##
## Saves directly in current working directory:
##   traj1_scFPCDE_GeneRatio_BgRatio_baseR.tiff
##   traj1_PseudotimeDE_GeneRatio_BgRatio_baseR.tiff
##   traj2_scFPCDE_GeneRatio_BgRatio_baseR.tiff
##   traj2_PseudotimeDE_GeneRatio_BgRatio_baseR.tiff
##   Reactome_GeneRatio_BgRatio_2x2_baseR.tiff
##
## Requires:
##   reactome_ora_results$traj1
##   reactome_ora_results$traj2
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

graphics.off()

## ============================================================
## Settings
## ============================================================

ora_padj_cut <- 0.05

output_dir <- getwd()

n_unique_labels <- 3
n_common_labels_traj1 <- 0
n_common_labels_traj2 <- 3

pathway_status_cols <- c(
  "Common" = "blue3",
  "Unique" = "red2"
)

point_cex_min <- 1.50
point_cex_max <- 4

single_width <- 12
single_height <- 9

combined_width <- 10
combined_height <- 8

tiff_dpi <- 300

title_cex <- 1.50
axis_label_cex <- 1.35
axis_tick_cex <- 1.10
legend_text_cex <- 1.10

## Pathway label sizes
label_cex_unique_main <- 1.4
label_cex_common_main <- 1.4

label_cex_unique_combined <- 1.25
label_cex_common_combined <- 1.25

cat("\nCurrent working directory:\n")
print(getwd())

if (file.access(getwd(), 2) != 0) {
  stop("Current working directory is not writable: ", getwd())
}

## ============================================================
## Convert Reactome ratio strings to numeric values
## ============================================================

ratio_to_numeric <- function(x) {
  
  x <- as.character(x)
  
  vapply(
    strsplit(
      x,
      split = "/",
      fixed = TRUE
    ),
    function(z) {
      
      if (length(z) != 2) {
        return(NA_real_)
      }
      
      numerator <- suppressWarnings(
        as.numeric(z[1])
      )
      
      denominator <- suppressWarnings(
        as.numeric(z[2])
      )
      
      if (
        !is.finite(numerator) ||
        !is.finite(denominator) ||
        denominator <= 0
      ) {
        return(NA_real_)
      }
      
      numerator / denominator
    },
    FUN.VALUE = numeric(1)
  )
}

## ============================================================
## Extract significant Reactome pathways
## ============================================================

get_significant_reactome_table <- function(trajectory_result,
                                           method = c(
                                             "scFPC-DE",
                                             "PseudotimeDE"
                                           ),
                                           padj_cut = 0.05) {
  
  method <- match.arg(method)
  
  pathway_tbl <- if (method == "scFPC-DE") {
    trajectory_result$scFPCDE$table
  } else {
    trajectory_result$PseudotimeDE$table
  }
  
  if (
    is.null(pathway_tbl) ||
    nrow(pathway_tbl) == 0
  ) {
    return(data.frame())
  }
  
  pathway_tbl <- as.data.frame(pathway_tbl)
  
  if (!"GeneRatio_numeric" %in% names(pathway_tbl)) {
    pathway_tbl$GeneRatio_numeric <- ratio_to_numeric(
      pathway_tbl$GeneRatio
    )
  }
  
  if (!"BgRatio_numeric" %in% names(pathway_tbl)) {
    pathway_tbl$BgRatio_numeric <- ratio_to_numeric(
      pathway_tbl$BgRatio
    )
  }
  
  if ("qvalue" %in% names(pathway_tbl)) {
    
    pathway_tbl$q_value_plot <- suppressWarnings(
      as.numeric(pathway_tbl$qvalue)
    )
    
    replace_q <- !is.finite(
      pathway_tbl$q_value_plot
    )
    
    pathway_tbl$q_value_plot[replace_q] <- suppressWarnings(
      as.numeric(
        pathway_tbl$p.adjust[replace_q]
      )
    )
    
  } else {
    
    pathway_tbl$q_value_plot <- suppressWarnings(
      as.numeric(pathway_tbl$p.adjust)
    )
  }
  
  pathway_tbl %>%
    mutate(
      ID = as.character(ID),
      Description = as.character(Description),
      p.adjust = suppressWarnings(
        as.numeric(p.adjust)
      ),
      Count = suppressWarnings(
        as.numeric(Count)
      )
    ) %>%
    filter(
      !is.na(ID),
      nzchar(ID),
      !is.na(Description),
      nzchar(Description),
      is.finite(p.adjust),
      p.adjust < padj_cut,
      is.finite(q_value_plot),
      is.finite(GeneRatio_numeric),
      is.finite(BgRatio_numeric),
      GeneRatio_numeric >= 0,
      BgRatio_numeric >= 0
    ) %>%
    arrange(
      p.adjust,
      desc(Count)
    ) %>%
    distinct(
      ID,
      .keep_all = TRUE
    )
}

## ============================================================
## Prepare plotting data
## ============================================================

prepare_ratio_plot_data <- function(trajectory_result,
                                    method = c(
                                      "scFPC-DE",
                                      "PseudotimeDE"
                                    ),
                                    padj_cut = 0.05,
                                    n_unique_labels = 5,
                                    n_common_labels = 0) {
  
  method <- match.arg(method)
  
  sc_tbl <- get_significant_reactome_table(
    trajectory_result = trajectory_result,
    method = "scFPC-DE",
    padj_cut = padj_cut
  )
  
  pt_tbl <- get_significant_reactome_table(
    trajectory_result = trajectory_result,
    method = "PseudotimeDE",
    padj_cut = padj_cut
  )
  
  shared_ids <- intersect(
    sc_tbl$ID,
    pt_tbl$ID
  )
  
  plot_tbl <- if (method == "scFPC-DE") {
    sc_tbl
  } else {
    pt_tbl
  }
  
  if (nrow(plot_tbl) == 0) {
    
    return(
      list(
        data = data.frame(),
        unique_labels = data.frame(),
        common_labels = data.frame(),
        shared_ids = shared_ids
      )
    )
  }
  
  plot_tbl <- plot_tbl %>%
    mutate(
      Pathway_status = ifelse(
        ID %in% shared_ids,
        "Common",
        "Unique"
      ),
      minus_log10_q = -log10(
        pmax(
          q_value_plot,
          .Machine$double.xmin
        )
      )
    )
  
  unique_label_tbl <- plot_tbl %>%
    filter(
      Pathway_status == "Unique"
    ) %>%
    arrange(
      desc(BgRatio_numeric),
      q_value_plot,
      desc(Count)
    ) %>%
    slice_head(
      n = n_unique_labels
    )
  
  common_label_tbl <- plot_tbl %>%
    filter(
      Pathway_status == "Common"
    ) %>%
    arrange(
      desc(BgRatio_numeric),
      q_value_plot,
      desc(Count)
    ) %>%
    slice_head(
      n = n_common_labels
    )
  
  list(
    data = plot_tbl,
    unique_labels = unique_label_tbl,
    common_labels = common_label_tbl,
    shared_ids = shared_ids
  )
}

## ============================================================
## Common axis limit within trajectory
## ============================================================

get_ratio_axis_limit <- function(trajectory_result,
                                 padj_cut = 0.05,
                                 expansion = 1.18) {
  
  sc_tbl <- get_significant_reactome_table(
    trajectory_result = trajectory_result,
    method = "scFPC-DE",
    padj_cut = padj_cut
  )
  
  pt_tbl <- get_significant_reactome_table(
    trajectory_result = trajectory_result,
    method = "PseudotimeDE",
    padj_cut = padj_cut
  )
  
  ratio_values <- c(
    sc_tbl$GeneRatio_numeric,
    sc_tbl$BgRatio_numeric,
    pt_tbl$GeneRatio_numeric,
    pt_tbl$BgRatio_numeric
  )
  
  ratio_values <- ratio_values[
    is.finite(ratio_values) &
      ratio_values >= 0
  ]
  
  if (length(ratio_values) == 0) {
    return(0.10)
  }
  
  max(
    0.05,
    max(ratio_values) * expansion
  )
}

## ============================================================
## Point-size scaling
## ============================================================

scale_ratio_point_cex <- function(x,
                                  reference = x,
                                  cex_min = point_cex_min,
                                  cex_max = point_cex_max) {
  
  x <- as.numeric(x)
  reference <- as.numeric(reference)
  
  reference_range <- range(
    reference,
    finite = TRUE
  )
  
  if (
    length(x) == 0 ||
    any(!is.finite(reference_range))
  ) {
    
    return(
      rep(
        cex_min,
        length(x)
      )
    )
  }
  
  if (diff(reference_range) == 0) {
    
    return(
      rep(
        mean(
          c(
            cex_min,
            cex_max
          )
        ),
        length(x)
      )
    )
  }
  
  x <- pmin(
    pmax(
      x,
      reference_range[1]
    ),
    reference_range[2]
  )
  
  cex_min +
    (
      x - reference_range[1]
    ) /
    diff(reference_range) *
    (
      cex_max - cex_min
    )
}

## ============================================================
## Label position helper
## ============================================================

get_ratio_label_positions <- function(label_tbl,
                                      ratio_limit) {
  
  if (nrow(label_tbl) == 0) {
    return(integer(0))
  }
  
  x <- label_tbl$BgRatio_numeric
  y <- label_tbl$GeneRatio_numeric
  
  pos <- rep(
    c(2, 4, 3, 1),
    length.out = nrow(label_tbl)
  )
  
  pos[x >= 0.72 * ratio_limit] <- 2
  pos[x <= 0.15 * ratio_limit] <- 4
  pos[y >= 0.82 * ratio_limit] <- 1
  pos[y <= 0.10 * ratio_limit] <- 3
  
  pos
}

## ============================================================
## Base R plot function
## ============================================================

plot_gene_bg_ratio_base <- function(trajectory_result,
                                    method = c(
                                      "scFPC-DE",
                                      "PseudotimeDE"
                                    ),
                                    trajectory_title,
                                    panel_label = NULL,
                                    padj_cut = 0.05,
                                    n_unique_labels = 5,
                                    n_common_labels = 0,
                                    ratio_limit = NULL,
                                    point_alpha = 0.82,
                                    point_cex_min_use = point_cex_min,
                                    point_cex_max_use = point_cex_max,
                                    unique_label_cex = label_cex_unique_main,
                                    common_label_cex = label_cex_common_main,
                                    compact_panel = FALSE) {
  
  method <- match.arg(method)
  
  plot_obj <- prepare_ratio_plot_data(
    trajectory_result = trajectory_result,
    method = method,
    padj_cut = padj_cut,
    n_unique_labels = n_unique_labels,
    n_common_labels = n_common_labels
  )
  
  plot_tbl <- plot_obj$data
  unique_label_tbl <- plot_obj$unique_labels
  common_label_tbl <- plot_obj$common_labels
  
  if (compact_panel) {
    
    par(
      mar = c(4.4, 4.7, 4.0, 1.0),
      pty = "s",
      xaxs = "i",
      yaxs = "i",
      bg = "white",
      fg = "black",
      col.axis = "black",
      col.lab = "black",
      col.main = "black",
      xpd = FALSE,
      cex.axis = 1.05,
      cex.lab = 1.18,
      cex.main = 1.12
    )
    
  } else {
    
    par(
      mar = c(5.3, 5.4, 4.8, 1.4),
      pty = "s",
      xaxs = "i",
      yaxs = "i",
      bg = "white",
      fg = "black",
      col.axis = "black",
      col.lab = "black",
      col.main = "black",
      xpd = FALSE,
      cex.axis = axis_tick_cex,
      cex.lab = axis_label_cex,
      cex.main = title_cex
    )
  }
  
  if (nrow(plot_tbl) == 0) {
    
    plot.new()
    
    text(
      x = 0.5,
      y = 0.5,
      labels = paste0(
        trajectory_title,
        "\n",
        method,
        "\nNo significant Reactome pathways"
      ),
      cex = 1.1
    )
    
    return(
      invisible(plot_obj)
    )
  }
  
  if (is.null(ratio_limit)) {
    
    ratio_limit <- max(
      plot_tbl$GeneRatio_numeric,
      plot_tbl$BgRatio_numeric,
      na.rm = TRUE
    ) * 1.18
  }
  
  if (
    !is.finite(ratio_limit) ||
    ratio_limit <= 0
  ) {
    ratio_limit <- 0.10
  }
  
  point_cex <- scale_ratio_point_cex(
    x = plot_tbl$minus_log10_q,
    reference = plot_tbl$minus_log10_q,
    cex_min = point_cex_min_use,
    cex_max = point_cex_max_use
  )
  
  point_cols <- grDevices::adjustcolor(
    pathway_status_cols[
      plot_tbl$Pathway_status
    ],
    alpha.f = point_alpha
  )
  
  plot(
    x = plot_tbl$BgRatio_numeric,
    y = plot_tbl$GeneRatio_numeric,
    type = "n",
    xlim = c(-0.05, ratio_limit),
    ylim = c(0, 0.365),
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    ann = FALSE
  )
  abline(v = 0, lty = 2)
  segments(
    x0 = 0,
    y0 = 0,
    x1 = ratio_limit,
    y1 = ratio_limit,
    lty = 2,
    lwd = 1.4,
    col = "black",
    xpd = FALSE
  )
  
  points(
    x = plot_tbl$BgRatio_numeric,
    y = plot_tbl$GeneRatio_numeric,
    pch = 21,
    bg = point_cols,
    col = point_cols,
    cex = point_cex
  )
  
  axis(
    side = 1,
    las = 1
  )
  
  axis(
    side = 2,
    las = 1
  )
  
  box(
    col = "black",
    lwd = 1
  )
  
  title(
    main = paste0(
      trajectory_title,
      "\n",
      method,
      " - GeneRatio vs BgRatio"
    ),
    xlab = "BgRatio",
    ylab = "GeneRatio",
    font.main = 2
  )
  
  if (!is.null(panel_label)) {
    
    mtext(
      panel_label,
      side = 3,
      adj = 0,
      line = 1.0,
      font = 2,
      cex = ifelse(
        compact_panel,
        1.15,
        1.35
      )
    )
  }
  
  if (nrow(unique_label_tbl) > 0) {
    
    unique_label_positions <- get_ratio_label_positions(
      label_tbl = unique_label_tbl,
      ratio_limit = ratio_limit
    )
    
    text(
      x = unique_label_tbl$BgRatio_numeric,
      y = unique_label_tbl$GeneRatio_numeric,
      labels = unique_label_tbl$Description,
      pos = unique_label_positions,
      offset = 0.48,
      cex = unique_label_cex,
      font = 2,
      col = "black",
      xpd = FALSE
    )
  }
  
  if (nrow(common_label_tbl) > 0) {
    
    common_label_positions <- get_ratio_label_positions(
      label_tbl = common_label_tbl,
      ratio_limit = ratio_limit
    )
    
    text(
      x = common_label_tbl$BgRatio_numeric,
      y = common_label_tbl$GeneRatio_numeric,
      labels = common_label_tbl$Description,
      pos = common_label_positions,
      offset = 0.55,
      cex = common_label_cex,
      font = 2,
      col = "blue4",
      xpd = FALSE
    )
  }
  
  legend(
    "topleft",
    legend = c(
      "Common",
      "Unique"
    ),
    pch = 21,
    pt.bg = pathway_status_cols[
      c(
        "Common",
        "Unique"
      )
    ],
    col = pathway_status_cols[
      c(
        "Common",
        "Unique"
      )
    ],
    pt.cex = 1.65,
    bty = "n",
    cex = ifelse(
      compact_panel,
      0.98,
      legend_text_cex
    ),
    title = "Pathway status",
    title.col = "black"
  )
  
  q_range <- range(
    plot_tbl$minus_log10_q,
    finite = TRUE
  )
  
  q_breaks <- pretty(
    q_range,
    n = 4
  )
  
  q_breaks <- q_breaks[
    is.finite(q_breaks) &
      q_breaks > 0 &
      q_breaks >= q_range[1] &
      q_breaks <= q_range[2]
  ]
  
  if (length(q_breaks) > 0) {
    
    legend_cex <- scale_ratio_point_cex(
      x = q_breaks,
      reference = plot_tbl$minus_log10_q,
      cex_min = point_cex_min_use,
      cex_max = point_cex_max_use
    )
    
    legend(
      "bottomright",
      legend = format(
        q_breaks,
        digits = 2,
        trim = TRUE
      ),
      pch = 21,
      pt.bg = "gray65",
      col = "gray35",
      pt.cex = legend_cex,
      bty = "n",
      cex = ifelse(
        compact_panel,
        1.02,
        1.18
      ),
      title = expression(-log[10]("q-value")),
      title.col = "black"
    )
  }
  
  invisible(
    list(
      data = plot_tbl,
      unique_labels = unique_label_tbl,
      common_labels = common_label_tbl,
      ratio_limit = ratio_limit
    )
  )
}

## ============================================================
## Common axis limits
## ============================================================

ratio_limit_traj1 <- get_ratio_axis_limit(
  trajectory_result = reactome_ora_results$traj1,
  padj_cut = ora_padj_cut,
  expansion = 1.18
)

ratio_limit_traj2 <- get_ratio_axis_limit(
  trajectory_result = reactome_ora_results$traj2,
  padj_cut = ora_padj_cut,
  expansion = 1.18
)

## ============================================================
## Display individual plots
## ============================================================

ratio_scfpcde_traj1 <- plot_gene_bg_ratio_base(
  trajectory_result = reactome_ora_results$traj1,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  padj_cut = ora_padj_cut,
  n_unique_labels = 4,
  n_common_labels = 0,
  ratio_limit = ratio_limit_traj1,
  compact_panel = FALSE,
  unique_label_cex = label_cex_unique_main,
  common_label_cex = label_cex_common_main
)

ratio_ptime_traj1 <- plot_gene_bg_ratio_base(
  trajectory_result = reactome_ora_results$traj1,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 0,
  ratio_limit = ratio_limit_traj1,
  compact_panel = FALSE,
  unique_label_cex = label_cex_unique_main,
  common_label_cex = label_cex_common_main
)

ratio_scfpcde_traj2 <- plot_gene_bg_ratio_base(
  trajectory_result = reactome_ora_results$traj2,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 5,
  ratio_limit = ratio_limit_traj2,
  compact_panel = FALSE,
  unique_label_cex = label_cex_unique_main,
  common_label_cex = label_cex_common_main
)

ratio_ptime_traj2 <- plot_gene_bg_ratio_base(
  trajectory_result = reactome_ora_results$traj2,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 0,
  ratio_limit = ratio_limit_traj2,
  compact_panel = FALSE,
  unique_label_cex = label_cex_unique_main,
  common_label_cex = label_cex_common_main
)

## ============================================================
## Save helper
## ============================================================

save_ratio_base_tiff <- function(filename,
                                 trajectory_result,
                                 method,
                                 trajectory_title,
                                 panel_label = NULL,
                                 padj_cut = 0.05,
                                 n_unique_labels = 5,
                                 n_common_labels = 0,
                                 ratio_limit = NULL,
                                 width = 7.4,
                                 height = 7.2,
                                 dpi = 300,
                                 compact_panel = FALSE,
                                 unique_label_cex = label_cex_unique_main,
                                 common_label_cex = label_cex_common_main) {
  
  filename <- file.path(
    output_dir,
    basename(filename)
  )
  
  if (file.exists(filename)) {
    unlink(filename)
  }
  
  if (capabilities("cairo")) {
    
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
    
  } else {
    
    grDevices::tiff(
      filename = filename,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white"
    )
  }
  
  plot_gene_bg_ratio_base(
    trajectory_result = trajectory_result,
    method = method,
    trajectory_title = trajectory_title,
    panel_label = panel_label,
    padj_cut = padj_cut,
    n_unique_labels = n_unique_labels,
    n_common_labels = n_common_labels,
    ratio_limit = ratio_limit,
    compact_panel = compact_panel,
    unique_label_cex = unique_label_cex,
    common_label_cex = common_label_cex
  )
  
  grDevices::dev.off()
  
  if (!file.exists(filename)) {
    stop("File was not created: ", filename)
  }
  
  if (is.na(file.info(filename)$size) || file.info(filename)$size == 0) {
    stop("File is empty: ", filename)
  }
  
  invisible(
    filename
  )
}

## ============================================================
## Save individual TIFF files
## ============================================================

file_traj1_scfpcde <- save_ratio_base_tiff(
  filename = "traj1_scFPCDE_GeneRatio_BgRatio_baseR.tiff",
  trajectory_result = reactome_ora_results$traj1,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  panel_label = NULL,
  padj_cut = ora_padj_cut,
  n_unique_labels = 4,
  n_common_labels = 0,
  ratio_limit = ratio_limit_traj1,
  width = single_width,
  height = single_height,
  dpi = tiff_dpi,
  compact_panel = FALSE,
  unique_label_cex = label_cex_unique_main,
  common_label_cex = label_cex_common_main
)

file_traj1_ptime <- save_ratio_base_tiff(
  filename = "traj1_PseudotimeDE_GeneRatio_BgRatio_baseR.tiff",
  trajectory_result = reactome_ora_results$traj1,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  panel_label = NULL,
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 0,
  ratio_limit = ratio_limit_traj1,
  width = single_width,
  height = single_height,
  dpi = tiff_dpi,
  compact_panel = FALSE,
  unique_label_cex = label_cex_unique_main,
  common_label_cex = label_cex_common_main
)

file_traj2_scfpcde <- save_ratio_base_tiff(
  filename = "traj2_scFPCDE_GeneRatio_BgRatio_baseR.tiff",
  trajectory_result = reactome_ora_results$traj2,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  panel_label = NULL,
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 5,
  ratio_limit = ratio_limit_traj2,
  width = single_width,
  height = single_height,
  dpi = tiff_dpi,
  compact_panel = FALSE,
  unique_label_cex = label_cex_unique_main,
  common_label_cex = label_cex_common_main
)

file_traj2_ptime <- save_ratio_base_tiff(
  filename = "traj2_PseudotimeDE_GeneRatio_BgRatio_baseR.tiff",
  trajectory_result = reactome_ora_results$traj2,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  panel_label = NULL,
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 0,
  ratio_limit = ratio_limit_traj2,
  width = single_width,
  height = single_height,
  dpi = tiff_dpi,
  compact_panel = FALSE,
  unique_label_cex = label_cex_unique_main,
  common_label_cex = label_cex_common_main
)

## ============================================================
## Save combined 2 x 2 multipanel TIFF
## ============================================================

file_combined <- file.path(
  output_dir,
  "Reactome_GeneRatio_BgRatio_2x2_baseR.tiff"
)

if (file.exists(file_combined)) {
  unlink(file_combined)
}

if (capabilities("cairo")) {
  
  grDevices::tiff(
    filename = file_combined,
    width = combined_width,
    height = combined_height,
    units = "in",
    res = tiff_dpi,
    compression = "lzw",
    type = "cairo",
    bg = "white"
  )
  
} else {
  
  grDevices::tiff(
    filename = file_combined,
    width = combined_width,
    height = combined_height,
    units = "in",
    res = tiff_dpi,
    bg = "white"
  )
}

old_par <- par(
  no.readonly = TRUE
)

par(
  mfrow = c(2, 2),
  oma = c(0.4, 0.4, 2.8, 0.4)
)

plot_gene_bg_ratio_base(
  trajectory_result = reactome_ora_results$traj1,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  panel_label = "(a)",
  padj_cut = ora_padj_cut,
  n_unique_labels = 4,
  n_common_labels = 0,
  ratio_limit = ratio_limit_traj1,
  point_cex_min_use = 1.00,
  point_cex_max_use = 3.30,
  unique_label_cex = label_cex_unique_combined,
  common_label_cex = label_cex_common_combined,
  compact_panel = TRUE
)

plot_gene_bg_ratio_base(
  trajectory_result = reactome_ora_results$traj1,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 1: Trans / Naive / C-mem1",
  panel_label = "(b)",
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 0,
  ratio_limit = ratio_limit_traj1,
  point_cex_min_use = 1.00,
  point_cex_max_use = 3.30,
  unique_label_cex = label_cex_unique_combined,
  common_label_cex = label_cex_common_combined,
  compact_panel = TRUE
)

plot_gene_bg_ratio_base(
  trajectory_result = reactome_ora_results$traj2,
  method = "scFPC-DE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  panel_label = "(c)",
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 5,
  ratio_limit = ratio_limit_traj2,
  point_cex_min_use = 1.00,
  point_cex_max_use = 3.30,
  unique_label_cex = label_cex_unique_combined,
  common_label_cex = label_cex_common_combined,
  compact_panel = TRUE
)

plot_gene_bg_ratio_base(
  trajectory_result = reactome_ora_results$traj2,
  method = "PseudotimeDE",
  trajectory_title = "Trajectory 2: M-mem1 / DN3 / DN2",
  panel_label = "(d)",
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 0,
  ratio_limit = ratio_limit_traj2,
  point_cex_min_use = 1.00,
  point_cex_max_use = 3.30,
  unique_label_cex = label_cex_unique_combined,
  common_label_cex = label_cex_common_combined,
  compact_panel = TRUE
)

mtext(
  "Reactome GeneRatio versus BgRatio by method and trajectory",
  outer = TRUE,
  side = 3,
  line = 0.9,
  font = 2,
  cex = 1.45
)

par(
  old_par
)

grDevices::dev.off()

## ============================================================
## Save labeled data tables
## ============================================================

label_summary_traj1_scfpcde <- prepare_ratio_plot_data(
  trajectory_result = reactome_ora_results$traj1,
  method = "scFPC-DE",
  padj_cut = ora_padj_cut,
  n_unique_labels = 4,
  n_common_labels = 0
)

label_summary_traj1_ptime <- prepare_ratio_plot_data(
  trajectory_result = reactome_ora_results$traj1,
  method = "PseudotimeDE",
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 0
)

label_summary_traj2_scfpcde <- prepare_ratio_plot_data(
  trajectory_result = reactome_ora_results$traj2,
  method = "scFPC-DE",
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 5
)

label_summary_traj2_ptime <- prepare_ratio_plot_data(
  trajectory_result = reactome_ora_results$traj2,
  method = "PseudotimeDE",
  padj_cut = ora_padj_cut,
  n_unique_labels = 0,
  n_common_labels = 0
)

write.csv(
  label_summary_traj1_scfpcde$unique_labels,
  file = file.path(
    output_dir,
    "traj1_scFPCDE_unique_pathways_labeled.csv"
  ),
  row.names = FALSE
)

write.csv(
  label_summary_traj2_scfpcde$common_labels,
  file = file.path(
    output_dir,
    "traj2_scFPCDE_common_pathways_labeled.csv"
  ),
  row.names = FALSE
)

## ============================================================
## Confirm files
## ============================================================

saved_files <- c(
  file_traj1_scfpcde,
  file_traj1_ptime,
  file_traj2_scfpcde,
  file_traj2_ptime,
  file_combined,
  file.path(
    output_dir,
    "traj1_scFPCDE_unique_pathways_labeled.csv"
  ),
  file.path(
    output_dir,
    "traj2_scFPCDE_common_pathways_labeled.csv"
  )
)

cat("\nSaved files:\n")

print(
  normalizePath(
    saved_files,
    mustWork = FALSE
  )
)

cat("\nFile exists check:\n")

print(
  file.exists(
    saved_files
  )
)

cat("\nFile sizes:\n")

print(
  file.info(
    saved_files
  )$size
)
