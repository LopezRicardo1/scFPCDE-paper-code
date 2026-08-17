# Simulation study and zero-inflation sensitivity analysis
#
# Extracted from the R code chunks in scFPCDE_Supplementary.Rmd.
# Analysis statements are preserved; only input paths are repository-relative.

library(knitr)
library(SingleCellExperiment)
library(scFPCDE)
library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(cowplot)
library(patchwork)

results <- readRDS("data/scFPCDE_simulation_results_final.rds")
options(digits = 4)

load("data/pseudotimeFPC.sim.data.rda")
# helper that returns requested summaries for one SCE object
summarise_sce <- function(sce, label) {
  counts_mat <- assay(sce, "counts")
  log_mat    <- log2(counts_mat + 1)

  tibble(
    Dataset         = label,
    Zero_Proportion = mean(counts_mat == 0),
    Min   = min(log_mat),
    Q1    = quantile(log_mat, 0.25),
    Median= median(log_mat),
    Q3    = quantile(log_mat, 0.75),
    Max   = max(log_mat),
    Mean  = mean(log_mat),
    SD    = sd(as.vector(log_mat))
  )
}

# names in sce.scd3.list are "y", "yz1", "yz2", "yz3"
labels <- c("ZI0", "ZI1", "ZI2", "ZI3")
sce_keys <- c("y",  "yz1",  "yz2",  "yz3")

zi_summary <- purrr::map2_dfr(
  sce_keys, labels,
  ~ summarise_sce(sce.scd3.list[[.x]]$scd3, .y)
)

knitr::kable(
  zi_summary,
  digits = c(NA, 3, rep(2, 7)),
  caption = "Summary of zero proportion and log transformed count distribution across the four simulated datasets (ZI0–ZI3)."
)

###############################################################################
#  scFPC-DE: gene smoothing & FPC curve visualization by D-statistic
###############################################################################
library(scFPCDE)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)
theme_set(theme_cowplot())
set.seed(1)
# Load counts & oracle pseudotime
sce    <- sce.scd3.list$y$scd3
Y_cnt  <- t(assay(sce, "counts"))
ptime  <- colData(sce)$pseudotime

# Log-transform & center
Y_log  <- log2(Y_cnt + 1)
Y_ctr  <- t(t(Y_log) - colMeans(Y_log))

# Filter & prepare
keep   <- scFPCDE_filter_genes(Y_cnt, qz = 0.10)
Y_use  <- Y_ctr[, keep, drop = FALSE]

# FPCA tuning and model fit
tun <- scFPCDE_tune_fpca(Y_use, ptime, L = 3, topvarper = .20)
fit <- scFPCDE_run(
  Y_use, ptime,
  L         = 3,
  nbasis    = tun$best_nbasis,
  r_pen     = tun$best_r_pen,
  n_perm    = 100,
  ncores = 2,
  topvarper = .20,
  center    = FALSE, scale = FALSE
)

# Rank by D-statistic
Dvals <- fit$D_test_result$D_obs
names(Dvals) <- fit$D_test_result$ID
top_genes <- names(sort(Dvals, decreasing = TRUE))[1:6]
bot_genes <- names(sort(Dvals))[1:6+2500]
topbot = c(top_genes, bot_genes)
Yt <- Y_log[, keep, drop = FALSE]
gene_means <- colMeans(Yt)
Xt_hat_centered <- fit$fpca_result$xt_hat
Xt_hat <- sweep(
  Xt_hat_centered,
  MARGIN = 2,
  STATS = gene_means[colnames(Xt_hat_centered)],
  FUN = "+"
)

library(ggplot2)
library(dplyr)
library(tidyr)

# Define pseudotime
tt <- seq(0, 1, length.out = 500)

# Define basis functions (signal generators)
phi1 <- exp(-tt/2) * cos(2 * pi * tt)
phi2 <- exp(-tt) * sin(2 * pi * tt)
phi3 <- exp(-2 * tt)

# Combine into a tidy data frame
basis_df <- data.frame(
  tt = tt,
  phi1 = phi1,
  phi2 = phi2,
  phi3 = phi3
) %>%
  pivot_longer(cols = starts_with("phi"),
               names_to = "Basis",
               values_to = "Value")

# Plot signal-generating basis functions
signals = ggplot(basis_df, aes(x = tt, y = Value, color = Basis)) +
  geom_line(size = 1.2) +
  theme_classic() +
  labs(
    x = "Pseudotime (t)",
    y = expression(phi[k](t)),
    title = "(a) DEGs Expression Basis Functions for Simulation Model"
  )+ theme_cowplot()



rankgenes = scFPCDE_gene_curves(
  tt            = ptime,
  yt            = Yt,
  yt_fit        = Xt_hat,
  cell_cluster  = rep(1, length(ptime)),
  facet_genes   = TRUE,
  subset        = topbot,
  ncol          = 4,
  nrow          = 3
) +
  theme_cowplot() +
  theme(legend.position = "none") +
  ggtitle("(b) Top Ranked by FPC-Distance DEGs vs Non DEGs", "Simulation")



signals/rankgenes



# Extract AUC values and reshape to long format
auc_df <- lapply(names(results), function(zilvl) {
  tibble(
    ZI     = zilvl,
    Method = names(results[[zilvl]]$auc),
    AUC    = as.numeric(results[[zilvl]]$auc)
  )
}) |>
  bind_rows() |>
  mutate(ZI = factor(ZI, levels = c("ZI0", "ZI1", "ZI2", "ZI3")),
         Method = factor(Method, levels = c("scFPC_DE", "FPC_F",
                                            "PseudotimeDE", "PseudotimeDE_ZINB", "PseudotimeDE_G")))
# 
# # Plot
# auc_plot <- ggplot(auc_df, aes(x = ZI, y = AUC, group = Method, color = Method, linetype = Method)) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   labs(title = "ROC-AUC vs. Zero-Inflation Level",
#        x = "Zero-Inflation Level", y = "Area Under Curve (AUC)") +
#   theme_cowplot() +
#   theme(legend.position = "right", aspect.ratio = 0.75)
# 
# print(auc_plot)

auc_plot <- ggplot(auc_df, aes(x = ZI, y = AUC, fill = Method)) +
  geom_col(position = "dodge", width = 0.7, alpha = .5, col = "black") +
  labs(
    title = "ROC-AUC vs. Zero-Inflation Level",
    x = "Zero-Inflation Level",
    y = "Area Under Curve (AUC)"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "right",
    aspect.ratio = 0.75
  )
print(auc_plot)

# Extract type I plots and suppress legends for the first three
typeI_list <- list(
  results$ZI0$typeI_plot + theme(legend.position = "none") + ggtitle("ZI0"),
  results$ZI1$typeI_plot + theme(legend.position = "none") + ggtitle("ZI1"),
  results$ZI2$typeI_plot + theme(legend.position = "bottom") + ggtitle("ZI2"),
  results$ZI3$typeI_plot + theme(legend.position = "none") + ggtitle("ZI3")
)

# Arrange into 2x2 grid
typeI_grid <- plot_grid(plotlist = typeI_list, ncol = 2, align = "hv")
# Display
print(typeI_grid)

# Extract power plots and suppress legends for the first three
power_list <- list(
  results$ZI0$power_plot + theme(legend.position = "none") + ggtitle("ZI0"),
  results$ZI1$power_plot + theme(legend.position = "none") + ggtitle("ZI1"),
  results$ZI2$power_plot + theme(legend.position = "bottom") + ggtitle("ZI2"),
  results$ZI3$power_plot + theme(legend.position = "none") + ggtitle("ZI3")
)

# Arrange into 2x2 grid
power_grid <- plot_grid(plotlist = power_list, ncol = 2, align = "hv")

# Display
print(power_grid)


library(patchwork)

p0 <- results$ZI0$score_plot + theme(legend.position = "none") +
  xlim(-3.5, 3.5) + ylim(-3.5, 3.5)

p1 <- results$ZI1$score_plot + theme(legend.position = "none") +
  xlim(-3.5, 3.5) + ylim(-3.5, 3.5)

p2 <- results$ZI2$score_plot + theme(legend.position = "none") +
  xlim(-3.5, 3.5) + ylim(-3.5, 3.5)

p3 <- results$ZI3$score_plot +
  xlim(-3.5, 3.5) + ylim(-3.5, 3.5)  # keep legend here

wrap_plots(p0, p1, p2, p3, ncol = 1)

## 3 ── run-time summary table ----------------------------------------------
timing_df <- bind_rows(
  results$ZI0$timing %>% mutate(ZI = "ZI0"),
  results$ZI1$timing %>% mutate(ZI = "ZI1"),
  results$ZI2$timing %>% mutate(ZI = "ZI2"),
  results$ZI3$timing %>% mutate(ZI = "ZI3")
) %>% 
  dplyr::select(ZI, Method, elapsed)

## Reorder methods for desired bar order
timing_df$Method <- factor(timing_df$Method,
  levels = c("scFPCDE", "PseudotimeDE_G", "PseudotimeDE", "PseudotimeDE_ZINB")
)

## 4 ── bar-plot of elapsed times -------------------------------------------
runtime_plot <- ggplot(timing_df, aes(x = Method, y = log(elapsed), fill = ZI)) +
  geom_col(position = position_dodge(), colour = "black", alpha = .5) +
  labs(title = "Elapsed time per method across Zero Inflation levels",
       x = NULL, y = "Log-Seconds") +
  theme_cowplot()

runtime_plot

timing_df |>
  group_by(Method) |>
  summarise(
    n     = n(),
    mean  = mean(elapsed, na.rm = TRUE),
    median = median(elapsed, na.rm = TRUE),
    sd    = sd(elapsed, na.rm = TRUE),
    min   = min(elapsed, na.rm = TRUE),
    max   = max(elapsed, na.rm = TRUE),
    iqr   = IQR(elapsed, na.rm = TRUE)
  ) |>
  kable(caption = "Summary stats in seconds across all zero inflation levels by method", digits = 2)
