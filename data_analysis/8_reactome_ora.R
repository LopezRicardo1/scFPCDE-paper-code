## ============================================================
## STEP 1: Reactome over-representation analysis
##
## Runs ReactomePA ORA for:
##
##   Trajectory 1:
##     - scFPC-DE
##     - PseudotimeDE
##
##   Trajectory 2:
##     - scFPC-DE
##     - PseudotimeDE
##
## Required objects:
##   traj_results$traj1
##   traj_results$traj2
##
## DEG definition:
##   Method-specific q-value < deg_q_cut
##
## ORA background:
##   Genes tested by both methods within each trajectory,
##   restricted to genes successfully mapped to Entrez IDs
##
## Main output:
##   reactome_ora_results
##
## No plots are generated in this step.
## ============================================================

suppressPackageStartupMessages({
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
})

## ============================================================
## User settings
## ============================================================

deg_q_cut <- 0.05

## Retain the complete ORA result and apply significance afterward
ora_pvalue_cut <- 1
ora_qvalue_cut <- 1

## Significant pathway threshold
ora_padj_cut <- 0.05

min_gs_size <- 10
max_gs_size <- 500

save_csv_files <- TRUE
output_dir <- "Reactome_ORA_results"

if (save_csv_files && !dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

## ============================================================
## Convert Reactome ratio strings to numeric values
##
## Example:
##   "12/150" -> 0.08
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
## Clean gene-symbol vector
## ============================================================

clean_gene_symbols <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  unique(
    x[
      !is.na(x) &
        nzchar(x)
    ]
  )
}

## ============================================================
## Map human gene symbols to Entrez IDs
## ============================================================

map_symbols_to_entrez <- function(symbols) {
  
  symbols <- clean_gene_symbols(
    symbols
  )
  
  if (length(symbols) == 0) {
    
    return(
      data.frame(
        SYMBOL = character(0),
        ENTREZID = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  mapped <- suppressMessages(
    AnnotationDbi::select(
      x = org.Hs.eg.db,
      keys = symbols,
      keytype = "SYMBOL",
      columns = c(
        "SYMBOL",
        "ENTREZID"
      )
    )
  )
  
  mapped %>%
    as.data.frame() %>%
    transmute(
      SYMBOL = toupper(
        trimws(
          as.character(SYMBOL)
        )
      ),
      ENTREZID = as.character(
        ENTREZID
      )
    ) %>%
    filter(
      !is.na(SYMBOL),
      nzchar(SYMBOL),
      !is.na(ENTREZID),
      nzchar(ENTREZID)
    ) %>%
    distinct(
      SYMBOL,
      ENTREZID
    )
}

## ============================================================
## Extract scFPC-DE gene-level results
## ============================================================

extract_scfpcde_table <- function(res) {
  
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
    "p_value",
    "q_value"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(sc_tbl)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      "Missing scFPC-DE columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  if (!"D_obs" %in% names(sc_tbl)) {
    sc_tbl$D_obs <- NA_real_
  }
  
  sc_tbl %>%
    transmute(
      Gene_ID = toupper(
        trimws(
          as.character(ID)
        )
      ),
      p_value = suppressWarnings(
        as.numeric(p_value)
      ),
      q_value = suppressWarnings(
        as.numeric(q_value)
      ),
      D_obs = suppressWarnings(
        as.numeric(D_obs)
      )
    ) %>%
    filter(
      !is.na(Gene_ID),
      nzchar(Gene_ID)
    ) %>%
    arrange(
      q_value,
      p_value,
      desc(D_obs)
    ) %>%
    distinct(
      Gene_ID,
      .keep_all = TRUE
    )
}

## ============================================================
## Extract PseudotimeDE gene-level results
## ============================================================

extract_pseudotimeDE_table <- function(res) {
  
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
      "Missing PseudotimeDE columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  ## Use previously stored q-values when aligned with res.gauss
  if (
    !is.null(res$gauss.qvals) &&
    length(res$gauss.qvals) == nrow(pt_tbl)
  ) {
    
    pt_qvalues <- suppressWarnings(
      as.numeric(res$gauss.qvals)
    )
    
  } else {
    
    pt_qvalues <- stats::p.adjust(
      suppressWarnings(
        as.numeric(pt_tbl$fix.pv)
      ),
      method = "BH"
    )
  }
  
  pt_tbl %>%
    transmute(
      Gene_ID = toupper(
        trimws(
          as.character(gene)
        )
      ),
      p_value = suppressWarnings(
        as.numeric(fix.pv)
      ),
      q_value = pt_qvalues
    ) %>%
    filter(
      !is.na(Gene_ID),
      nzchar(Gene_ID)
    ) %>%
    arrange(
      q_value,
      p_value
    ) %>%
    distinct(
      Gene_ID,
      .keep_all = TRUE
    )
}

## ============================================================
## Extract significant genes and common tested-gene universe
## ============================================================

extract_ora_gene_sets <- function(res,
                                  q_cut = 0.05) {
  
  sc_tbl <- extract_scfpcde_table(
    res = res
  )
  
  pt_tbl <- extract_pseudotimeDE_table(
    res = res
  )
  
  ## Common tested genes ensure both methods use the same universe
  common_universe_symbols <- intersect(
    sc_tbl$Gene_ID,
    pt_tbl$Gene_ID
  )
  
  if (length(common_universe_symbols) == 0) {
    
    stop(
      "No common tested genes were found between ",
      "scFPC-DE and PseudotimeDE"
    )
  }
  
  sc_deg_symbols <- sc_tbl %>%
    filter(
      Gene_ID %in% common_universe_symbols,
      is.finite(q_value),
      q_value < q_cut
    ) %>%
    pull(Gene_ID) %>%
    unique()
  
  pt_deg_symbols <- pt_tbl %>%
    filter(
      Gene_ID %in% common_universe_symbols,
      is.finite(q_value),
      q_value < q_cut
    ) %>%
    pull(Gene_ID) %>%
    unique()
  
  list(
    sc_table = sc_tbl,
    pt_table = pt_tbl,
    universe_symbols = common_universe_symbols,
    sc_deg_symbols = sc_deg_symbols,
    pt_deg_symbols = pt_deg_symbols
  )
}

## ============================================================
## Run Reactome ORA for one method
## ============================================================

run_reactome_ora <- function(deg_symbols,
                             universe_symbols,
                             trajectory_name,
                             method_name,
                             min_gs_size = 10,
                             max_gs_size = 500,
                             pvalue_cutoff = 1,
                             qvalue_cutoff = 1,
                             pathway_padj_cut = 0.05) {
  
  deg_symbols <- clean_gene_symbols(
    deg_symbols
  )
  
  universe_symbols <- clean_gene_symbols(
    universe_symbols
  )
  
  ## ----------------------------------------------------------
  ## Map the common tested-gene universe
  ## ----------------------------------------------------------
  
  universe_mapping <- map_symbols_to_entrez(
    universe_symbols
  )
  
  if (nrow(universe_mapping) == 0) {
    
    stop(
      trajectory_name,
      ": no background genes mapped to Entrez IDs"
    )
  }
  
  universe_entrez <- unique(
    universe_mapping$ENTREZID
  )
  
  ## ----------------------------------------------------------
  ## Map significant genes using the same universe mapping
  ## ----------------------------------------------------------
  
  deg_mapping <- universe_mapping %>%
    filter(
      SYMBOL %in% deg_symbols
    )
  
  deg_entrez <- unique(
    deg_mapping$ENTREZID
  )
  
  mapping_summary <- data.frame(
    Trajectory = trajectory_name,
    Method = method_name,
    Significant_symbols = length(deg_symbols),
    Mapped_significant_Entrez = length(deg_entrez),
    Tested_symbols = length(universe_symbols),
    Mapped_background_Entrez = length(universe_entrez),
    stringsAsFactors = FALSE
  )
  
  print(
    mapping_summary
  )
  
  if (length(deg_entrez) < 2) {
    
    warning(
      trajectory_name,
      " | ",
      method_name,
      ": fewer than two significant genes mapped to Entrez IDs"
    )
    
    return(
      list(
        enrichment = NULL,
        table = data.frame(),
        significant_table = data.frame(),
        deg_symbols = deg_symbols,
        deg_entrez = deg_entrez,
        universe_symbols = universe_symbols,
        universe_entrez = universe_entrez,
        deg_mapping = deg_mapping,
        universe_mapping = universe_mapping,
        mapping_summary = mapping_summary
      )
    )
  }
  
  ## ----------------------------------------------------------
  ## Reactome over-representation analysis
  ## ----------------------------------------------------------
  
  enrichment <- ReactomePA::enrichPathway(
    gene = as.character(deg_entrez),
    universe = as.character(universe_entrez),
    organism = "human",
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalue_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    readable = TRUE
  )
  
  enrichment_table <- as.data.frame(
    enrichment
  )
  
  if (nrow(enrichment_table) > 0) {
    
    enrichment_table <- enrichment_table %>%
      mutate(
        Trajectory = trajectory_name,
        Method = method_name,
        GeneRatio_numeric = ratio_to_numeric(
          GeneRatio
        ),
        BgRatio_numeric = ratio_to_numeric(
          BgRatio
        )
      ) %>%
      arrange(
        p.adjust,
        pvalue,
        desc(Count)
      )
  }
  
  significant_table <- enrichment_table %>%
    filter(
      !is.na(p.adjust),
      is.finite(p.adjust),
      p.adjust < pathway_padj_cut
    )
  
  list(
    enrichment = enrichment,
    table = enrichment_table,
    significant_table = significant_table,
    deg_symbols = deg_symbols,
    deg_entrez = deg_entrez,
    universe_symbols = universe_symbols,
    universe_entrez = universe_entrez,
    deg_mapping = deg_mapping,
    universe_mapping = universe_mapping,
    mapping_summary = mapping_summary
  )
}

## ============================================================
## Run both methods for one trajectory
## ============================================================

run_trajectory_reactome_ora <- function(res,
                                        trajectory_name,
                                        q_cut = 0.05,
                                        pathway_padj_cut = 0.05,
                                        min_gs_size = 10,
                                        max_gs_size = 500) {
  
  gene_sets <- extract_ora_gene_sets(
    res = res,
    q_cut = q_cut
  )
  
  sc_result <- run_reactome_ora(
    deg_symbols = gene_sets$sc_deg_symbols,
    universe_symbols = gene_sets$universe_symbols,
    trajectory_name = trajectory_name,
    method_name = "scFPC-DE",
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size,
    pvalue_cutoff = ora_pvalue_cut,
    qvalue_cutoff = ora_qvalue_cut,
    pathway_padj_cut = pathway_padj_cut
  )
  
  pt_result <- run_reactome_ora(
    deg_symbols = gene_sets$pt_deg_symbols,
    universe_symbols = gene_sets$universe_symbols,
    trajectory_name = trajectory_name,
    method_name = "PseudotimeDE",
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size,
    pvalue_cutoff = ora_pvalue_cut,
    qvalue_cutoff = ora_qvalue_cut,
    pathway_padj_cut = pathway_padj_cut
  )
  
  list(
    trajectory_name = trajectory_name,
    gene_sets = gene_sets,
    scFPCDE = sc_result,
    PseudotimeDE = pt_result
  )
}

## ============================================================
## Run Reactome ORA for both trajectories
## ============================================================

reactome_ora_results <- list(
  
  traj1 = run_trajectory_reactome_ora(
    res = traj_results$traj1,
    trajectory_name = "Trajectory 1: Trans / Naive / C-mem1",
    q_cut = deg_q_cut,
    pathway_padj_cut = ora_padj_cut,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size
  ),
  
  traj2 = run_trajectory_reactome_ora(
    res = traj_results$traj2,
    trajectory_name = "Trajectory 2: M-mem1 / DN3 / DN2",
    q_cut = deg_q_cut,
    pathway_padj_cut = ora_padj_cut,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size
  )
)

## ============================================================
## Reactome ORA summary
## ============================================================

make_ora_summary_row <- function(trajectory_result,
                                 trajectory_label,
                                 method = c(
                                   "scFPC-DE",
                                   "PseudotimeDE"
                                 )) {
  
  method <- match.arg(method)
  
  method_result <- if (method == "scFPC-DE") {
    trajectory_result$scFPCDE
  } else {
    trajectory_result$PseudotimeDE
  }
  
  significant_symbols <- if (method == "scFPC-DE") {
    trajectory_result$gene_sets$sc_deg_symbols
  } else {
    trajectory_result$gene_sets$pt_deg_symbols
  }
  
  data.frame(
    Trajectory = trajectory_label,
    Method = method,
    Common_tested_symbols = length(
      trajectory_result$gene_sets$universe_symbols
    ),
    Significant_gene_symbols = length(
      significant_symbols
    ),
    Mapped_background_Entrez = length(
      method_result$universe_entrez
    ),
    Mapped_significant_Entrez = length(
      method_result$deg_entrez
    ),
    Reactome_pathways_returned = nrow(
      method_result$table
    ),
    Significant_Reactome_pathways = nrow(
      method_result$significant_table
    ),
    stringsAsFactors = FALSE
  )
}

reactome_ora_summary <- bind_rows(
  
  make_ora_summary_row(
    trajectory_result = reactome_ora_results$traj1,
    trajectory_label = "Trajectory 1",
    method = "scFPC-DE"
  ),
  
  make_ora_summary_row(
    trajectory_result = reactome_ora_results$traj1,
    trajectory_label = "Trajectory 1",
    method = "PseudotimeDE"
  ),
  
  make_ora_summary_row(
    trajectory_result = reactome_ora_results$traj2,
    trajectory_label = "Trajectory 2",
    method = "scFPC-DE"
  ),
  
  make_ora_summary_row(
    trajectory_result = reactome_ora_results$traj2,
    trajectory_label = "Trajectory 2",
    method = "PseudotimeDE"
  )
)

print(
  reactome_ora_summary
)

## ============================================================
## Inspect significant pathways
## ============================================================

reactome_ora_results$traj1$scFPCDE$significant_table
reactome_ora_results$traj1$PseudotimeDE$significant_table

reactome_ora_results$traj2$scFPCDE$significant_table
reactome_ora_results$traj2$PseudotimeDE$significant_table

## ============================================================
## Save Reactome ORA results
## ============================================================

if (save_csv_files) {
  
  ## Complete ORA results
  
  write.csv(
    reactome_ora_results$traj1$scFPCDE$table,
    file = file.path(
      output_dir,
      "traj1_scFPCDE_Reactome_ORA_all.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    reactome_ora_results$traj1$PseudotimeDE$table,
    file = file.path(
      output_dir,
      "traj1_PseudotimeDE_Reactome_ORA_all.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    reactome_ora_results$traj2$scFPCDE$table,
    file = file.path(
      output_dir,
      "traj2_scFPCDE_Reactome_ORA_all.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    reactome_ora_results$traj2$PseudotimeDE$table,
    file = file.path(
      output_dir,
      "traj2_PseudotimeDE_Reactome_ORA_all.csv"
    ),
    row.names = FALSE
  )
  
  ## Significant pathways
  
  write.csv(
    reactome_ora_results$traj1$scFPCDE$significant_table,
    file = file.path(
      output_dir,
      "traj1_scFPCDE_Reactome_ORA_significant.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    reactome_ora_results$traj1$PseudotimeDE$significant_table,
    file = file.path(
      output_dir,
      "traj1_PseudotimeDE_Reactome_ORA_significant.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    reactome_ora_results$traj2$scFPCDE$significant_table,
    file = file.path(
      output_dir,
      "traj2_scFPCDE_Reactome_ORA_significant.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    reactome_ora_results$traj2$PseudotimeDE$significant_table,
    file = file.path(
      output_dir,
      "traj2_PseudotimeDE_Reactome_ORA_significant.csv"
    ),
    row.names = FALSE
  )
  
  ## Gene-symbol-to-Entrez mappings
  
  write.csv(
    reactome_ora_results$traj1$scFPCDE$universe_mapping,
    file = file.path(
      output_dir,
      "traj1_symbol_to_entrez_mapping.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    reactome_ora_results$traj2$scFPCDE$universe_mapping,
    file = file.path(
      output_dir,
      "traj2_symbol_to_entrez_mapping.csv"
    ),
    row.names = FALSE
  )
  
  ## Analysis summary
  
  write.csv(
    reactome_ora_summary,
    file = file.path(
      output_dir,
      "Reactome_ORA_summary.csv"
    ),
    row.names = FALSE
  )
}
