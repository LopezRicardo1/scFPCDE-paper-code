## ============================================================
## Resolve the HB6 analysis input
##
## Modes are selected with SCFPCDE_PAPER_DATA_MODE:
##   auto    - use local full objects when both exist; otherwise package data
##   local   - require data/scPure2_HB6_UMAP3D.rds and data/cds_HB6.rds
##   package - require the preprocessed scFPCDE_hb6 package dataset
## ============================================================

suppressPackageStartupMessages(library(scFPCDE))

if (file.exists(file.path("data_analysis", "0_load_data.R"))) {
  scFPCDE_paper_root <- normalizePath(".", mustWork = TRUE)
} else if (file.exists("0_load_data.R")) {
  scFPCDE_paper_root <- normalizePath("..", mustWork = TRUE)
} else {
  stop(
    "Run this script from the repository root or data_analysis directory.",
    call. = FALSE
  )
}

local_data_files <- c(
  seurat = file.path(scFPCDE_paper_root, "data", "scPure2_HB6_UMAP3D.rds"),
  monocle = file.path(scFPCDE_paper_root, "data", "cds_HB6.rds")
)
local_data_available <- all(file.exists(local_data_files))
if (any(file.exists(local_data_files)) && !local_data_available) {
  warning(
    "Only one local HB6 input file is present; package mode will be used ",
    "unless local mode is explicitly requested.",
    call. = FALSE
  )
}

requested_mode <- tolower(Sys.getenv(
  "SCFPCDE_PAPER_DATA_MODE",
  unset = "auto"
))
allowed_modes <- c("auto", "local", "package")
if (!requested_mode %in% allowed_modes) {
  stop(
    "SCFPCDE_PAPER_DATA_MODE must be one of: ",
    paste(allowed_modes, collapse = ", "),
    call. = FALSE
  )
}

analysis_data_source <- if (requested_mode == "auto") {
  if (local_data_available) "local" else "package"
} else {
  requested_mode
}

if (analysis_data_source == "local") {
  if (!local_data_available) {
    stop(
      "Local mode requires both files:\n- ",
      paste(local_data_files, collapse = "\n- "),
      call. = FALSE
    )
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Local mode requires the SeuratObject package.", call. = FALSE)
  }

  seu <- readRDS(local_data_files[["seurat"]]) |>
    SeuratObject::UpdateSeuratObject()
  cds <- readRDS(local_data_files[["monocle"]])
  mst <- t(cds@principal_graph_aux$UMAP$pr_graph_cell_proj_dist)
  hb6_package_data <- NULL
  analysis_data_source <- "local_full_objects"

  message(
    "HB6 source: local full Seurat and Monocle objects in ",
    file.path(scFPCDE_paper_root, "data"), "."
  )
} else {
  utils::data(
    "scFPCDE_hb6",
    package = "scFPCDE",
    envir = environment()
  )
  if (!exists("scFPCDE_hb6", inherits = FALSE)) {
    stop(
      "The installed scFPCDE package does not provide scFPCDE_hb6. ",
      "Install the current package version and try again.",
      call. = FALSE
    )
  }

  hb6_package_data <- scFPCDE_hb6
  required_fields <- c("yt", "counts", "tt", "clusters")
  invalid <- vapply(hb6_package_data, function(x) {
    !all(required_fields %in% names(x)) ||
      !identical(dimnames(x$counts), dimnames(x$yt))
  }, logical(1))
  if (any(invalid)) {
    stop(
      "The installed scFPCDE_hb6 dataset lacks aligned raw counts. ",
      "Reinstall the current scFPCDE package before running paper code.",
      call. = FALSE
    )
  }

  seu <- NULL
  cds <- NULL
  mst <- NULL
  analysis_data_source <- "package_preprocessed"

  message(
    "HB6 source: scFPCDE::scFPCDE_hb6 preprocessed trajectory subsets."
  )
}
