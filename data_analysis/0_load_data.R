## ============================================================
## Resolve the HB6 analysis input
##
## Modes are selected with SCFPCDE_PAPER_DATA_MODE:
##   auto         - use local full objects when both exist; otherwise use the
##                  preprocessed paper dataset
##   local        - require data/scPure2_HB6_UMAP3D.rds and data/cds_HB6.rds
##   preprocessed - require data/scFPCDE_hb6.rda from this repository
##   package      - deprecated alias for preprocessed
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
preprocessed_data_file <- file.path(
  scFPCDE_paper_root,
  "data",
  "scFPCDE_hb6.rda"
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
allowed_modes <- c("auto", "local", "preprocessed", "package")
if (!requested_mode %in% allowed_modes) {
  stop(
    "SCFPCDE_PAPER_DATA_MODE must be one of: ",
    paste(allowed_modes, collapse = ", "),
    call. = FALSE
  )
}
if (requested_mode == "package") {
  warning(
    "SCFPCDE_PAPER_DATA_MODE=package is deprecated; using preprocessed data ",
    "from the paper-code repository.",
    call. = FALSE
  )
  requested_mode <- "preprocessed"
}

analysis_data_source <- if (requested_mode == "auto") {
  if (local_data_available) "local" else "preprocessed"
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
  hb6_preprocessed_data <- NULL
  analysis_data_source <- "local_full_objects"

  message(
    "HB6 source: local full Seurat and Monocle objects in ",
    file.path(scFPCDE_paper_root, "data"), "."
  )
} else {
  if (!file.exists(preprocessed_data_file)) {
    stop(
      "The preprocessed HB6 dataset was not found at: ",
      preprocessed_data_file,
      call. = FALSE
    )
  }

  hb6_data_environment <- new.env(parent = emptyenv())
  loaded_objects <- load(
    preprocessed_data_file,
    envir = hb6_data_environment
  )
  if (!"scFPCDE_hb6" %in% loaded_objects) {
    stop(
      "The preprocessed HB6 data file does not contain scFPCDE_hb6.",
      call. = FALSE
    )
  }

  hb6_preprocessed_data <- hb6_data_environment$scFPCDE_hb6
  required_fields <- c("yt", "counts", "tt", "clusters")
  invalid <- vapply(hb6_preprocessed_data, function(x) {
    !all(required_fields %in% names(x)) ||
      !identical(dimnames(x$counts), dimnames(x$yt))
  }, logical(1))
  if (any(invalid)) {
    stop(
      "The preprocessed HB6 dataset lacks required fields or aligned counts.",
      call. = FALSE
    )
  }

  seu <- NULL
  cds <- NULL
  mst <- NULL
  analysis_data_source <- "paper_preprocessed"

  message(
    "HB6 source: data/scFPCDE_hb6.rda preprocessed trajectory subsets."
  )
}
