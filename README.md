# scFPCDE-paper-code

This repository contains the analysis code used to reproduce results from the manuscript *scFPC-DE: Differential Expression Analysis Along Single-Cell Trajectories via Functional Principal Component Analysis*.

## Repository structure

- `data_analysis/` contains the HB6 B-cell pseudotime reconstruction and the numbered real-data analysis scripts. Run the numbered scripts in order.
- `simulation/` contains the simulation-study code extracted from the supplementary material.

## HB6 input modes

Install the current `scFPCDE` package before running the numbered scripts. The
HB6 preprocessed data are stored in this repository, keeping the R package
focused on its reusable methods and simulation example. The loader chooses
between two reproducible input modes:

1. **Local full objects.** If both `data/scPure2_HB6_UMAP3D.rds` and
   `data/cds_HB6.rds` exist, `0_load_data.R` loads them and script 1 reproduces
   the cell-quality, pseudotime, variable-gene, and detection-filtering steps.
2. **Repository preprocessed data.** If those files are absent,
   `0_load_data.R` loads `data/scFPCDE_hb6.rda`. This contains the exact two
   already-subsetted paper trajectories, with aligned raw counts, uncentered
   logcounts, pseudotime, and B-cell labels. Script 1 starts from those inputs
   and reruns scFPCDE and PseudotimeDE.

Run from the repository root:

```r
source("data_analysis/0_load_data.R")
source("data_analysis/1_trajs_analysis_workflow.R")
```

Then run scripts 2 through 12 in numerical order. Script 1 can also be sourced
directly; it will invoke the loader when needed.

To choose a mode explicitly before starting R:

```sh
SCFPCDE_PAPER_DATA_MODE=preprocessed R
SCFPCDE_PAPER_DATA_MODE=local R
```

`auto` is the default. Preprocessed mode begins after trajectory reconstruction
and cell/gene selection, so `bcell_hb6_pseudotime.R` is needed only when
rebuilding the full Monocle trajectory from the original data. All required R
packages must still be installed before running the analyses.
