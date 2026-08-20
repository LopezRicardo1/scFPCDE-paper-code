# HB6 analysis data

`scFPCDE_hb6.rda` contains the two preprocessed HB6 B-cell trajectories used
by the numbered manuscript-analysis scripts. Each trajectory contains aligned
cell-by-gene raw counts and uncentered log-normalized expression, pseudotime,
B-cell state labels, cell and gene identifiers, and analysis provenance.

The object is intentionally stored in the paper-code repository rather than in
the `scFPCDE` R package. The original Seurat and Monocle objects are not bundled.
When both full objects are placed in this directory, `0_load_data.R` can instead
reconstruct the analysis inputs from them by using local mode.

SHA-256 for the distributed file:

```text
e8eb71ceac203bf18e24a534d8d048975b41837ad7e724ca669fb4fd56cf969c  scFPCDE_hb6.rda
```
