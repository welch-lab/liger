## rliger Next

- Standardized H5 writing specification that can be shared with other platforms.
  - Currently we allow analysis with 10X cellranger output H5 file and H5AD file from anndata>=0.8.0
  - Writing to H5AD file should follow anndata specification otherwise the file cannot be read back to a Python seesion.
  - Writing to 10X H5 file should be carefully investigated.
- Ability to reorganize datasets
  - Allow doing something like `reorganize(ligerObj, variable = "somethingNotDataset")` and resulting in a new liger object with different ligerDataset grouping.
- Ability to do downstream analysis on H5 data
  - Pseudo-bulk should be easy because we are just aggregating cells.
  - Wilcoxon might be a bit harder because ranks are calculated per gene but the H5 sparse data is column majored. Might need to find a fast on-disk transposition method.

## rliger 2.0.0

- Added `ligerDataset` class for per-dataset information storage, with inheritance for specific modalities
- Added a number of plotting functions with clear function names and useful functionality
- Added Leiden clustering method, now as default rather than Louvain
- Added pseudo-bulk DEG method
- Added DEG analysis with one-vs-rest marker detection in `runMarkerDEG()` and pairwise comparison in `runPairwiseDEG()`
- Added gene name pattern for expression percentage QC metric
- Added native Seurat object support for the core integration workflow
- Added a documentation website built with pkgdown
- Added new iNMF variant method, consensus iNMF (c-iNMF), in `runCINMF()`. Not stable.
- Changed `liger` object class structure
- Moved iNMF (previously `optimizeALS()`), UINMF (previously `optimizeALS(unshared = TRUE)`) and online iNMF (previously `online_iNMF()`) implementation to new package *RcppPlanc* with vastly improved performance. Now wrapped in `runINMF()`, `runUINMF()` and `runOnlineINMF()` respectively, and all can be invoked with `runIntegration()`.
- Updated H5AD support to match up with Python anndata package 0.8.0 specs
- Renamed many function/argument names to follow camelCase style, original names are still available while deprecation warnings are issued

## rliger 1.0.1

- Allow setting mito pattern in `getMitoProportion()` #271
- Fix efficiency issue when taking the log of norm.data (e.g. `runWilcoxon`)
- Add runable examples to all exported functions when possible
- Fix typo in online_iNMF matrix initialization
- Adapt to Seurat5
- Other minor fixes

