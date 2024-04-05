## rliger Next

- Standardized H5 IO specification that can be shared with other platforms.
  - Will move to use HDF5Array (TENxMatrix, H5ADMatrix)/ or BPCells for backed data representation.
  - Read feature metadata (e.g. id, name, ...) if available; Allow setting "id" as rownames, "name" for visualization.
  - rawData - coming from the original input, read only (qc filtering should be just stored in the object, no IO)
  - preprocessing metrics - nUMI, nGene and etc, still go "chunkApply" so the file is read only once
  - normData - delayed computed data from rawData, no on disk representation
  - scaleData - new on-disk file and then create object back, because RcppPlanc won't be able to handle delayed computation
- Ability to reorganize datasets
  - Allow doing something like `reorganize(ligerObj, variable = "somethingNotDataset")` and resulting in a new liger object with different ligerDataset grouping.
- Ability to do downstream analysis on H5 data
  - Pseudo-bulk should be easy because we are just aggregating cells.
  - Wilcoxon might be a bit harder because ranks are calculated per gene but the H5 sparse data is column majored. Might need to find a fast on-disk transposition method.

## rliger 2.0.1

- Fixed wrong UINMF aborting criteria
- Fixed example/test skipping criteria for non-existing dependencies
- Fixed file access issue when checking on CRAN
- Updated installed data file `system.file("extdata/ctrl.h5", "extdata/stim.h5")` to be of standard 10X H5 format
- Updated `quantileNorm()` automatic reference selection according to #297
- Other minor fixes (including #308)

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
- Added GO enrichment dowsntream analysis in `runGOEnrich()`
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

