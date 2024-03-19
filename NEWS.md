## rliger 2.0.0

- Added `ligerDataset` class for per-dataset information storage, with inheritance for specific modalities
- Added a number of plotting functions with clear function names and useful functionality
- Added Leiden clustering method, now as default rather than Louvain
- Added pseudo-bulk DEG method
- Added DEG analysis with one-vs-rest marker detection in `runMarkerDEG()` and pairwise comparison in `runPairwiseDEG()`
- Added gene name pattern for expression percentage QC metric
- Added native Seurat object support for the core integration workflow
- Added a documentation website built with pkgdown
- Added new iNMF variant method, consensus iNMF (c-iNMF), in `runCINMF()`
- Changed `liger` object class structure
- Changed iNMF (previously `optimizeALS()`), UINMF (previously `optimizeALS(unshared = TRUE)`) and online iNMF (previously `online_iNMF()`) implementation with vastly improved performance. Now named by `runINMF()`, `runUINMF()` and `runOnlineINMF()` respectively, and wrapped in `runIntegration()`.
- Updated H5AD support to match up with Python anndata package 0.8.0 specs

## rliger 1.0.1

- Allow setting mito pattern in `getMitoProportion()` #271
- Fix efficiency issue when taking the log of norm.data (e.g. `runWilcoxon`)
- Add runable examples to all exported functions when possible
- Fix typo in online_iNMF matrix initialization
- Adapt to Seurat5
- Other minor fixes

