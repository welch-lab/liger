## rliger2 1.9.9

- Added `ligerDataset` class for per-dataset information storage, with inheritance for specific modalities
- Added a number of plotting functions with clear function names and useful functionality
- Added Leiden clustering method, now as default rather than Louvain
- Added pseudo-bulk DEG method
- Added gene name pattern for expression percentage QC metric
- Added native Seurat object support for the core integration workflow
- Added a documentation website built with pkgdown
- Changed `liger` object class structure
- Changed iNMF (previously `optimizeALS()`), UINMF (previously `optimizeALS(unshared = TRUE)`) and online iNMF (previously `online_iNMF()`) implementation with vastly improved performance.
Now named by `runINMF()`, `runUINMF()` and `runOnlineINMF()` respectively, and wrapped in 
`runIntegration()`.
