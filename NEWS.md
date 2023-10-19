## Next release (rliger 1.9.9)

### Added

- `ligerDataset` class for per-dataset information storage, with inheritance for specific modalities
- A number of plotting functions with clear function names and useful functionality
- Leiden clustering method, now as default rather than Louvain
- Pseudo-bulk DEG method
- Gene name pattern for expression percentage QC metric
- Native Seurat object support for the core integration workflow
- A documentation website built with pkgdown

### Changed

- `liger` object class structure
- iNMF (previously `optimizeALS()`), UINMF (previously `optimizeALS(unshared = TRUE)`) and online iNMF (previously `online_iNMF()`) implementation with vastly improved performance.
Now named by `runINMF()`, `runUINMF()` and `runOnlineINMF()` respectively, and wrapped in 
`runIntegration()`.
