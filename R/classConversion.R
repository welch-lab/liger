setClass("ligerDataset")
setClass("SingleCellExperiment")

setAs("ligerDataset", "SingleCellExperiment", function(from) {
    requireNamespace("SingleCellExperiment")
    assays <- list()
    if (!is.null(raw.data(from)))
        assays <- c(assays, list(counts = raw.data(from)))
    if (!is.null(norm.data(from)))
        assays <- c(assays, list(normcounts = norm.data(from)))
    SingleCellExperiment::SingleCellExperiment(
        assays = assays
    )
})



#' @importClassesFrom SeuratObject Seurat
setAs("Seurat", "ligerDataset", function(from) {
    requireNamespace("Seurat")
    raw.data <- Seurat::GetAssayData(from, slot = "counts")
})
