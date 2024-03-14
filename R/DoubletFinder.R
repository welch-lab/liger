#' Doublet detection with DoubletFinder
#' @description Detect doublet with DoubletFinder. Package "Seurat" and
#' "DoubletFinder" would be required to run this function.
#'
#' This wrapper runs Seurat PCA workflow (NormalizeData,
#' FindVariableFeatures, ScaleData, RunPCA) with all default settings on each
#' dataset, and then calls \code{DoubletFinder::doubletFinder}. Users that
#' prefer having more control on the preprocessing part might consider creating
#' single-sample Seurat object with
#' \code{CreateSeuratObject(rawData(object, "datasetName"))}.
#' @param object A \linkS4class{liger} object.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to run
#' \code{DoubletFinder::doubletFinder} with. Default \code{NULL}
#' applies to all datasets.
#' @param PCs Specific principal components to use. Default \code{1:10}.
#' @param nNeighbors Number of the PC neighborhood size used to compute pANN.
#' See "See Also". Scalar for all used datasets or vector for each. Default
#' \code{20}.
#' @param nExp The total number of doublet predictions produced. Scalar for all
#' used datasets or vector for each. Default \code{NULL} sets a 0.15 proportion.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param ... Additional arguments passed to
#' \code{DoubletFinder::doubletFinder}.
#' @return Updated \code{object} with variables \code{DoubletFinder_pANN} and
#' \code{DoubletFinder_classification} updated in \code{cellMeta} slot
#' @export
#' @examples
#' if (requireNamespace("DoubletFinder", quietly = TRUE)) {
#'     pbmc <- runDoubletFinder(pbmc)
#'     print(cellMeta(pbmc))
#' }
runDoubletFinder <- function(
        object,
        useDatasets = NULL,
        PCs = 1:10,
        nNeighbors = 20,
        nExp = NULL,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    if (!requireNamespace("DoubletFinder", quietly = TRUE)) { # nocov start
        cli::cli_abort(
            "Package {.pkg DoubletFinder} is needed for this function to work.
            Please install it by command:
            {.code remotes::install_github('DoubletFinder')}")
    }
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort(
            "Package {.pkg Seurat} is needed for this function to work.
            Please install it by command:
            {.code install.packages('Seurat')}")
    } # nocov end
    useDatasets <- .checkUseDatasets(object, useDatasets = useDatasets)
    nNeighbors <- .checkArgLen(nNeighbors, length(useDatasets), repN = TRUE, class = "numeric")
    nExp <- .checkArgLen(nExp, length(useDatasets), repN = TRUE, class = "numeric")
    if (is.null(nExp))
        nExp <- sapply(useDatasets, function(d)
            round(0.15 * ncol(dataset(object, d)))
        )
    object <- recordCommand(object, ...,
                            dependencies = c("Seurat", "DoubletFinder"))
    for (i in seq_along(useDatasets)) {
        d <- useDatasets[i]
        if (isTRUE(verbose)) cliID <- cli::cli_process_start("Running DoubletFinder on dataset {.val {d}}")
        seu <- Seurat::CreateSeuratObject(rawData(object, d)) %>%
            Seurat::NormalizeData(verbose = FALSE) %>%
            Seurat::FindVariableFeatures(verbose = FALSE) %>%
            Seurat::ScaleData(verbose = FALSE) %>%
            Seurat::RunPCA(verbose = FALSE)
        seu <- DoubletFinder::doubletFinder(seu, PCs = PCs, pK = nNeighbors[i],
                                            nExp = nExp[i], sct = FALSE, ...)
        seuMeta <- seu[[]]
        pANNCol <- grep(pattern = "pANN", colnames(seuMeta))
        DFCol <- grep(pattern = "DF.classifications", colnames(seuMeta))
        cellMeta(object, "DoubletFinder_pANN", useDatasets = d) <- seuMeta[,pANNCol]
        cellMeta(object, "DoubletFinder_classification", useDatasets = d) <- seuMeta[,DFCol]
        if (isTRUE(verbose)) cli::cli_process_done(id = cliID)
    }
    return(object)
}
