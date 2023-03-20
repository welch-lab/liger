#' Doublet detection with DoubletFinder
#' @description Detect doublet with DoubletFinder. Package "Seurat" and
#' "DoubletFinder" would be required to run this function. "DoubletFinder" can
#' only be installed via
#' \href{https://github.com/chris-mcginnis-ucsf/DoubletFinder}{GitHub}.
#'
#' Seurat default PCA workflow is performed prior to calling DoubletFinder.
#' @param object A \linkS4class{liger} object.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to run
#' \code{\link[DoubletFinder]{doubletFinder_v3}} with. Default \code{NULL}
#' applies to all datasets.
#' @param nPCs Number of top PC to use. Default \code{10}.
#' @param nNeighbors Number of the PC neighborhood size used to compute pANN.
#' See "See Also". Scalar for all used datasets or vector for each. Default
#' \code{20}.
#' @param nExp The total number of doublet predictions produced. Scalar for all
#' used datasets or vector for each. Default \code{NULL} sets a 0.15 proportion.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Additional arguments passed to
#' \code{\link[DoubletFinder]{doubletFinder_v3}}.
#' @return Updated \code{object} with variables \code{DoubletFinder_pANN} and
#' \code{DoubletFinder_classification} updated in \code{cellMeta} slot
#' @export
#' @seealso \code{\link[DoubletFinder]{doubletFinder_v3}}
#' @examples
#' pbmc <- runDoubletFinder(pbmc)
#' cellMeta(pbmc)
runDoubletFinder <- function(
        object,
        useDatasets = NULL,
        nPCs = 10,
        nNeighbors = 20,
        nExp = NULL,
        verbose = getOption("ligerVerbose"),
        ...
) {
    if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
        stop("DoubletFinder need to be installed")
    }
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat need to be installed")
    }
    useDatasets <- .checkUseDatasets(object, useDatasets = useDatasets)
    if (length(nNeighbors) == 1) {
        nNeighbors <- rep(nNeighbors, length(useDatasets))
    }
    if (length(nNeighbors) != length(useDatasets)) {
        stop("`nNeighbors` should be a single value for all used datasets ",
             "or an integer vector with number for each used dataset.")
    }
    if (!is.null(nExp)) {
        if (length(nExp) == 1) {
            nExp <- rep(nExp, length(useDatasets))
        }
        if (length(nExp) != length(useDatasets)) {
            stop("`nExp` should be a single value for all used datasets ",
                 "or an integer vector with number for each used dataset.")
        }
    } else {
        nExp <- sapply(useDatasets, function(d) {
            round(0.15 * ncol(dataset(object, d)))
        })
    }
    recordCommand(object, dependencies = c("Seurat", "DoubletFinder"))

    cellMeta(object, "DoubletFinder_pANN") <- NA
    cellMeta(object, "DoubletFinder_classification") <- NA
    for (i in seq_along(useDatasets)) {
        d <- useDatasets[i]
        if (isTRUE(verbose)) .log("Running DoubletFinder on dataset: ", d)
        ld <- dataset(object, d)
        seu <- Seurat::CreateSeuratObject(rawData(ld)) %>%
            Seurat::NormalizeData(verbose = FALSE) %>%
            Seurat::FindVariableFeatures(verbose = FALSE) %>%
            Seurat::ScaleData(verbose = FALSE) %>%
            Seurat::RunPCA(verbose = FALSE)
        suppressMessages(
            seu <- DoubletFinder::doubletFinder_v3(seu, PCs = seq(nPCs),
                                                   pK = nNeighbors[i],
                                                   nExp = nExp[i], ...)
        )
        seuMeta <- seu[[]]
        pANNCol <- grep(pattern = "pANN", colnames(seuMeta))
        DFCol <- grep(pattern = "DF.classifications", colnames(seuMeta))
        object$DoubletFinder_pANN[object$dataset == d] <- seuMeta[,pANNCol]
        object$DoubletFinder_classification[object$dataset == d] <-
            seuMeta[,DFCol]
    }
    object$DoubletFinder_classification <-
        factor(object$DoubletFinder_classification)
    return(object)
}
