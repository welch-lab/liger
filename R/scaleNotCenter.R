#' Scale genes by root-mean-square across cells
#' @description This function scales normalized gene expression data after
#' variable genes have been selected. Note that the data is not mean-centered
#' before scaling because expression values must remain positive (NMF only
#' accepts positive values).
#' @param object \linkS4class{liger} object, \code{link{selectGenes}} should
#' have been run in advance.
#' @param chunk Integer. Number of maximum number of cells in each chunk, when
#' scaling is applied to any HDF5 based dataset. Default \code{1000}.
#' @param verbose Logical. Whether to show information of the progress.
#' @return Updated \code{object}, where the \code{scale.data} slot of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot is updated.
#' @export
scaleNotCenter <- function(
    object,
    chunk = 1000,
    verbose = TRUE
) {
    if (is.null(var.features(object)) ||
        length(var.features(object)) == 0) {
        stop("No variable feature found. ",
             "Please check the result of `selectGenes()`")
    }
    for (d in names(object)) {
        if (isTRUE(verbose)) message(date(), " ... Scaling dataset: ", d)
        ld <- dataset(object, d)
        if (isH5Liger(ld)) {
            # Scale H5 based data
            featureIdx <- rownames(ld) %in% var.features(object)
            features <- rownames(ld)[featureIdx]
            geneSumSq <- feature.meta(ld)$geneSumSq[featureIdx]
            nCells <- ncol(ld)
            geneRootMeanSumSq = sqrt(geneSumSq / (nCells - 1))
            h5file <- getH5File(ld)
            resultH5Path <- "scale.data"
            safeH5Create(
                ld,
                dataPath = resultH5Path,
                dims = c(length(features), nCells),
                dtype = "double",
                chunkSize = c(length(features), chunk)
            )
            H5Apply(
                ld,
                useData = "norm.data",
                chunkSize = chunk,
                verbose = verbose,
                FUN = function(chunk, sparseXIdx, cellIdx, values) {
                    chunk <- chunk[featureIdx,]
                    chunk = as.matrix(chunk)
                    chunk = sweep(chunk, 1, geneRootMeanSumSq, "/")
                    rownames(chunk) <- features
                    chunk = chunk[var.features(object),]
                    chunk[is.na(chunk)] = 0
                    chunk[chunk == Inf] = 0
                    h5file[[resultH5Path]][seq_along(features),
                                           cellIdx] <- chunk
                }
            )
            h5file.info(ld, "scale.data", check = FALSE) <- resultH5Path
            safeH5Create(
                ld,
                dataPath = "scale.data.featureIdx",
                dims = length(features),
                dtype = "double"
            )
            h5file[["scale.data.featureIdx"]][1:length(features)] <-
                which(featureIdx)
        } else {
            # Scale in memory data
            # IMPORTANT: scale data is cell (row) by gene (column)
            norm.subset <- t(norm.data(ld)[var.features(object),])
            scaled <- scaleNotCenterFast(norm.subset)
            scaled <- as.matrix(t(scaled))
            scaled[is.na(scaled)] <- 0
            scaled[scaled == Inf] = 0
            rownames(scaled) <- var.features(object)
            colnames(scaled) <- colnames(ld)
            scale.data(ld, check = FALSE) <- scaled
        }
        datasets(object)[[d]] <- ld
    }
    object
}
