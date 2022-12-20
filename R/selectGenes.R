#' Select a subset of informative genes
#' @description This function identifies highly variable genes from each dataset
#' and combines these gene sets (either by union or intersection) for use in
#' downstream analysis. Assuming that gene expression approximately follows a
#' Poisson distribution, this function identifies genes with gene expression
#' variance above a given variance threshold (relative to mean gene expression).
#' @param object ligerObject
selectGenes <- function(
    object,
    var.thresh = 0.1,
    alpha.thresh = 0.99,
    num.genes = NULL,
    tol = 0.0001,
    capitalize = FALSE,
    combine = c("intersection", "union"),
    verbose = TRUE,
    chunkSize = 1000,
    ...
) {
    combine <- match.arg(combine)
    if (length(var.thresh) != 1 & length(var.thresh) != length(object)) {
        stop("Wrong length of `var.thresh`. Use 1 or `length(object)` values.")
    }
    if (length(var.thresh) == 1) var.thresh <- rep(var.thresh, length(object))
    gene.use <- list()
    for (d in names(object)) {
        i <- names(object) == d
        if (isTRUE(verbose))
            message(date(), " ... Selecting HVG for dataset: ", d)
        ld <- dataset(object, d)
        # Make sure that all required feature meta values exist
        if (!"geneVars" %in% names(feature.meta(ld))) {
            if (isH5Liger(ld)) {
                ld <- calcGeneVars.H5(ld, chunkSize = chunkSize,
                                      verbose = verbose)
            } else {
                feature.meta(ld)$geneMeans <- rowMeansFast(norm.data(ld))
                feature.meta(ld)$geneVars <-
                    rowVarsFast(norm.data(ld), feature.meta(ld)$geneMeans)
            }
            datasets(object)[[d]] <- ld
        }
        geneMeans <- feature.meta(ld)$geneMeans
        geneVars <- feature.meta(ld)$geneVars
        trx_per_cell <- cell.meta(object)[object$dataset == d, "nUMI"]
        nolan_constant <- mean((1 / trx_per_cell))
        alphathresh.corrected <- alpha.thresh / nrow(ld)
        geneMeanUpper <- geneMeans + qnorm(1 - alphathresh.corrected / 2) *
            sqrt(geneMeans * nolan_constant / ncol(ld))
        basegenelower <- log10(geneMeans * nolan_constant)
        num_varGenes <- function(x, num.genes.des) {
            # This function returns the difference between the desired number of
            # genes and the number actually obtained when thresholded on x
            y <- length(
                which(
                    geneVars / nolan_constant > geneMeanUpper &
                        log10(geneVars) > basegenelower + x
                )
            )
            return(abs(num.genes.des - y))
        }
        if (!is.null(num.genes)) {
            # Optimize to find value of x which gives the desired number of
            # genes for this dataset if very small number of genes requested,
            # `var.thresh` may need to exceed 1
            optimized <- stats::optimize(num_varGenes, c(0, 1.5), tol = tol,
                                         num.genes.des = num.genes[i])
            var.thresh[i] <- optimized$minimum
            if (is.na(optimized$objective)) {
                warning("Cannot optimize the number of selected genes for ",
                        "dataset \"", d, "\"")
            } else
            if (optimized$objective > 1) {
                warning("Returned number of genes for dataset ", d,
                        " differs from requested by ", optimized$objective,
                        ". Lower `tol` or `alpha.thresh` for better results.")
            }
        }
        selected <- rownames(ld)[
            geneVars / nolan_constant > geneMeanUpper &
                log10(geneVars) > basegenelower + var.thresh[i]
        ]
        if (isTRUE(verbose))
            message(date(), " ...   ", length(selected), " features selected")
        gene.use[[d]] <- selected
    }
    if (combine == "intersection") gene.use <- Reduce(intersect, gene.use)
    else gene.use <- Reduce(union, gene.use)

    for (d in names(object)) {
        gene.use <-
            gene.use[gene.use %in% rownames(dataset(object, d))]
    }

    if (length(gene.use) == 0) {
        warning("No genes were selected. Lower `var.thresh` values or choose",
                " \"union\" for `combine`", immediate. = TRUE)
    }
    object@var.features <- gene.use
    return(object)
}


#' Calculate Gene Variance for ligerDataset object
#' @param object ligerDataset object
#' @param chunkSize Integer for the maximum number of cells in each chunk.
#' Default \code{1000}.
#' @param verbose Logical. Whether to show a progress bar. Default \code{TRUE}.
#' @return The input \code{object} with calculated var updated in the H5 file.
#' @export
#' @useDynLib rliger, .registration = TRUE
#' @importFrom Rcpp sourceCpp
calcGeneVars.H5 <- function(object, chunkSize = 1000, verbose = TRUE) {
    h5file <- getH5File(object)
    if (h5file$exists("gene_vars")) {
        geneVars <- h5file[["gene_vars"]][]
    } else {
        geneVars <- rep(0, nrow(object))
        geneMeans <- h5file[["gene_means"]][]
        safeH5Create(
            object = object,
            dataPath = "gene_vars",
            dims = nrow(object),
            dtype = "double"
        )
        geneVars <- H5Apply(
            object,
            function(chunk, sparseXIdx, cellIdx, values) {
                values + sumSquaredDeviations(chunk, geneMeans)
            },
            init = geneVars,
            useData = "norm.data",
            chunkSize = chunkSize,
            verbose = verbose
        )
        geneVars <- geneVars / (ncol(object) - 1)
        h5file[["gene_vars"]][seq_along(geneVars)] <- geneVars
    }
    feature.meta(object)$geneVars <- geneVars
    object
}

