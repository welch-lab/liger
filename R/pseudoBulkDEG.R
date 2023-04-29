#' Run DESeq2 differential expression test on pseudo-bulk expression
#' @description
#' This function aggregate raw counts data of specified comparison groups of
#' cells basing on either given biological replicate annotation or pseudo
#' replicates. According to user prior knowledge, existing replicate annotation
#' should be used when it is true that the cells involved in comparison come
#' from multiple biological replicates, while when all cells in both groups
#' come from the same replicate, pseudo replicates should be created.
#'
#' Known replicate annotation can be specified by \code{useCellMetaVar} to use
#' metadata stored in the \code{cellMeta} slot, or specify external annotation
#' through \code{replicateAnn}. To use pseudo replicates, users should leave the
#' two mentioned argument as default \code{NULL}.
#' @param object A \linkS4class{liger} object.
#' @param comparison A named list with two elements, where each is a vector of
#' valid cell specification (numeric or logical indexing should be based on the
#' whole object).
#' @param method Choose from \code{"DESeq2"}, \code{"limma"} or \code{"edgeR"}.
#' Default \code{"DESeq2"}
#' @param useCellMetaVar Name(s) of \code{cellMeta} variable that specify the
#' replicate annotation. Default \code{NULL}.
#' @param replicateAnn A data.frame or factor to pass external replicate
#' annotation. See Details for requirements. Default \code{NULL}.
#' @param nRep Number of pseudo replicates per comparison group, applied only
#' when running pseudo replicate mode (i.e. both \code{useCellMetaVar} and
#' \code{replicateAnn} are \code{NULL}). Default \code{5}.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @details
#' For using \code{replicateAnn}. When passing a data.frame, rownames have to be
#' available and it should cover all cells involved in both comparison group.
#' When using a named factor/vector, the names should as well cover all involved
#' cells. When using a un-named factor/vector, it is recommended to have its
#' length matched with the total number of cells of the whole object. When the
#' length of the un-named factor/vector equals to \code{length(comparison[[1]])
#' + length(comparison[[2]])}, it is assumed that the order matches to
#' \code{c(comparison[[1]], comparison[[2]])}.
#' @return A data.frame as the output of DESeq2. Ordered by adjusted p-value
#' (\code{padj}).
#' @export
#' @examples
#' rawData(datasets(pbmcPlot)[[1]]) <- rawData(dataset(pbmc, 1))
#' rawData(datasets(pbmcPlot)[[2]]) <- rawData(dataset(pbmc, 2))
#' res <- runPseudoBulkDEG(
#'     pbmcPlot,
#'     comparison = list(
#'         c1 = pbmcPlot$leiden_cluster == 1,
#'         c2 = pbmcPlot$leiden_cluster == 2
#'     )
#' )
#' print(head(res))
runPseudoBulkDEG <- function(
        object,
        comparison,
        method = c("DESeq2", "limma", "edgeR"),
        useCellMetaVar = NULL,
        replicateAnn = NULL,
        nRep = 5,
        seed = 1,
        verbose = getOption("ligerVerbose")
) {
    if (!inherits(object, "liger")) {
        stop("Please use a `liger` object.")
    }
    method <- match.arg(method)
    comparison <- .setupComparison(object, comparison)

    set.seed(seed)
    if (is.null(replicateAnn) && is.null(useCellMetaVar)) {
        if (isTRUE(verbose)) {
            .log("No replicate annotation specified. ",
                 "Creating pseudo-replicates...")
        }
        replicateAnn <- setupPseudoRep(object, comparison, nRep = nRep)
    } else {
        replicateAnn <- .checkReplicateAnn(object, comparison, useCellMetaVar,
                                           replicateAnn)
    }
    pbs <- makePseudoBulk(object, comparison, replicateAnn, verbose = verbose)
    comparison <- factor(rep(names(pbs), sapply(pbs, ncol)))

    switch(EXPR = method,
        DESeq2 = .callDESeq2(pbs, comparison, verbose),
        limma = stop("Not implemented yet")
    )

}

makePseudoBulk <- function(object, comparison, replicateAnn, verbose = TRUE) {
    # Fetch and order the raw data of involved cells
    allCellIdx <- Reduce(c, comparison)
    rawSubset <- subsetLiger(object, cellIdx = allCellIdx,
                             useSlot = "rawData",
                             returnObject = FALSE, verbose = verbose)
    rawList <- lapply(rawSubset, function(dlist) dlist$rawData)
    rawdata <- mergeSparseAll(rawList)
    #rawdata <- rawdata[, rownames(replicateAnn)]
    aggList <- list()
    # Check whether enough replicates per condition
    useSingleCell <- FALSE
    splitRawData <- list()
    splitLabels <- list()
    for (gr in names(comparison)) {
        cellInvolved <- rownames(replicateAnn)[replicateAnn$comparison == gr]
        splitRawData[[gr]] <- rawdata[, cellInvolved]
        subrep <- replicateAnn[replicateAnn$comparison == gr,]
        # Aggregate by replicate annotation
        splitLabel <- interaction(subrep, drop = TRUE)
        if (length(levels(splitLabel)) < 2) {
            warning("Too few replicate labels for condition \"", gr, "\". ",
                    "will not create pseudo-bulks but test at single cell ",
                    "level. `runWilcoxon()` is recommended instead.",
                    immediate. = TRUE)
            useSingleCell <- TRUE
        }
        splitLabels[[gr]] <- splitLabel
    }
    if (isTRUE(useSingleCell)) {
        return(splitRawData)
    }
    for (gr in names(comparison)) {
        if (isTRUE(verbose)) {
            .log("Generating pseudo-bulks for condition \"", gr, "\":\n",
                 paste(levels(splitLabel), collapse = ", "))
        }
        splitLabel <- splitLabels[[gr]]
        agg <- lapply(levels(splitLabel), function(rep) {
            Matrix::rowSums(splitRawData[[gr]][, splitLabel == rep])
        })
        names(agg) <- levels(splitLabel)
        agg <- Matrix::Matrix(as.matrix(data.frame(agg)), sparse = TRUE)
        rownames(agg) <- rownames(splitRawData[[gr]])
        aggList[[gr]] <- agg
    }
    return(aggList)
}

setupPseudoRep <- function(object, comparison, nRep = 3) {
    # The output data.frame should be cell per row by variable per col
    # Cells should be ordered by c(comparison[[1]], comparison[[2]])
    allCellIdx <- Reduce(c, comparison)
    allCellBarc <- colnames(object)[allCellIdx]
    grVars <- lapply(seq_along(comparison), function(i) {
        # "+1" for converting 0-base number to 1-base
        repVar <- sample(seq_along(comparison[[i]]) %% nRep) + 1 + (i - 1)*nRep
        cpVar <- rep(names(comparison)[i], length(comparison[[i]]))
        data.frame(comparison = cpVar, pseudoRep = repVar)
    })
    comparison <- Reduce(rbind, grVars)
    rownames(comparison) <- allCellBarc
    return(comparison)
}

.setupComparison <- function(object, comparison) {
    if (!is.list(comparison)) stop("Please use a named list for `comparison`")
    if (length(comparison) < 2) {
        stop("Please use at least 2 elements in `comparison` list")
    }
    cpnames <- names(comparison)
    comparison <- lapply(comparison, .idxCheck,
                         object = object, orient = "cell")
    if (!is.null(cpnames)) names(comparison) <- cpnames
    else names(comparison) <- paste0("group", seq_along(comparison))
    return(comparison)
}

.checkReplicateAnn <- function(object, comparison, useCellMetaVar = NULL,
                               replicateAnn = NULL) {
    # The output data.frame should be cell per row by variable per col
    # Cells should be ordered by c(comparison[[1]], comparison[[2]])
    allCellIdx <- Reduce(c, comparison)
    allCellBarc <- colnames(object)[allCellIdx]

    cpVar <- factor(Reduce(c, lapply(seq_along(comparison), function(i) {
        rep(names(comparison)[i], length(comparison[[i]]))
    })))
    comparison <- data.frame(comparison = cpVar, row.names = allCellBarc)
    useCellMetaVar <- .fetchCellMetaVar(object, useCellMetaVar,
                                        cellIdx = allCellIdx, drop = FALSE,
                                        checkCategorical = TRUE,
                                        droplevels = TRUE)
    if (is.null(useCellMetaVar)) {
        useCellMetaVar <- data.frame(row.names = allCellBarc)
    }

    if (is.null(replicateAnn)) {
        replicateAnn <- data.frame(row.names = allCellBarc)
    } else if (is.data.frame(replicateAnn)) {
        if (!all(allCellBarc %in% rownames(replicateAnn))) {
            stop("Not all cells involved in `comparison` are annotated in ",
                 "`replicateAnn`. Missing cells: ",
                 .nfstr(allCellBarc, rownames(replicateAnn)))
        }
        replicateAnn <- replicateAnn[allCellBarc, , drop = FALSE]
    } else if (is.factor(replicateAnn)) {
        if (!is.null(names(replicateAnn))) {
            # When names available, select basing on names,
            # can ignore the length
            if (!all(allCellBarc %in% names(replicateAnn))) {
                stop("Not all cells involved in `comparison` are annotated",
                     " in `replicateAnn`. Missing cells: ",
                     .nfstr(allCellBarc, names(replicateAnn)))
            }
            replicateAnn <- replicateAnn[allCellBarc]
        } else {
            if (length(replicateAnn) != ncol(object) &&
                length(replicateAnn) != length(allCellIdx)) {
                stop("Unable to format replicate annotation with given
                     vector/factor of length ", length(replicateAnn), ". ",
                     "See ?runPseudoBulkDEG")
            } else if (length(replicateAnn) == ncol(object)) {
                # Using all cell annotation, need to subset for involved cells
                replicateAnn <- replicateAnn[allCellIdx]
            }
        }

        replicateAnn <- droplevels(replicateAnn)
        replicateAnn <- data.frame(replicateAnn, row.names = allCellBarc)
    } else {
        stop("Please use a `data.frame` or a `factor` to specify replicate ",
             "annotation")
    }
    replicateAnn <- cbind(comparison, useCellMetaVar, replicateAnn)
    return(replicateAnn)
}

.callDESeq2 <- function(pseudoBulks, comparison,
                        verbose = getOption("ligerVerbose")) {
    if (!requireNamespace("DESeq2", quietly = TRUE))
        stop("Package \"DESeq2\" needed for this function to work. ",
             "Please install it by command:\n",
             "BiocManager::install('DESeq2')",
             call. = FALSE)
    # DESeq2 workflow
    if (isTRUE(verbose)) .log("Calling DESeq2 Wald test")
    des <- DESeq2::DESeqDataSetFromMatrix(
        mergeSparseAll(pseudoBulks),
        colData = data.frame(comparison = comparison),
        design = ~comparison
    )
    des <- DESeq2::DESeq(des, test = "Wald", quiet = !verbose)
    res <- DESeq2::results(des)
    res <- as.data.frame(res)
    res <- res[order(res$padj),]
    return(res)
}
