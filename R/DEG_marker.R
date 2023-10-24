#' @title Find DEG between two groups
#' @description Find DEG between two groups. Two methods are supported:
#' \code{"wilcoxon"} and \code{"pseudoBulk"}. Wilcoxon rank sum test is
#' performed on single-cell level, while pseudo-bulk method aggregates cells
#' basing on biological replicates and calls bulk RNAseq DE methods, DESeq2 wald
#' test. When real biological replicates are not available, pseudo replicates
#' can be generated. Please see below for detailed scenario usage.
#' @section Pairwise DEG Scenarios:
#' Users can select classes of cells from a variable in \code{cellMeta}.
#' \code{variable1} and \code{variable2} are used to specify a column in
#' \code{cellMeta}, and \code{group1} and \code{group2} are used to specify
#' existing classes from \code{variable1} and \code{variable2}, respectively.
#' When \code{variable2} is missing, \code{group2} will be considered from
#' \code{variable1}.
#'
#' For example, when \code{variable1 = "celltype"} and \code{variable2 = NULL},
#' \code{group1} and \code{group2} should be valid cell types in
#' \code{object$celltype}.
#'
#' When \code{variable1} is "celltype" and \code{variable2} is "gender",
#' \code{group1} should be a valid cell type from \code{object$celltype} and
#' \code{group2} should be a valid class from \code{object$gender}.
#'
#' When both \code{variable1} and \code{variable2} are missing, \code{group1}
#' and \code{group2} should be valid index of cells in \code{object}.
#' @param object A \linkS4class{liger} object, with normalized data available
#' @param group1,group2,variable1,variable2 Condition specification. See
#' \code{?runPairwiseDEG} section \bold{Pairwise DEG Scenarios} for detail.
#' @param method DEG test method to use. Choose from \code{"wilcoxon"} or
#' \code{"pseudoBulk"}. Default \code{"wilcoxon"}
#' @param usePeak Logical. Whether to use peak count instead of gene count.
#' Only supported when ATAC datasets are involved. Default \code{FALSE}.
#' @param useReplicate \code{cellMeta} variable of biological replicate
#' annotation. Only used with \code{method = "pseudoBulk"}. Default \code{NULL}
#' will create \code{nPsdRep} pseudo replicates per group.
#' @param nPsdRep Number of pseudo replicates to create. Only used when
#' \code{method = "pseudoBulk", useReplicate = NULL}. Default \code{5}.
#' @param seed Random seed to use for pseudo-replicate generation. Default
#' \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return A data.frame with DEG information
#' @rdname liger-DEG
#' @export
#' @examples
#' # Compare between cluster "0" and cluster "1"
#' degStats <- runPairwiseDEG(pbmcPlot, group1 = 0, group2 = 1,
#'                            variable1 = "leiden_cluster")
#' # Compare between all cells from cluster "5" and
#' # all cells from dataset "stim"
#' degStats <- runPairwiseDEG(pbmcPlot, group1 = "5", group2 = "stim",
#'                            variable1 = "leiden_cluster",
#'                            variable2 = "dataset")
runPairwiseDEG <- function(
        object,
        group1,
        group2,
        variable1 = NULL,
        variable2 = NULL,
        method = c("wilcoxon", "pseudoBulk"),
        usePeak = FALSE,
        useReplicate = NULL,
        nPsdRep = 5,
        seed = 1,
        verbose = getOption("ligerVerbose")
) {
    method <- match.arg(method)
    if (is.null(variable1) && is.null(variable2)) {
        # Directly using cell index
        groups <- list(
            .idxCheck(object, group1, "cell"),
            .idxCheck(object, group2, "cell")
        )
        group1Name <- "group1"
        group2Name <- "group2"
        names(groups) <- c("group1", "group2")
    } else if (!is.null(variable1)) {
        var1 <- .fetchCellMetaVar(object, variable1,
                                  checkCategorical = TRUE, drop = TRUE,
                                  droplevels = TRUE)
        group1Idx <- which(var1 %in% group1)
        group1Name <- paste(group1, collapse = ".")
        if (is.null(variable2)) {
            variable2 <- variable1
            var2 <- var1
        } else {
            var2 <- .fetchCellMetaVar(object, variable2,
                                           checkCategorical = TRUE, drop = TRUE,
                                           droplevels = TRUE)
        }
        group2Idx <- which(var2 %in% group2)
        group2Name <- paste(group2, collapse = ".")
        groups <- list(group1Idx, group2Idx)
        names(groups) <- c(group1Name, group2Name)
    } else {
        stop("Please see `?runPairwiseDEG` for usage.")
    }
    result <- .runDEG(object, groups = groups, method = method,
                      usePeak = usePeak, useReplicate = useReplicate,
                      nPsdRep = nPsdRep, seed = seed, verbose = verbose)
    result <- result[result$group == group1Name,]
    result$group <- NULL
    attributes(result)$meta <- list(
        group1 = group1,
        variable1 = variable1,
        group2 = group2,
        variable2 = variable2
    )
    return(result)
}

#' @rdname liger-DEG
#' @export
#' @param conditionBy \code{cellMeta} variable(s). Marker detection will be
#' performed for each level of this variable. Multiple variables will be
#' combined.
#' @param splitBy Split data by \code{cellMeta} variable(s) here and identify
#' markers for \code{conditionBy} within each chunk. Default \code{NULL}.
#' @param useDatasets Datasets to perform marker detection within. Default
#' \code{NULL} will use all datasets.
#' @section Marker Detection Scenarios:
#' Marker detection is generally performed in a one vs. rest manner. The
#' grouping of such condition is specified by \code{conditionBy}, which should
#' be a column name in \code{cellMeta}. When \code{splitBy} is specified as
#' another variable name in \code{cellMeta}, the marker detection will be
#' iteratively done for each level of \code{splitBy} variable.
#'
#' For example, when \code{conditionBy = "celltype"} and \code{splitBy = NULL},
#' marker detection will be performed by comparing all cells of "celltype_i"
#' against all other cells, and etc.
#'
#' When \code{conditionBy = "celltype"} and \code{splitBy = "gender"}, marker
#' detection will be performed by comparing "celltype_i" cells from "gender_i"
#' against other cells from "gender_i", and etc.
#' @examples
#' # Identify markers for each cluster
#' markerStats <- runMarkerDEG(pbmcPlot, conditionBy = "leiden_cluster")
#' # Identify dataset markers within each cluster
#' markerStatsList <- runMarkerDEG(pbmcPlot, conditionBy = "dataset",
#'                                 splitBy = "leiden_cluster")
runMarkerDEG <- function(
        object,
        conditionBy,
        splitBy = NULL, # The previous by dataset strategy
        method = c("wilcoxon", "pseudoBulk"),
        useDatasets = NULL,
        usePeak = FALSE,
        useReplicate = NULL,
        nPsdRep = 5,
        seed = 1,
        verbose = getOption("ligerVerbose")
) {
    useDatasets <- .checkUseDatasets(object, useDatasets)
    allCellIdx <- seq(ncol(object))[object$dataset %in% useDatasets]
    conditionBy <- .fetchCellMetaVar(
        object, conditionBy, cellIdx = allCellIdx,
        checkCategorical = TRUE, drop = FALSE, droplevels = TRUE
    )
    conditionBy <- interaction(conditionBy, drop = TRUE)
    splitBy <- .fetchCellMetaVar(
        object, splitBy, cellIdx = allCellIdx,
        checkCategorical = TRUE, drop = FALSE, droplevels = TRUE
    )
    splitBy <- interaction(splitBy, drop = TRUE)
    if (nlevels(splitBy) <= 1) {
        groups <- split(allCellIdx, conditionBy)
        result <- .runDEG(object, groups = groups, method = method,
                          usePeak = usePeak, useReplicate = useReplicate,
                          nPsdRep = nPsdRep, seed = seed, verbose = verbose)
    } else {
        result <- list()
        for (i in seq_along(levels(splitBy))) {
            subIdx <- splitBy == levels(splitBy)[i]
            subCellIdx <- allCellIdx[subIdx]
            groups <- split(subCellIdx, conditionBy[subIdx])
            result[[levels(splitBy)[i]]] <- .runDEG(
                object, groups = groups, method = method, usePeak = usePeak,
                useReplicate = useReplicate, nPsdRep = nPsdRep, seed = seed,
                verbose = verbose
            )
        }
    }

    return(result)
}
# groups - As underlying function, this must be organized into list of numeric
# cell index vectors.
.runDEG <- function(
        object,
        groups,
        method = c("wilcoxon", "pseudoBulk"),
        # byDataset = FALSE,
        usePeak = FALSE,
        useReplicate = NULL,
        nPsdRep = 5,
        seed = 1,
        verbose = getOption("ligerVerbose")
) {
    method <- match.arg(method)
    allCellIdx <- unlist(groups)
    allCellBC <- colnames(object)[allCellIdx]
    datasetInvolve <- levels(object$dataset[allCellIdx, drop = TRUE])
    var <- factor(unlist(lapply(names(groups), function(n) {
        rep(n, length(groups[[n]]))
    })))
    if (isTRUE(usePeak)) {
        useDatasets <- .checkUseDatasets(object, useDatasets = datasetInvolve,
                                         modal = "atac")
    } else {
        useDatasets <- .checkUseDatasets(object, useDatasets = datasetInvolve)
    }
    slot <- .DE.checkDataAvail(object, datasetInvolve, method, usePeak)
    dataList <- getMatrix(object, slot, datasetInvolve, returnList = TRUE)
    features <- Reduce(intersect, lapply(dataList, rownames))
    dataList <- lapply(dataList, function(x) x[features,])
    mat <- Reduce(cbind, dataList)
    mat <- mat[, allCellBC]
    if (method == "wilcoxon") {
        mat <- log1p(1e10*mat)
        result <- wilcoxauc(mat, var)
    } else if (method == "pseudoBulk") {
        if (is.null(useReplicate)) {
            replicateAnn <- setupPseudoRep2(var, nRep = nPsdRep,
                                            seed = seed)
        } else {
            replicateAnn <- .fetchCellMetaVar(
                object, useReplicate,
                cellIdx = allCellIdx,
                drop = FALSE,
                checkCategorical = TRUE,
                droplevels = TRUE
            )
            replicateAnn$groups <- var
        }
        pbs <- makePseudoBulk2(mat, replicateAnn, verbose = verbose)
        var <- sapply(levels(replicateAnn$groups), function(x) {
            nlevels(interaction(replicateAnn[replicateAnn$groups == x,],
                                drop = TRUE))
        })
        var <- factor(rep(names(var), var))
        result <- .callDESeq22(pbs, var, verbose)
    }
    return(result)
}

.DE.checkDataAvail <- function(object, useDatasets, method, usePeak) {
    if (isH5Liger(object, useDatasets)) {
        stop("HDF5 based datasets detected but is not supported. \n",
             "Try `object.sub <- downsample(object, useSlot = ",
             "'normData')` to create ANOTHER object with in memory data.")
    }
    if (method == "wilcoxon") {
        if (!isTRUE(usePeak)) slot <- "normData"
        else slot <- "normPeak"
    } else if (method == "pseudoBulk") {
        if (!requireNamespace("DESeq2", quietly = TRUE))
            stop("Package \"DESeq2\" needed for this function to work. ",
                 "Please install it by command:\n",
                 "BiocManager::install('DESeq2')",
                 call. = FALSE)
        if (!isTRUE(usePeak)) slot <- "rawData"
        else slot <- "rawPeak"
    }
    allAvail <- all(sapply(useDatasets, function(d) {
        ld <- dataset(object, d)
        !is.null(methods::slot(ld, slot))
    }))
    if (!allAvail)
        stop(slot, " not all available for involved datasets. [method = \"",
             method, "\", usePeak = ", usePeak, "]")
    return(slot)
}

setupPseudoRep2 <- function(groups, nRep = 3, seed = 1) {
    # The output data.frame should be cell per row by variable per col
    set.seed(seed)
    psdRep <- c()
    for (i in seq_along(levels(groups))) {
        groupSize <- sum(groups == levels(groups)[i])
        repVar <- sample(seq_len(groupSize) %% nRep) + 1 + (i - 1)*nRep
        psdRep <- c(psdRep, repVar)
    }
    return(data.frame(
        groups = groups,
        pseudoRep = factor(psdRep)
    ))
}

makePseudoBulk2 <- function(mat, replicateAnn, verbose = TRUE) {
    # mat - Extracted and contatenated matrix. intersection of genes by
    #       c(group1, group2) cells
    # groups - list of groups
    # replicateAnn - data.frame of replicate annotation, with rownames as
    #                barcodes and columns as variables

    # Check whether enough replicates per condition
    for (gr in levels(replicateAnn$groups)) {
        subrep <- replicateAnn[replicateAnn$groups == gr,]
        splitLabel <- interaction(subrep, drop = TRUE)
        if (length(levels(splitLabel)) < 2) {
            stop("Too few replicate labels for condition \"", gr, "\". ",
                 "Cannot not create pseudo-bulks. Please use ",
                 "`method = \"wilcoxon\"` instead.")
        }
    }
    splitLabel <- interaction(replicateAnn, drop = TRUE)
    pseudoBulks <- colAggregateSums_sparse(mat, as.integer(splitLabel) - 1,
                                           nlevels(splitLabel))
    dimnames(pseudoBulks) <- list(rownames(mat), levels(splitLabel))
    return(pseudoBulks)
}

.callDESeq22 <- function(pseudoBulks, groups,
                         verbose = getOption("ligerVerbose")) {
    # DESeq2 workflow
    if (isTRUE(verbose)) .log("Calling DESeq2 Wald test")

    des <- DESeq2::DESeqDataSetFromMatrix(
        pseudoBulks,
        colData = data.frame(groups = groups),
        design = ~groups
    )
    des <- DESeq2::DESeq(des, test = "Wald", quiet = !verbose)
    res <- DESeq2::results(des)
    # res <- as.data.frame(res)
    # res$feature <- rownames(res)
    # rownames(res) <- NULL
    # res <- res[, c(7, 2, 5, 6)]
    # colnames(res) <- c("feature", "logFC", "pval", "padj")
    # res <- res[!is.na(res$padj),]
    # res <- res[order(res$padj),]
    return(res)
}
