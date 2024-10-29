#' @title Find DEG between groups
#' @description Two methods are supported: \code{"pseudoBulk"} and
#' \code{"wilcoxon"}. Pseudo-bulk method aggregates cells basing on biological
#' replicates and calls bulk RNAseq DE methods, DESeq2 wald test, while
#' Wilcoxon rank sum test is performed on single-cell level.
#' \code{runPairwiseDEG()} is generally used for flexibly comparing two specific
#' groups of cells, while \code{runMarkerDEG()} is used for a one-vs-rest marker
#' test strategy.
#'
#' While using pseudo-bulk method, it is generally recommended that you have
#' these variables available in your object:
#'
#' \enumerate{
#'   \item{The cell type or cluster labeling. This can be obtained from prior
#'   study or computed with \code{\link{runCluster}}}
#'   \item{The biological replicate labeling, most of the time the
#'   \code{"dataset"} variable automatically generated when the
#'   \linkS4class{liger} object is created. Users may use other variables if
#'   a "dataset" is merged from multiple replicates.}
#'   \item{The condition labeling that reflects the study design, such as the
#'   treatment or disease status for each sample/dataset.}
#' }
#'
#' Please see below for detailed scenarios.
#'
#' @section Using Wilcoxon rank-sum test:
#' Wilcoxon rank-sum test works for each gene and is based on the rank of the
#' expression in each cell. LIGER provides dataset integration but does not
#' "correct" the expression values. Projects with strong batch effects or
#' integrate drastically different modalities should be cautious when using
#' this method.
#'
#' @section Comparing difference between/across cell types:
#' Most of times, people would want to know what cell types are for each cluster
#' after clustering. This can be done with a marker detection method that test
#' each cluster against all the other cells. This can be done with a command
#' like \code{runMarkerDEG(object, conditionBy = "cluster_var")}. When using
#' default pseudo-bulk method, users should additionaly determine the
#' pseudo-bulk setup parameters. If the real biological replicate variable is
#' available, it should be supplied to argument \code{useReplicate}, otherwise,
#' pseudo-replicates should be created. See "Pseudo-Replicate" section for more.
#'
#' @section Compare between conditions:
#' It is frequently needed to identify the difference between conditions. Users
#' can simply set \code{conditionBy = "condition_var"}. However, most of time,
#' such comparisons should be ideally done in a per-cluster manner. This can be
#' done by setting \code{splitBy = "cluster_var"}. This will run a loop for each
#' cluster, and within the group of cells, compare each condition against all
#' other cells in the cluster.
#'
#' In the scenario when users only need to compare two conditions for each
#' cluster, running \code{runPairwiseDEG(object, groupTest = "condition1",
#' groupCtrl = "condition2", variable1 = "condition_var",
#' splitBy = "cluster_var")} would address the need.
#'
#' For both use case, if pseudo-bulk (default) method is used, users should
#' determine the pseudo-bulk setup parameters as mentioned in the previous
#' section.
#'
#' @section Detailed \code{runMarkerDEG} usage:
#' Marker detection is performed in a one vs. rest manner. The grouping of such
#' condition is specified by \code{conditionBy}, which should be a column name
#' in \code{cellMeta}. When \code{splitBy} is specified as another variable
#' name in \code{cellMeta}, the marker detection will be iteratively done for
#' within each level of \code{splitBy} variable.
#'
#' For example, when \code{conditionBy = "celltype"} and \code{splitBy = NULL},
#' marker detection will be performed by comparing all cells of "celltype_i"
#' against all other cells, and etc. This is analogous to the old version when
#' running \code{runWilcoxon(method = "cluster")}.
#'
#' When \code{conditionBy = "gender"} and \code{splitBy = "leiden_cluster"},
#' marker detection will be performed by comparing "gender_i" cells from "cluster_j"
#' against other cells from "cluster_j", and etc. This is analogous to the old
#' version when running \code{runWilcoxon(method = "dataset")}.
#'
#' @section Detailed \code{runPairwiseDEG} usage:
#' Users can select classes of cells from a variable in \code{cellMeta}.
#' \code{variable1} and \code{variable2} are used to specify a column in
#' \code{cellMeta}, and \code{groupTest} and \code{groupCtrl} are used to specify
#' existing classes from \code{variable1} and \code{variable2}, respectively.
#' When \code{variable2} is missing, \code{groupCtrl} will be considered from
#' \code{variable1}.
#'
#' For example, when \code{variable1 = "celltype"} and \code{variable2 = NULL},
#' \code{groupTest} and \code{groupCtrl} should be valid cell types in
#' \code{object$celltype}.
#'
#' When \code{variable1} is "celltype" and \code{variable2} is "gender",
#' \code{groupTest} should be a valid cell type from \code{object$celltype} and
#' \code{groupCtrl} should be a valid class from \code{object$gender}.
#'
#' When both \code{variable1} and \code{variable2} are missing, \code{groupTest}
#' and \code{groupCtrl} should be valid index of cells in \code{object}.
#'
#' @section Pseudo-Replicate:
#' Pseudo-replicate assignment is a technique to complement the lack of real
#' biological replicates when using pseudo-bulk DE methods. LIGER's pseudo-bulk
#' method generally requires that each comparison group has at least 3
#' replicates each composed of at least 3 cells, in order to ensure the
#' statistic power. When less than 3 real replicates are found for a comparison,
#' the default setting (\code{nPsdRep = NULL}) splits each into 3
#' pseudo-replicates, otherwise no pseudo-replicates are automatically
#' generated. When \code{nPsdRep} is given a number, LIGER will always go
#' through each comparison group and split each real replicate into the given
#' number of pseudo-replicates.
#'
#' @param object A \linkS4class{liger} object, with normalized data available
#' @param groupTest,groupCtrl,variable1,variable2 Condition specification. See
#' \code{?runPairwiseDEG} section \bold{Pairwise DEG Scenarios} for detail.
#' @param splitBy Name(s) of the variable(s) in \code{cellMeta} to split the
#' comparison. See Details. Default \code{NULL}.
#' @param method DEG test method to use. Choose from \code{"pseudoBulk"} or
#' \code{"wilcoxon"}. Default \code{"pseudoBulk"}
#' @param usePeak Logical. Whether to use peak count instead of gene count.
#' Only supported when ATAC datasets are involved. Default \code{FALSE}.
#' @param useReplicate \code{cellMeta} variable of biological replicate
#' annotation. Only used with \code{method = "pseudoBulk"}. Default
#' \code{"dataset"}.
#' @param nPsdRep Number of pseudo-replicates to create. Only used when
#' \code{method = "pseudoBulk"}. Default \code{NULL}. See Details.
#' @param minCellPerRep Numeric, will not make pseudo-bulk for replicate with
#' less than this number of cells. Default \code{3}.
#' @param printDiagnostic Logical. Whether to show more detail when
#' \code{verbose = TRUE}. Default \code{FALSE}.
#' @param chunk Number of features to process at a time during Wilcoxon test.
#' Useful when memory is limited. Default \code{NULL} will process all features
#' at once.
#' @param seed Random seed to use for pseudo-replicate generation. Default
#' \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @return A data.frame with DEG information with the all or some of the
#' following fields:
#'  \item{feature}{Gene names}
#'  \item{group}{Test group name. Multiple tests might be present for each
#'    function call. This is the main variable to distinguish the tests. For a
#'    pairwise test, a row with a certain group name represents the test result
#'    between the this group against the other control group; When split by a
#'    variable, it would be presented in "split.group" format, meaning the stats
#'    is by comparing the group in the split level against the control group in
#'    the same split level. When running marker detection without splitting,
#'    a row with group "a" represents the stats of the gene in group "a" against
#'    all other cells. When running split marker detection, the group name would
#'    be in "split.group" format, meaning the stats is by comparing the group in
#'    the split level against all other cells in the same split level.}
#'  \item{logFC}{Log fold change}
#'  \item{pval}{P-value}
#'  \item{padj}{Adjusted p-value}
#'  \item{avgExpr}{Mean expression in the test group indicated by the "group"
#'    field. Only available for wilcoxon tests.}
#'  \item{statistic}{Wilcoxon rank-sum test statistic. Only available for
#'    wilcoxon tests.}
#'  \item{auc}{Area under the ROC curve. Only available for wilcoxon tests.}
#'  \item{pct_in}{Percentage of cells in the test group, indicated by the
#'    "group" field, that express the feature. Only available for wilcoxon
#'    tests.}
#'  \item{pct_out}{Percentage of cells in the control group or other cells, as
#'    explained for the "group" field, that express the feature. Only available
#'    for wilcoxon tests.}
#' @rdname liger-DEG
#' @export
#' @examples
#' \donttest{
#' pbmc$leiden_cluster <- pbmcPlot$leiden_cluster
#'
#' # Identify cluster markers
#' degStats1 <- runMarkerDEG(pbmc, conditionBy = "leiden_cluster")
#'
#' # Compare "stim" data against "ctrl" data within each cluster
#' degStats3 <- runPairwiseDEG(pbmc, groupTest = "stim", groupCtrl = "ctrl",
#'                             variable1 = "dataset",
#'                             splitBy = "leiden_cluster",
#'                             minCellPerRep = 4)
#' }
runPairwiseDEG <- function(
        object,
        groupTest,
        groupCtrl,
        variable1 = NULL,
        variable2 = NULL,
        splitBy = NULL,
        method = c("pseudoBulk", "wilcoxon"),
        usePeak = FALSE,
        useReplicate = "dataset",
        nPsdRep = NULL,
        minCellPerRep = 3,
        printDiagnostic = FALSE,
        chunk = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose", TRUE)
) {
    method <- match.arg(method)
    if (!is.null(splitBy)) {
        splitVar <- .fetchCellMetaVar(object, splitBy,
                                      checkCategorical = TRUE, drop = TRUE,
                                      droplevels = TRUE)
        splitVar <- interaction(splitVar, drop = TRUE)
        splitGroup <- lapply(levels(splitVar), function(x) {
            which(splitVar == x)
        })
        names(splitGroup) <- levels(splitVar)
    } else {
        splitGroup <- list(seq(ncol(object)))
    }

    if (is.null(variable1) && is.null(variable2)) {
        # Directly using cell index
        groups <- list(
            .idxCheck(object, groupTest, "cell"),
            .idxCheck(object, groupCtrl, "cell")
        )
        group1Name <- "test"
        group2Name <- "control"
        names(groups) <- c("test", "control")
    } else if (!is.null(variable1)) {
        var1 <- .fetchCellMetaVar(object, variable1,
                                  checkCategorical = TRUE, drop = TRUE,
                                  droplevels = TRUE)
        group1Idx <- which(var1 %in% groupTest)
        group1Name <- paste(groupTest, collapse = ".")
        if (is.null(variable2)) {
            variable2 <- variable1
            var2 <- var1
        } else {
            var2 <- .fetchCellMetaVar(object, variable2,
                                      checkCategorical = TRUE, drop = TRUE,
                                      droplevels = TRUE)
        }
        group2Idx <- which(var2 %in% groupCtrl)
        group2Name <- paste(groupCtrl, collapse = ".")
        groups <- list(group1Idx, group2Idx)
        names(groups) <- c(group1Name, group2Name)
    } else {
        cli::cli_abort("Please see {.code ?runPairwiseDEG} for usage.")
    }
    resultList <- list()
    groupsMeta <- list()
    for (i in seq_along(splitGroup)) {
        splitName <- names(splitGroup)[i]
        splitIdx <- splitGroup[[i]]

        if (isTRUE(verbose)) {
            if (length(splitGroup) > 1) {
                cli::cli_alert_info(
                    "Running DEG within: {.val {splitName}}"
                )
            }
        }

        groups.sub <- lapply(groups, function(x) {
            intersect(x, splitIdx)
        })
        names(groups.sub) <- sapply(names(groups), function(x) {
            paste(c(splitName, x), collapse = ".")
        })
        if (length(groups.sub[[1]]) == 0) {
            cli::cli_alert_warning("No cell selected for group {.val {names(groups.sub)[1]}} when split by {.val {splitBy}} in level {.val {splitName}}. Skipping.")
            next
        }
        if (length(groups.sub[[2]]) == 0) {
            cli::cli_alert_warning("No cell selected for group {.val {names(groups.sub)[1]}} when split by {.val {splitBy}} in level {.val {splitName}}. Skipping.")
            next
        }
        groupsMeta[[names(groups.sub)[1]]] <- groups.sub
        result <- .runDEG(object, groups = groups.sub, method = method,
                          usePeak = usePeak, useReplicate = useReplicate,
                          nPsdRep = nPsdRep,
                          minCellPerRep = minCellPerRep, seed = seed,
                          printDiagnostic = printDiagnostic,
                          skipTwoGroup = TRUE, chunk = chunk,
                          verbose = verbose)
        if (nrow(result) > 0) {
            result <- result[result$group == names(groups.sub)[1],]
        }

        # attributes(result)$meta <- list(
        #     groupTest = groupTest,
        #     variable1 = variable1,
        #     groupCtrl = groupCtrl,
        #     variable2 = variable2
        # )
        resultList[[i]] <- result
    }
    resultList <- Reduce(rbind, resultList)
    attr(resultList, "groups") <- groupsMeta
    return(resultList)
}

#' @rdname liger-DEG
#' @export
#' @param conditionBy \code{cellMeta} variable(s). Marker detection will be
#' performed for each level of this variable. Multiple variables will be
#' combined. Default \code{NULL} uses default cluster.
#' @param useDatasets Datasets to perform marker detection within. Default
#' \code{NULL} will use all datasets.
runMarkerDEG <- function(
        object,
        conditionBy = NULL,
        splitBy = NULL, # The previous by dataset strategy
        method = c("pseudoBulk", "wilcoxon"),
        useDatasets = NULL,
        usePeak = FALSE,
        useReplicate = "dataset",
        nPsdRep = NULL,
        minCellPerRep = 3,
        printDiagnostic = FALSE,
        chunk = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose", TRUE)
) {
    useDatasets <- .checkUseDatasets(object, useDatasets)
    allCellIdx <- seq(ncol(object))[object$dataset %in% useDatasets]
    conditionBy <- conditionBy %||% object@uns$defaultCluster
    if (is.null(conditionBy)) {
        cli::cli_abort("No {.var conditionBy} given or default cluster not set.")
    }
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
    groups <- split(allCellIdx, conditionBy)
    if (nlevels(splitBy) <= 1) {
        result <- .runDEG(object, groups = groups, method = method,
                          usePeak = usePeak, useReplicate = useReplicate,
                          nPsdRep = nPsdRep,
                          minCellPerRep = minCellPerRep,
                          printDiagnostic = printDiagnostic,
                          skipTwoGroup = FALSE, chunk = chunk,
                          seed = seed, verbose = verbose)
        attr(result, "groups") <- groups
    } else {
        resultList <- list()
        groupsMeta <- list()
        for (i in seq_along(levels(splitBy))) {
            if (isTRUE(verbose)) {
                cli::cli_alert_info(
                    "Running DEG within: {.val {levels(splitBy)[i]}}"
                )
            }
            subIdx <- which(splitBy == levels(splitBy)[i])
            subGroups <- lapply(groups, intersect, subIdx)
            names(subGroups) <- paste0(levels(splitBy)[i], '.', names(groups))
            groupsMeta <- c(groupsMeta, subGroups)
            resultList[[levels(splitBy)[i]]] <- .runDEG(
                object, groups = subGroups, method = method,
                usePeak = usePeak, useReplicate = useReplicate,
                nPsdRep = nPsdRep,
                minCellPerRep = minCellPerRep,
                printDiagnostic = printDiagnostic, skipTwoGroup = FALSE,
                chunk = chunk, seed = seed, verbose = verbose
            )
        }
        result <- Reduce(rbind, resultList)
        attr(result, "groups") <- groupsMeta
    }

    return(result)
}

#' @rdname liger-DEG
#' @export
#' @param data.use Same as \code{useDatasets}.
#' @param compare.method Choose from \code{"clusters"} (default) or
#' \code{"datasets"}. \code{"clusters"} compares each cluster against all other
#' cells, while \code{"datasets"} run within each cluster and compare each
#' dataset against all other datasets.
runWilcoxon <- function(
        object,
        data.use = NULL,
        compare.method = c("clusters", "datasets")
) {
    lifecycle::deprecate_warn(
        "1.99.0", "runWilcoxon()",
        details = "Please use `runMarkerDEG()` with `method = 'wilcoxon'` instead."
    )
    compare.method <- match.arg(compare.method)
    if (compare.method == "clusters") {
        res <- runMarkerDEG(object, conditionBy = object@uns$defaultCluster,
                            splitBy = NULL, method = "wilcoxon")
    } else if (compare.method == "datasets") {
        res <- runMarkerDEG(object, conditionBy = "dataset",
                            splitBy = object@uns$defaultCluster,
                            method = "wilcoxon")
    }
    return(res)
}

# groups - As underlying function, this must be organized into list of numeric
# cell index vectors.
.runDEG <- function(
        object,
        groups,
        method = c("pseudoBulk", "wilcoxon"),
        usePeak = FALSE,
        useReplicate = "dataset",
        nPsdRep = NULL,
        minCellPerRep = 10,
        printDiagnostic = FALSE,
        skipTwoGroup = TRUE,
        chunk = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose", TRUE)
) {
    method <- match.arg(method)
    allCellIdx <- Reduce(c, groups)
    if (length(allCellIdx) == 0)
        cli::cli_abort(c(x = "No cell selected"))
    allCellBC <- colnames(object)[allCellIdx]
    datasetInvolve <- levels(object$dataset[allCellIdx, drop = TRUE])
    var <- factor(rep(names(groups), lengths(groups)), levels = names(groups))
    if (isTRUE(usePeak)) {
        useDatasets <- .checkUseDatasets(object, useDatasets = datasetInvolve,
                                         modal = "atac")
    } else {
        useDatasets <- .checkUseDatasets(object, useDatasets = datasetInvolve)
    }
    slot <- .DE.checkDataAvail(object, datasetInvolve, method, usePeak)
    dataList <- getMatrix(object, slot, datasetInvolve, returnList = TRUE)
    # mat <- mergeSparseAll(dataList, mode = "intersection")
    features <- Reduce(intersect, lapply(dataList, rownames))
    if (length(features) == 0) {
        cli::cli_abort("No shared feature available from datasets involved for the test.
                       Datasets involved: {.val {datasetInvolve}}")
    }
    featureOrder <- stats::setNames(seq_along(features), features)
    # dataList <- lapply(dataList, function(x) x[features, , drop = FALSE])

    # mat <- Reduce(cbind, dataList)

    # mat <- mat[, allCellBC, drop = FALSE]
    if (method == "wilcoxon") {
        chunk <- chunk %||% length(features)
        nchunk <- ceiling(length(features)/chunk)
        resultList <- list()
        cli::cli_progress_bar("Wilcoxon rank-sum test", total = nchunk, clear = FALSE)
        for (i in seq_len(nchunk)) {
            start <- (i - 1) * chunk + 1
            end <- min(i * chunk, length(features))
            mat <- extractMergedNormData(
                object,
                cellIdx = allCellBC,
                featureIdx = features[start:end]
            )
            mat <- log1p(1e10*mat)
            resultList[[i]] <- wilcoxauc(mat, var)
            if (nchunk > 1) gc()
            cli::cli_progress_update(set = i)
        }
        result <- Reduce(rbind, resultList)
        result <- result[order(result$group, featureOrder[result$feature]), ]
        if (nchunk > 1) {
            result %<>%
                dplyr::group_by(.data[['group']]) %>%
                dplyr::mutate(padj = stats::p.adjust(.data[['pval']], method = "BH")) %>%
                as.data.frame()
        }

        # result$group <- factor(result$group, levels = levels(var))
    } else if (method == "pseudoBulk") {
        resultList <- list()
        if (isTRUE(verbose)) {
            if (nlevels(var) <= 2) {
                cliID <- cli::cli_process_start("Calling pairwise DESeq2 Wald test")
            } else if (nlevels(var) > 2) {
                if (isFALSE(printDiagnostic)) {
                    cliID_pb <- cli::cli_progress_bar(name = "DESeq2 Wald test",
                                                      total = nlevels(var))
                }
            }
        }

        for (i in seq_along(levels(var))) {
            testName <- levels(var)[i]
            if (isTRUE(verbose) &&
                isTRUE(printDiagnostic) &&
                nlevels(var) > 2) {
                cliID_perLevel <- cli::cli_process_start("Working on {.val {testName}}")
            }
            result <- tryCatch({
                subVar <- factor(ifelse(var == testName, testName, "others"),
                                 levels = c(testName, "others"))
                aux <- calcPctInOut(
                    object,
                    cellIdx = allCellIdx,
                    features = features,
                    groups = subVar
                )
                # `useReplicate` can be a vector of multiple variables
                replicateAnn <- .fetchCellMetaVar(
                    object, useReplicate,
                    cellIdx = allCellIdx,
                    drop = FALSE,
                    checkCategorical = TRUE,
                    droplevels = TRUE
                )
                # Collapse multiple variables into one factor
                replicateAnn <- interaction(replicateAnn, drop = TRUE)

                replicateAnn <- setupRepVar(
                    groups = subVar,
                    existing = replicateAnn,
                    nPsdRep = nPsdRep,
                    seed = seed,
                    verbose = printDiagnostic
                )
                pbs <- makePseudoBulk(
                    object,
                    cellIdx = allCellIdx,
                    features = features,
                    groups = subVar,
                    replicateAnn = replicateAnn,
                    minCellPerRep = minCellPerRep,
                    verbose = printDiagnostic
                )
                pb <- pbs[[1]]
                subVar <- pbs[[2]]


                # resultList[[testName]] <- .callDESeq2(pb, subVar, printDiagnostic)
                result <- .callDESeq2(pb, subVar, printDiagnostic)
                if (length(levels(var)) <= 2) {
                    if (isTRUE(verbose)) cli::cli_process_done(id = cliID)
                } else {
                    if (isTRUE(verbose)) {
                        if (isTRUE(printDiagnostic)) {
                            cli::cli_process_done(id = cliID_perLevel)
                        } else {
                            cli::cli_progress_update(set = i)
                        }
                    }
                }
                cbind(result, aux)
            }, error = function(e) {
                cli::cli_alert_danger(
                    "Error when computing on {.val {testName}}: {e$message}"
                )
                cli::cli_alert_warning(
                    "Empty result returned for this test."
                )
                if (length(levels(var)) <= 2) {
                    if (isTRUE(verbose)) cli::cli_process_failed(id = cliID)
                } else {
                    if (isTRUE(verbose)) {
                        if (isTRUE(printDiagnostic)) {
                            cli::cli_process_failed(id = cliID_perLevel)
                        } else {
                            cli::cli_progress_update(set = i, id = cliID_pb)
                        }
                    }
                }
                return(data.frame(
                    feature = character(0),
                    group = character(0),
                    logFC = numeric(0),
                    pval = numeric(0),
                    padj = numeric(0)
                ))
            })
            resultList[[testName]] <- result
            if (length(levels(var)) <= 2 && isTRUE(skipTwoGroup)) break
        }
        result <- Reduce(rbind, resultList)

    }
    rownames(result) <- NULL
    result$group <- factor(result$group, levels = levels(var)[levels(var) %in% result$group])
    return(result)
}

.DE.checkDataAvail <- function(object, useDatasets, method, usePeak) {
    if (isH5Liger(object, useDatasets)) { # nocov start
        cli::cli_abort(
            c("HDF5 based datasets detected but is not supported. ",
              "i" = "Try {.code object.sub <- downsample(object, useSlot = 'normData')} to create another object with in memory data")
        )
    } # nocov end
    if (method == "wilcoxon") {
        slot <- ifelse(usePeak, "normPeak", "normData")
    } else if (method == "pseudoBulk") {
        if (!requireNamespace("DESeq2", quietly = TRUE)) # nocov start
            cli::cli_abort(
                "Package {.pkg DESeq2} is needed for this function to work.
                Please install it by command:
                {.code BiocManager::install('DESeq2')}"
            ) # nocov end
        slot <- ifelse(usePeak, "rawPeak", "rawData")
    }
    allAvail <- all(sapply(useDatasets, function(d) {
        ld <- dataset(object, d)
        !is.null(methods::slot(ld, slot))
    }))
    if (!allAvail)
        cli::cli_abort(
            c("{.field {slot}} not all available for involved datasets: {.val {useDatasets}}",
              "i" = "{.code method = '{method}'}; {.code usePeak = {usePeak}}")
        )
    return(slot)
}

###################### Pseudo-bulk Method helper ###############################

# This function checks for real replicate availability for each comparison group
# and settings for pseudo-replicates in situations where real replicates are not
# available or too few.
# groups - a factor with two levels defining the comparison groups
# existing - a factor defining real replicates, or NULL if real replicate is not
#            available
# nRep - number of pseudo-replicates to generate. If NULL, check for number of
# existing replicates and prompt warning if too few. Pseudo-replicates will be
# generated for each group and each replicate if `existing` and `nRep` are not
# NULL.
# seed - random seed to use for generating pseudo-replicates
# verbose - whether to show not-too-important diagnostic messages
# Returns a factor with the same length as `groups`, guiding how pseudo-bulk
# should be finally aggregated. Value string is formatted by
# 1. `group.realRep`, when real replicates are available and enough
# 2. `group.psdRep`, when real replicates are not available and peudo-replicates
#    are generated accordingly
# 3. `group.realRep.psdRep`, when real replicates are available but not enough
#    and pseudo-replicates are generated within per group per real-rep
#    accordingly
setupRepVar <- function(
        groups,
        existing,
        nPsdRep = NULL,
        seed = 1,
        verbose = TRUE
) {
    set.seed(seed)
    # Initialize final returned value
    allRepVar <- rep(NA, length(groups))

    for (i in seq_along(levels(groups))) {
        # Need to check whether there are enough real replicates for each group first
        # Otherwise check whether to generate pseudo-replicates
        groupIdx <- groups == levels(groups)[i]
        if (length(existing) == 0) {
            # No real replicates given, must do pseudo-replicates
            nPsdRepGroup <- nPsdRep %||% 3
            repVarSub <- sample(seq_len(sum(groupIdx)) %% nPsdRepGroup) + 1
            repVarSub <- paste0(levels(groups)[i], ".rep", repVarSub)
        } else {
            existingSub <- droplevels(existing[groups == levels(groups)[i]])
            if (nlevels(existingSub) < 3) {
                # Too few real replicates
                if (is.null(nPsdRep)) {
                    if (isTRUE(verbose)) {
                        cli::cli_alert_warning(
                            "Too few real replicates for condition {.val {levels(groups)[i]}}. Creating 3 pseudo-replicates among each real replicate."
                        )
                    }
                    nPsdRepGroup <- 3
                } else {
                    nPsdRepGroup <- nPsdRep
                }
            } else {
                nPsdRepGroup <- nPsdRep %||% 1
            }
            if (nPsdRepGroup > 1) {
                repVarSub <- rep(NA, length(existingSub))
                for (j in seq_along(levels(existingSub))) {
                    existingSubIdx <- existingSub == levels(existingSub)[j]
                    repVarSubSub <- sample(seq_len(sum(existingSubIdx)) %% nPsdRepGroup) + 1
                    repVarSubSub <- sprintf(
                        '%s.%s.rep%s',
                        levels(groups)[i], levels(existingSub)[j], repVarSubSub
                    )
                    repVarSub[existingSubIdx] <- repVarSubSub
                }
            } else {
                # Don't psdrep anyway
                repVarSub <- paste0(levels(groups)[i], ".", as.character(existingSub))
            }
        }
        allRepVar[groupIdx] <- repVarSub
    }
    return(factor(allRepVar))
}

# This function aggregate raw counts data directly from un-merged matrices in
# a liger object, guided by `replicateAnn`. Returns aggregated pseudo-bulk
# matrix together with the comparison group assignment for each pseudo-bulk.
# The key idea of making this quite-complexed function is to avoid merging
# involved matrix into a giant one, but instead only take what we need and
# reduce (add) it into the final pseudo-bulk matrix.
# object - liger object
# cellIdx - cell index being involved for the test
# features - character vector for features to be included, should already be
#            generated in upstream steps.
# groups - factor defining the comparison groups for each cell specified by
#          `cellIdx`
# replicateAnn - factor defining the replicate assignment for each cell
#                specified by `cellIdx`.
# minCellPerRep - minimum number of cells required for each replicate. Replicates
#                 with less cells will be ignored.
# verbose - whether to show not-too-important diagnostic messages
# Returns a list of two elements, the first being a matrix of length(features)
# rows and `nlevels(replicateAnn)` columns. Could be less for ncols due to
# `minCellPerRep` filtering. The second element is the comparison group
# belonging factor for each pseudo-bulk.
makePseudoBulk <- function(
        object,
        cellIdx,
        features,
        groups,
        replicateAnn,
        minCellPerRep,
        verbose
) {
    # First check if each replicate is having a healthy (large enough) size
    repTab <- table(replicateAnn)
    if (isTRUE(verbose)) {
        cli::cli_alert_info("Replicate sizes:")
        print(repTab)
    }

    ignored <- names(repTab)[repTab < minCellPerRep]
    keep <- names(repTab)[repTab >= minCellPerRep]
    if (length(ignored) > 0) {
        cli::cli_alert_warning(
            "Ignoring replicates (size in bracket) with too few cells: {.val {paste0(ignored, ' (', repTab[ignored], ')')}}"
        )
        cli::cli_alert_info(
            "Consider decrease {.field minCellPerRep} to exclude less replicates or/and lower {.field nPsdRep} to generate larger pseudo-replicates."
        )
    }


    # Then check if after removal, each condition group still has enough
    # replicates
    for (gr in levels(groups)) {
        repPerGroup <- levels(droplevels(replicateAnn[groups == gr]))
        repPerGroup <- repPerGroup[repPerGroup %in% keep]
        if (length(repPerGroup) < 3) {
            cli::cli_abort(
                "Too few replicates with more than {.val {minCellPerRep}} cells ({.var minCellPerRep}) for condition {.val {gr}}."
            )
        }
    }
    cellIdx <- cellIdx[replicateAnn %in% keep]
    groups <- groups[replicateAnn %in% keep]
    replicateAnn <- droplevels(replicateAnn[replicateAnn %in% keep])

    datasetInvolved <- as.character(unique(object$dataset[cellIdx]))
    # Initialize the output matrix
    pseudoBulk <- matrix(
        data = 0,
        nrow = length(features),
        ncol = nlevels(replicateAnn),
        dimnames = list(features, levels(replicateAnn))
    )
    # This won't work if `groups` is not super-categorical of `replicateAnn`
    # But this should not happen if no bug in upstream
    groupOut <- data.frame(rep = replicateAnn, groups = groups) %>%
        dplyr::group_by(.data[['rep']], .data[['groups']]) %>%
        dplyr::count() %>%
        `[[`('groups')

    # Go through each dataset, find out the cells involved and aggregate
    # them into the initialized pseudo-bulk following the guidance of
    # `replicateAnn`

    # Meanwhile, calculate auxiliary metrics `pct_in` and `pct_out` during the
    # for loop for each dataset.

    # repAnnExpand - broad cast the repAnn for only involved cells to all cells
    repAnnExpand <- rep(NA, ncol(object))
    repAnnExpand[cellIdx] <- replicateAnn
    repAnnExpand <- factor(repAnnExpand)
    for (i in seq_along(datasetInvolved)) {
        # dn - dataset name
        dn <- datasetInvolved[i]
        raw <- rawData(object, dn)
        repAnnExpandSub <- repAnnExpand[object$dataset == dn]
        updatePseudoBulkRcpp(
            psdBulk = pseudoBulk,
            sparseRaw = raw,
            featureIdx = match(rownames(raw), features) - 1,
            repIdx = as.integer(repAnnExpandSub) - 1
        )
    }
    return(list(pseudoBulk, groupOut))
}

calcPctInOut <- function(
    object,
    cellIdx,
    features,
    groups
) {
    datasetInvolved <- as.character(unique(object$dataset[cellIdx]))
    # Initialize the output matrix
    # `nCellExpr` is a matrix of n features rows and 2 cols. The first column
    # stores number of cells in test group that express each feature, and the
    # second column stores number of cells in control group that express each
    # feature.
    # After counting, we will calculate `pct_in` and `pct_out` for each feature
    # by dividing the first column by total number of cells in test group and
    # the second column by total number of cells in control group.
    nCellExpr <- matrix(
        data = 0,
        nrow = length(features),
        ncol = 2,
        dimnames = list(features, c("pct_in", "pct_out"))
    )
    groupExpand <- rep(NA, ncol(object))
    # When inserting a factor (`group`) in to the NA vector, the factor
    # is automatically converted to 1-based integers.
    groupExpand[cellIdx] <- groups
    groupExpand <- groupExpand - 1

    # Go through each dataset, find out the cells involved and calculate
    # `pct_in` and `pct_out` for each feature
    for (i in seq_along(datasetInvolved)) {
        # dn - dataset name
        dn <- datasetInvolved[i]
        raw <- rawData(object, dn)
        # The `updateNCellExprRcpp` Rcpp function in-place updates `nCellExpr`
        # using each raw matrix so no redundant memory allocation is needed.
        # Information needed:
        # - index that tells which cells in `raw` are for the test group and
        #   which are for the control group
        updateNCellExprRcpp(
            out = nCellExpr,
            sparseRaw = raw,
            featureIdx = match(rownames(raw), features) - 1,
            groupVar = groupExpand[object$dataset == dn]
        )
    }
    # Finally, calculate the percentage
    nCellExpr[,1] <- nCellExpr[,1]/sum(groups == levels(groups)[1]) * 100
    nCellExpr[,2] <- nCellExpr[,2]/sum(groups == levels(groups)[2]) * 100
    return(nCellExpr)
}

# makePseudoBulkOld <- function(mat, replicateAnn, minCellPerRep, verbose = TRUE) {
#     # mat - Extracted and contatenated matrix. intersection of genes by
#     #       c(groupTest, groupCtrl) cells
#     # groups - list of groups
#     # replicateAnn - data.frame of replicate annotation, with rownames as
#     #                barcodes and columns as variables
#
#     # Check whether enough replicates per condition
#     for (gr in levels(replicateAnn$groups)) {
#         subrep <- replicateAnn[replicateAnn$groups == gr,]
#         splitLabel <- interaction(subrep, drop = TRUE)
#         if (nlevels(splitLabel) < 2) {
#             cli::cli_abort(
#                 c("Too few replicates for condition {.val {gr}}. Cannot create pseudo-bulks.",
#                   "i" = "Please consider creating pseudo-replicates or using {.code method = 'wilcoxon'} instead.")
#             )
#         }
#     }
#     splitLabel <- interaction(replicateAnn, drop = TRUE)
#     repSizeTab <- table(splitLabel)
#     if (verbose) {
#         cli::cli_alert_info("Replicate sizes:")
#         print(repSizeTab)
#     }
#     labelCounts <- table(splitLabel)
#     ignored <- names(labelCounts)[labelCounts < minCellPerRep]
#     if (length(ignored) > 0) {
#         cli::cli_alert_warning(
#             "Ignoring replicates (size in bracket) with too few cells: {.val {paste0(ignored, ' (', repSizeTab[ignored], ')')}}"
#         )
#         cli::cli_alert_info(
#             "Consider decrease {.field minCellPerRep} to exclude less replicates or/and lower {.field nPsdRep} to generate larger pseudo-replicates."
#         )
#     }
#     keep <- names(labelCounts)[labelCounts >= minCellPerRep]
#     idx <- splitLabel %in% keep
#     splitLabel <- splitLabel[idx, drop = TRUE]
#     mat <- mat[, idx, drop = FALSE]
#     replicateAnn <- replicateAnn[idx, , drop = FALSE]
#
#     pseudoBulks <- colAggregateSums_sparse(mat, as.integer(splitLabel) - 1,
#                                            nlevels(splitLabel))
#     dimnames(pseudoBulks) <- list(rownames(mat), levels(splitLabel))
#     pseudoBulks <- pseudoBulks[rowSums(pseudoBulks) > 0,]
#     return(list(pseudoBulks, replicateAnn))
# }

.callDESeq2 <- function(pseudoBulks, groups,
                         verbose = getOption("ligerVerbose", TRUE)) {
    # DESeq2 workflow

    ## NOTE: DESeq2 wishes that the contrast/control group is the first level
    ## whereas we required it as the second in upstream input. So we need to
    ## reverse it here.
    groups <- stats::relevel(groups, ref = levels(groups)[2])
    groupsDropped <- droplevels(groups)
    if (nlevels(groupsDropped) < 2) {
        cli::cli_abort("No enough replicates for conditions being compared.")
    }
    ## Now levels(groups)[1] becomes control and levels(groups)[2] becomes
    ## the test group
    suppressMessages({
        des <- DESeq2::DESeqDataSetFromMatrix(
            countData = pseudoBulks,
            colData = data.frame(groups = groups),
            design = stats::formula("~groups")
        )
    })
    des <- DESeq2::DESeq(des, test = "Wald", quiet = !verbose)
    res <- DESeq2::results(des, contrast = c("groups", levels(groups)[2],
                                             levels(groups)[1]))
    res <- .DataFrame.as.data.frame(res)
    res$feature <- rownames(res)
    rownames(res) <- NULL
    res$group <- levels(groups)[2]
    res <- res[, c(7, 8, 2, 5, 6)]
    colnames(res) <- c("feature", "group", "logFC", "pval", "padj")

    return(res)
}


####################### Wilcoxon rank-sum test helper ##########################

extractMergedNormData <- function(
        object,
        cellIdx = NULL,
        featureIdx = NULL
) {
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    datasetInvolved <- unique(object$dataset[cellIdx])
    cellID <- colnames(object)[cellIdx]
    if (is.null(featureIdx)) {
        # Use intersection by default
        featureIdx <- Reduce(
            intersect,
            lapply(
                normData(object, datasetInvolved),
                rownames
            )
        )
    }
    out <- NULL
    for (i in seq_along(datasetInvolved)) {
        dn <- datasetInvolved[i]
        ldFeatureIdx <- .idxCheck(dataset(object, dn), featureIdx, "feature")
        ldCellIdx <- match(cellID, colnames(dataset(object, dn)))
        ldCellIdx <- ldCellIdx[!is.na(ldCellIdx)]
        out <- cbind(out, normData(object, dn)[ldFeatureIdx, ldCellIdx, drop = FALSE])
    }
    out[, cellID]
}

# X: matrix of data to be tested
# y: grouping label of columns of X
# Rcpp source code located in src/wilcoxon.cpp
wilcoxauc <- function(x, clusterVar) {
    if (methods::is(x, 'dgTMatrix')) x <- methods::as(x, 'CsparseMatrix') # nocov start
    if (methods::is(x, 'TsparseMatrix')) x <- methods::as(x, 'CsparseMatrix')
    if (is.null(row.names(x))) {
        rownames(x) <- paste0('Feature', seq(nrow(x)))
    } # nocov end
    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    clusterVar <- droplevels(clusterVar)
    groupSize <- as.numeric(table(clusterVar))

    ## Compute primary statistics
    n1n2 <- groupSize * (ncol(x) - groupSize)
    # rankRes - list(X_ranked, ties), where X_ranked is obs x feature
    xRanked <- Matrix::t(x)
    # This computes the ranking of non-zero values and the ties
    ties <- cpp_rank_matrix_dgc(xRanked@x, xRanked@p,
                                nrow(xRanked), ncol(xRanked))
    # ranksRes <- list(X_ranked = xT, ties = ties)

    # rankRes <- colRanking(x)
    ustat <- computeUstat(xRanked, clusterVar, n1n2, groupSize)
    auc <- t(ustat / n1n2)
    pvals <- computePval(ustat, ties, ncol(x), n1n2)
    fdr <- apply(pvals, 2, function(p) stats::p.adjust(p, 'BH'))

    ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
    groupSums <- colAggregateSum_sparse(x, as.integer(clusterVar) - 1, length(unique(clusterVar)))
    # groupSums <- colAggregateSum(x, clusterVar)
    group_nnz <- colNNZAggr_sparse(x, as.integer(clusterVar) - 1, length(unique(clusterVar)))
    # group_nnz <- colNNZAggr(x, clusterVar)
    group_pct <- t(sweep(group_nnz, 1, as.numeric(table(clusterVar)), "/"))

    group_pct_out <- sweep(-group_nnz, 2, colSums(group_nnz), "+")
    group_pct_out <- sweep(group_pct_out, 1,
                           as.numeric(length(clusterVar) - table(clusterVar)),
                           "/")
    group_pct_out <- t(group_pct_out)

    groupMeans <- t(sweep(groupSums, 1, as.numeric(table(clusterVar)), "/"))

    cs <- colSums(groupSums)
    gs <- as.numeric(table(clusterVar))
    lfc <- Reduce(cbind, lapply(seq_along(levels(clusterVar)), function(g) {
        groupMeans[, g] - (cs - groupSums[g, ])/(length(clusterVar) - gs[g])
    }))

    data.frame(
        feature = rep(row.names(x), times = length(levels(clusterVar))),
        group = factor(rep(levels(clusterVar), each = nrow(x)),
                       levels = levels(clusterVar)),
        avgExpr = as.numeric(groupMeans),
        logFC = as.numeric(lfc),
        statistic = as.numeric(t(ustat)),
        auc = as.numeric(auc),
        pval = as.numeric(pvals),
        padj = as.numeric(fdr),
        pct_in = as.numeric(100 * group_pct),
        pct_out = as.numeric(100 * group_pct_out)
    )
}

computeUstat <- function(Xr, cols, n1n2, groupSize) {
    grs <- rowAggregateSum_sparse(Xr, as.integer(cols) - 1, length(unique(cols)))
    # grs <- rowAggregateSum(Xr, cols)

    # if (inherits(Xr, 'dgCMatrix')) {
    # With the ranking of only non-zero features, here the tie-ranking of
    # zeros need to be added.
    nnz <- rowNNZAggr_sparse(Xr, as.integer(cols) - 1, length(unique(cols)))
    gnz <- groupSize - nnz
    zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
    ustat <- t((t(gnz) * zero.ranks)) + grs - groupSize*(groupSize + 1)/2
    # } else {
    #     ustat <- grs - groupSize * (groupSize + 1) / 2
    # }
    return(ustat)
}

computePval <- function(ustat, ties, N, n1n2) {
    z <- ustat - .5 * n1n2
    z <- z - sign(z) * .5
    .x1 <- N ^ 3 - N
    .x2 <- 1 / (12 * (N ^ 2 - N))
    rhs <- unlist(lapply(ties, function(tvals) {
        (.x1 - sum(tvals ^ 3 - tvals)) * .x2
    }))
    usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
    z <- t(z / usigma)
    pvals <- matrix(2 * stats::pnorm(-abs(as.numeric(z))), ncol = ncol(z))
    return(pvals)
}





######################## Visualization #########################################

#' Create heatmap for showing top marker expression in conditions
#' @export
#' @param object A \linkS4class{liger} object, with normalized data and metadata
#' to annotate available.
#' @param result The data.frame returned by \code{\link{runMarkerDEG}}.
#' @param topN Number of top features to be plot for each group. Default
#' \code{5}.
#' @param lfcThresh Hard threshold on logFC value. Default \code{1}.
#' @param padjThresh Hard threshold on adjusted P-value. Default \code{0.05}.
#' @param pctInThresh,pctOutThresh Threshold on expression percentage. These
#' mean that a feature will only pass the filter if it is expressed in more than
#' \code{pctInThresh} percent of cells in the corresponding cluster. Similarly
#' for \code{pctOutThresh}. Only applied when these metrics are available.
#' Default \code{50} percent for both.
#' @param dedupBy When ranking by padj and logFC and a feature is ranked as top
#' for multiple clusters, assign this feature as the marker of a cluster when
#' it has the largest \code{"logFC"} in the cluster or has the lowest
#' \code{"padj"}. Default \code{"logFC"}.
#' @param groupBy Cell metadata variable names for cell grouping. Downsample
#' balancing will also be aware of this. Default \code{"dataset"} and the
#' default cluster.
#' @param groupSize Maximum number of cells in each group to be downsampled for
#' plotting. Default \code{50}.
#' @param column_title Title on the column. Default \code{NULL}.
#' @inheritDotParams plotGeneHeatmap cellAnnotation
#' @inheritDotParams .plotHeatmap transpose showCellLabel showCellLegend showFeatureLabel cellAnnColList featureAnnColList scale trim baseSize cellTextSize featureTextSize cellTitleSize featureTitleSize legendTextSize legendTitleSize viridisOption viridisDirection RColorBrewerOption
#' @return A \link[ComplexHeatmap]{HeatmapList-class} object.
#' @examples
#' defaultCluster(pbmc) <- pbmcPlot$leiden_cluster
#' pbmc <- normalize(pbmc)
#' plotMarkerHeatmap(pbmc, deg.marker)
plotMarkerHeatmap <- function(
        object,
        result,
        topN = 5,
        lfcThresh = 1,
        padjThresh = 0.05,
        pctInThresh = 50,
        pctOutThresh = 50,
        dedupBy = c("logFC", "padj"),
        groupBy = NULL,
        groupSize = 50,
        column_title = NULL,
        ...
) {
    dedupBy <- match.arg(dedupBy)
    if (dedupBy == "logFC") {
        result <- result[order(result[[dedupBy]], decreasing = TRUE), ]
    } else if (dedupBy == "padj") {
        result <- result[order(result[[dedupBy]], decreasing = FALSE), ]
    }
    pctInThresh <- pctInThresh %||% 0
    pctOutThresh <- pctOutThresh %||% 100
    groupBy <- groupBy %||% c("dataset", object@uns$defaultCluster)
    result <- result[!duplicated(result$feature), ]
    result <- result %>%
        dplyr::filter(
            .data$logFC > lfcThresh,
            .data$padj < padjThresh
        )
    if ("pct_in" %in% colnames(result) &&
        "pct_out" %in% colnames(result)) {
        result <- result %>%
            dplyr::filter(
                .data$pct_in > pctInThresh,
                .data$pct_out < pctOutThresh
            )
    }
    result <- result %>%
        dplyr::group_by(.data[["group"]]) %>%
        dplyr::arrange(
            .data[["padj"]],
            -.data[["logFC"]],
            .by_group = TRUE
        ) %>%
        dplyr::slice_head(n = topN) %>%
        as.data.frame()
    cellIdx <- downsample(object, maxCells = groupSize, balance = groupBy,
                          returnIndex = TRUE)
    featureAnn <- result[, "group", drop = FALSE]

    rownames(featureAnn) <- result$feature
    colnames(featureAnn) <- "marker"
    plotGeneHeatmap(object, features = result$feature,
                    cellIdx = cellIdx,
                    useCellMeta = groupBy,
                    featureAnnotation = featureAnn,
                    cellSplitBy = rev(groupBy),
                    featureSplitBy = "marker",
                    showFeatureLegend = FALSE,
                    cluster_columns = FALSE,
                    cluster_column_slices = FALSE,
                    cluster_rows = FALSE,
                    cluster_row_slices = FALSE,
                    column_title = column_title,
                    ...)
}


#' Create heatmap for pairwise DEG analysis result
#' @export
#' @param object A \linkS4class{liger} object, with normalized data and metadata
#' to annotate available.
#' @param result The data.frame returned by \code{\link{runPairwiseDEG}}.
#' @param group The test group name among the result to be shown. Must specify
#' only one if multiple tests are available (i.e. split test). Default
#' \code{NULL} works with single-test result and raises error with split-test
#' result.
#' @param topN Maximum number of top significant features to be plot for up- and
#' down-regulated genes. Default \code{20}.
#' @param absLFCThresh Hard threshold on absolute logFC value. Default \code{1}.
#' @param padjThresh Hard threshold on adjusted P-value. Default \code{0.05}.
#' @param pctInThresh,pctOutThresh Threshold on expression percentage. These
#' mean that a feature will only pass the filter if it is expressed in more than
#' \code{pctInThresh} percent of cells in the corresponding cluster. Similarly
#' for \code{pctOutThresh}. Only applied when these metrics are available.
#' Default \code{50} percent for both.
#' @param downsampleSize Maximum number of downsampled cells to be shown in the
#' heatmap. The downsampling is balanced on the cells involved in the test
#' specified. Default \code{200}.
#' @param useCellMeta Cell metadata variable names for cell grouping. Default
#' \code{NULL} includes dataset source and the default cluster.
#' @param column_title Title on the column. Default \code{NULL}.
#' @param seed Random seed for reproducibility. Default \code{1}.
#' @inheritDotParams .plotHeatmap transpose showCellLabel showCellLegend showFeatureLabel cellAnnColList featureAnnColList scale trim baseSize cellTextSize featureTextSize cellTitleSize featureTitleSize legendTextSize legendTitleSize viridisOption viridisDirection RColorBrewerOption
#' @return A \link[ComplexHeatmap]{HeatmapList-class} object.
#' @examples
#' defaultCluster(pbmc) <- pbmcPlot$leiden_cluster
#' pbmc <- normalize(pbmc)
#' plotPairwiseDEGHeatmap(pbmc, deg.pw, '4.stim')
plotPairwiseDEGHeatmap <- function(
        object,
        result,
        group = NULL,
        topN = 20,
        absLFCThresh = 1,
        padjThresh = 0.05,
        pctInThresh = 50,
        pctOutThresh = 50,
        downsampleSize = 200,
        useCellMeta = NULL,
        column_title = NULL,
        seed = 1,
        ...
) {
    resultMeta <- attr(result, 'groups')
    if (is.null(resultMeta)) {
        cli::cli_abort(
            c(x = "{.var result} object is corrupted, no meta-information available for selecting cells.",
              i = "Please re-run {.func runPairwiseDEG} with the latest {.pkg rliger} package (>2.0.1).")
        )
    }
    if (is.null(group)) {
        if (length(resultMeta) == 1) {
            group <- names(resultMeta)
        } else {
            cli::cli_abort(
                c("Please specify {.field group} when multiple tests are available",
                  i = "Available one{?s} {?is/are}: {.val {names(resultMeta)}}")
            )
        }
    }
    if (length(group) != 1) cli::cli_abort("Only 1 group allowed at a time.")
    if (!group %in% names(resultMeta)) {
        cli::cli_abort(
            c(x = "Group {.val {group}} not found in the result object.",
              i = "Available one{?s} {?is/are}: {.val {names(resultMeta)}}")
        )
    }
    cellIdx <- resultMeta[[group]]
    groupVar <- factor(rep(names(cellIdx), lengths(cellIdx)), levels = names(cellIdx))
    cellIdx <- Reduce(c, cellIdx)
    pctInThresh <- pctInThresh %||% 0
    pctOutThresh <- pctOutThresh %||% 100
    useCellMeta <- useCellMeta %||% c("dataset", object@uns$defaultCluster)
    # Idk why `filter(.data[['group']] == group)` won't work there lol
    result <- result[result$group == group,]
    result <- result %>%
        dplyr::filter(!is.na(.data[['padj']])) %>%
        dplyr::filter(
            abs(.data[['logFC']]) >= absLFCThresh,
            .data[['padj']] <= padjThresh
        ) %>%
        dplyr::mutate(regulation = factor(
            ifelse(.data[['logFC']] > 0, "up", "down"),
            levels = c("up", "down")
        ))
    if ("pct_in" %in% colnames(result) &&
        "pct_out" %in% colnames(result)) {
        result <- result %>%
            dplyr::filter(
                dplyr::case_when(
                    .data[['logFC']] > 0 ~ .data[['pct_in']] > pctInThresh & .data[['pct_out']] < pctOutThresh,
                    .data[['logFC']] < 0 ~ .data[['pct_out']] > pctInThresh & .data[['pct_in']] < pctOutThresh,
                )
            )
    }
    result <- result %>%
        dplyr::group_by(.data[['regulation']]) %>%
        dplyr::arrange(
            .data[['padj']],
            -.data[['logFC']],
            .by_group = TRUE
        ) %>%
        dplyr::slice_head(n = topN) %>%
        as.data.frame()

    set.seed(seed)
    downsampleSize <- min(length(cellIdx), downsampleSize)
    downsampleIdx <- seq_along(cellIdx) %in% sample(length(cellIdx), downsampleSize)
    cellIdxSub <- cellIdx[downsampleIdx]
    groupVar <- groupVar[downsampleIdx]

    cellAnn <- data.frame(
        group = groupVar,
        row.names = colnames(object)[cellIdxSub]
    )
    featureAnn <- data.frame(
        regulation = result$regulation,
        row.names = result$feature
    )
    plotGeneHeatmap(object, features = result$feature,
                    cellIdx = cellIdxSub,
                    useCellMeta = useCellMeta,
                    cellAnnotation = cellAnn,
                    featureAnnotation = featureAnn,
                    cellSplitBy = "group",
                    featureSplitBy = "regulation",
                    showFeatureLegend = FALSE,
                    cluster_columns = FALSE,
                    cluster_column_slices = FALSE,
                    cluster_rows = FALSE,
                    cluster_row_slices = FALSE,
                    column_title = column_title,
                    ...)
}

