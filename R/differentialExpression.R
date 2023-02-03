#' Perform Wilcoxon rank-sum test
#' @description Perform Wilcoxon rank-sum tests on specified dataset using
#' given method.
#' @param object A \linkS4class{liger} object with cluster assignment available.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be normalized. Default
#' \code{NULL} uses all datasets.
#' @param method Choose from \code{"clusters"} or \code{"datasets"}. Default
#' \code{"clusters"} compares between clusters across all datasets, while
#' \code{"datasets"} compares between datasets within each cluster.
#' @param useCluster The name of the column in \code{cell.meta} slot storing the
#' cluster assignment variable. Default \code{"louvain_cluster"}
#' @param verbose Logical. Whether to show information of the progress.
#' Default \code{TRUE}.
#' @param data.use,compare.method \bold{Deprecated}. See Usage section for
#' replacement.
#' @return A 10-columns data.frame with test results.
#' @export
runWilcoxon <- function(
        object,
        useDatasets = NULL,
        method = c("clusters", "datasets"),
        useCluster = "louvain_cluster",
        verbose = TRUE,
        # Deprecated coding style,
        data.use = useDatasets,
        compare.method = method
) {
    .deprecateArgs(list(data.use = "useDatasets", compare.method = "method"),
                   call = rlang::call_args(match.call()))
    method <- match.arg(method)
    # Input checks ####
    useDatasets <- .checkUseDatasets(object, useDatasets)
    if (!useCluster %in% colnames(cell.meta(object)))
        # TODO: other clustering methods?
        stop("Cannot find `", useCluster, "` cluster assignment in `cell.meta`",
             " please run `louvainCluster()` in advance.")
    if (method == "datasets" & length(useDatasets) < 2)
        stop("Should have at least 2 datasets as input ",
             "when compare between datasets")
    if (isH5Liger(object, useDatasets)) {
        stop("HDF5 based datasets detected but is not supported. \n",
             "Try `object.sub <- downsample(object, useSlot = 'norm.data')`",
             " to create ANOTHER object with in memory data.")
    }
    allNormed <- all(sapply(datasets(object),
                            function(ld) !is.null(norm.data(ld))))
    if (!allNormed)
        stop("All datasets being involved has to be normalized")

    if (isTRUE(verbose))
        .log("Performing Wilcoxon test on ", length(useDatasets), " datasets: ",
             paste(useDatasets, collapse = ", "))

    # Data preperation ####
    ## get all shared genes of datasets involved
    geneNameList <- lapply(datasets(object)[useDatasets], function(ld) {
        rownames(ld)
    })
    genes <- Reduce(intersect, geneNameList)
    normDataList <- getMatrix(object, "norm.data", dataset = useDatasets)
    normDataList <- lapply(normDataList, function(x) x[genes,])
    featureMatrix <- Reduce(cbind, normDataList)
    ## Subset metadata involved
    cellIdx <- object$dataset %in% useDatasets
    cellSource <- object$dataset[cellIdx]
    clusters <- cell.meta(object)[[useCluster]][cellIdx]

    # perform wilcoxon test ####
    if (method == "clusters") {
        # compare between clusters across datasets
        nGenes <- length(genes)
        if (nGenes > 100000) {
            if (isTRUE(verbose)) .log("Calculating Large-scale Input...")
            results <- Reduce(rbind, lapply(
                suppressWarnings(split(seq(nGenes), seq(nGenes / 100000))),
                function(index) {
                    wilcoxauc(log(featureMatrix[index, ] + 1e-10), clusters)
                }))
        } else {
            results <- wilcoxauc(log(featureMatrix + 1e-10), clusters)
        }
    } else {
        # compare between datasets within each cluster
        results <- Reduce(rbind, lapply(levels(clusters), function(cluster) {
            clusterIdx <- clusters == cluster
            subLabel <- paste0(cluster, "-", cellSource[clusterIdx])
            if (length(unique(subLabel)) == 1) {
                # if cluster has only 1 data source
                warning("Skipped Cluster ", cluster,
                        " since it has only one dataset source.")
                return()
            } else {
                subMatrix <- featureMatrix[, clusterIdx]
                return(wilcoxauc(log(subMatrix + 1e-10), subLabel))
            }
        }))
    }
    return(results)
}
