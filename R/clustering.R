#' SNN Graph Based Community Detection
#' @description
#' After quantile normalization, users can additionally run the Leiden or
#' Louvain algorithm for community detection, which is widely used in
#' single-cell analysis and excels at merging small clusters into broad cell
#' classes.
#'
#' While using quantile normalized factor loadings (result from
#' \code{\link{quantileNorm}}) is recommended, this function looks for
#' unnormalized factor loadings (result from \code{\link{runIntegration}}) when
#' the former is not available.
#' @param object A \linkS4class{liger} object. Should have valid factorization
#' result available.
#' @param nNeighbors Integer, the maximum number of nearest neighbors to
#' compute. Default \code{20}.
#' @param resolution Numeric, value of the resolution parameter, a larger value
#' results in a larger number of communities with smaller sizes. Default
#' \code{1.0}.
#' @param prune Numeric. Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the stringency of pruning. \code{0} for no pruning,
#' while \code{1} prunes everything. Default \code{1/15}.
#' @param eps Numeric, the error bound of the nearest neighbor search. Default
#' \code{0.1}.
#' @param nRandomStarts Integer number of random starts. Will pick the
#' membership with highest quality to return. Default \code{10}.
#' @param nIterations Integer, maximal number of iterations per random start.
#' Default \code{5}.
#' @param method Community detection algorithm to use. Choose from
#' \code{"leiden"} or \code{"louvain"}. Default \code{"leiden"}.
#' @param useRaw Whether to use un-aligned cell factor loadings (\eqn{H}
#' matrices). Default \code{NULL} search for quantile-normalized loadings first
#' and un-aligned loadings then.
#' @param useDims Indices of factors to use for clustering. Default \code{NULL}
#' uses all available factors.
#' @param groupSingletons Whether to group single cells that make up their own
#' cluster in with the cluster they are most connected to. Default \code{TRUE},
#' if \code{FALSE}, assign all singletons to a \code{"singleton"} group.
#' @param clusterName Name of the variable that will store the clustering result
#' in \code{cellMeta} slot of \code{object}. Default \code{"leiden_cluster"} and
#' \code{"louvain_cluster"}.
#' @param seed Seed of the random number generator. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return \code{object} with refined cluster assignment updated in
#' \code{clusterName} variable in \code{cellMeta} slot. Can be fetched
#' with \code{object[[clusterName]]}
#' @rdname runCluster
#' @export
#' @examples
#' pbmcPlot <- runCluster(pbmcPlot)
#' head(pbmcPlot$leiden_cluster)
#' pbmcPlot <- runCluster(pbmcPlot, method = "louvain")
#' head(pbmcPlot$louvain_cluster)
runCluster <- function(
        object,
        resolution = 1.0,
        nNeighbors = 20,
        prune = 1 / 15,
        eps = 0.1,
        nRandomStarts = 10,
        nIterations = 5,
        method = c("leiden", "louvain"),
        useRaw = NULL,
        useDims = NULL,
        groupSingletons = TRUE,
        clusterName = paste0(method, "_cluster"),
        seed = 1,
        verbose = getOption("ligerVerbose", TRUE)
) {
    method <- match.arg(method)
    object <- switch(method,
        leiden = recordCommand(object, dependencies = c("RANN", "leidenAlg")),
        louvain = recordCommand(object, dependencies = c("RANN"))
    )
    Hsearch <- searchH(object, useRaw)
    H <- Hsearch$H
    useRaw <- Hsearch$useRaw
    type <- ifelse(useRaw, "unnormalized", "quantile normalized")

    if (!is.null(useDims)) H <- H[, useDims, drop = FALSE]

    if (isTRUE(verbose))
        cli::cli_alert_info("{method} clustering on {type} cell factor loadings...")
    knn <- RANN::nn2(H, k = nNeighbors, eps = eps)
    snn <- ComputeSNN(knn$nn.idx, prune = prune)
    if (!is.null(seed)) set.seed(seed)
    if (method == "leiden") {
        snnSummary <- summary(snn)
        edgelist <- as.vector(t(snnSummary[, c(2, 1)])) - 1
        edgelist_length <- length(edgelist)
        edge_weights <- snnSummary[,3]
        clusts <- leidenAlg::find_partition_with_rep_rcpp(
            edgelist = edgelist, edgelist_length = edgelist_length,
            num_vertices = nrow(H), direction = FALSE,
            edge_weights = edge_weights, resolution = resolution,
            niter = nIterations, nrep = nRandomStarts
        )
    } else {
        edgeOutPath <- paste0("edge_", sub("\\s", "_", Sys.time()), '.txt')
        edgeOutPath <- gsub("-", "", edgeOutPath)
        edgeOutPath <- gsub(":", "", edgeOutPath)
        WriteEdgeFile(snn, edgeOutPath, display_progress = FALSE)
        clusts <- RunModularityClusteringCpp(
            snn,
            modularityFunction = 1,
            resolution = resolution,
            nRandomStarts = nRandomStarts,
            nIterations = nIterations,
            algorithm = 1,
            randomSeed = seed,
            printOutput = verbose,
            edgefilename = edgeOutPath
        )
        unlink(edgeOutPath)
    }
    clusts <- .labelClustBySize(clusts)
    names(clusts) <- colnames(object)
    rownames(snn) <- colnames(object)
    colnames(snn) <- colnames(object)
    clusts <- groupSingletons(ids = clusts, SNN = snn,
                              groupSingletons = groupSingletons,
                              verbose = verbose)
    cellMeta(object, clusterName, check = FALSE) <- clusts
    if (isTRUE(verbose))
        cli::cli_alert_success("Found {nlevels(clusts)} clusters.")
    object@uns$defaultCluster <- object@uns$defaultCluster %||% clusterName
    if (isTRUE(verbose))
        cli::cli_alert_info("cellMeta variable {.val {clusterName}} is now set as default.")
    return(object)
}

#' [Deprecated] Louvain algorithm for community detection
#' @description
#' After quantile normalization, users can additionally run the Louvain
#' algorithm for community detection, which is widely used in single-cell
#' analysis and excels at merging small clusters into broad cell classes.
#' @param object \code{liger} object. Should run quantile_norm before calling.
#' @param k The maximum number of nearest neighbours to compute. (default 20)
#' @param resolution Value of the resolution parameter, use a value above
#' (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' (default 1.0)
#' @param prune Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything). (default 1/15)
#' @param eps The error bound of the nearest neighbor search. (default 0.1)
#' @param nRandomStarts Number of random starts. (default 10)
#' @param nIterations Maximal number of iterations per random start. (default
#' 100)
#' @param random.seed Seed of the random number generator. (default 1)
#' @param verbose Print messages (TRUE by default)
#' @param dims.use Indices of factors to use for clustering. Default \code{NULL}
#' uses all available factors.
#' @return \code{object} with refined cluster assignment updated in
#' \code{"louvain_cluster"} variable in \code{cellMeta} slot. Can be fetched
#' with \code{object$louvain_cluster}
#' @name louvainCluster-deprecated
#' @seealso \code{\link{rliger2-deprecated}}
NULL

#' @rdname rliger2-deprecated
#' @section \code{louvainCluster}:
#' For \code{louvainCluster}, use \code{\link{runCluster}(method = "louvain")}
#' as the replacement, while \code{\link{runCluster}} with default
#' \code{method = "leiden"} is more recommended.
#' @export
louvainCluster <- function( # nocov start
        object,
        resolution = 1.0,
        k = 20,
        prune = 1 / 15,
        eps = 0.1,
        nRandomStarts = 10,
        nIterations = 100,
        random.seed = 1,
        verbose = getOption("ligerVerbose"),
        dims.use = NULL
) {
    lifecycle::deprecate_warn(
        "1.99.0", "louvainCluster()",
        details = "Please use `runCluster()` with `method = 'louvain'` instead.")
    runCluster(
        object, method = "louvain", resolution = resolution, nNeighbors = k,
        prune = prune, eps = eps, nRandomStarts = nRandomStarts,
        nIterations = nIterations, useDims = dims.use, groupSingletons = TRUE,
        clusterName = "louvain_cluster", seed = random.seed, verbose = verbose
    )
} # nocov end

# Group single cells that make up their own cluster in with the cluster they are
# most connected to. (Adopted from Seurat v3)
#
# @param ids Named vector of cluster ids
# @param SNN SNN graph used in clustering
# @param groupSingletons Group singletons into nearest cluster (TRUE by
# default). If FALSE, assign all singletons to a "singleton" group.
# @param verbose Print message
# @return Returns updated cluster assignment (vector) with all singletons merged
# with most connected cluster
groupSingletons <- function(
        ids,
        SNN,
        groupSingletons = TRUE,
        verbose = FALSE
) {
    # identify singletons
    ids <- droplevels(ids)
    singletons <- names(which(table(ids) == 1))
    singletons <- intersect(levels(ids), singletons)
    if (length(singletons) == 0) {
        return(ids)
    }
    if (!isTRUE(groupSingletons)) {
        if (isTRUE(verbose)) {
            cli::cli_alert_info("{length(singletons)} singletons identified.")
        }
        ids <- as.character(ids)
        ids[ids %in% singletons] <- "singleton"
        return(factor(ids))
    }
    # calculate connectivity of singletons to other clusters, add singleton
    # to cluster it is most connected to
    clusterNames <- levels(ids)
    clusterNames <- setdiff(clusterNames, singletons)
    connectivity <- vector(mode = "numeric", length = length(clusterNames))
    names(connectivity) <- clusterNames
    ids <- as.character(ids)
    for (i in singletons) {
        i.cells <- which(ids == i)
        for (j in clusterNames) {
            j.cells <- which(ids == j)
            subSNN <- SNN[i.cells, j.cells, drop = FALSE]
            # to match previous behavior, random seed being set in WhichCells
            set.seed(1)
            if (is.object(subSNN))
                connectivity[j] <- sum(subSNN) / (nrow(subSNN) * ncol(subSNN))
            else
                connectivity[j] <- mean(subSNN)
        }
        m <- max(connectivity, na.rm = TRUE)
        mi <- which(connectivity == m, arr.ind = TRUE)
        closest_cluster <- sample(names(connectivity[mi]), 1)
        ids[i.cells] <- closest_cluster
    }
    ids <- factor(ids)
    if (isTRUE(verbose))
        cli::cli_alert_info("{length(singletons)} singletons identified.")
    return(ids)
}

.labelClustBySize <- function(clusts) {
    clusts <- as.character(clusts)
    count <- data.frame(table(clusts))
    count <- count[order(count$Freq, decreasing = TRUE),]
    map <- seq(0, nrow(count) - 1)
    names(map) <- count[[1]]
    newClusts <- map[clusts]
    newClusts <- factor(newClusts)
    return(newClusts)
}



#' Create new variable from categories in cellMeta
#' @description
#' Designed for fast variable creation when a new variable is going to be
#' created from existing variable. For example, multiple samples can be mapped
#' to the same study design condition, clusters can be mapped to cell types.
#' @param object A \linkS4class{liger} object.
#' @param from The name of the original variable to be mapped from.
#' @param newTo The name of the new variable to store the mapped result. Default
#' \code{NULL} returns the new variable (factor class).
#' @param ... Mapping criteria, argument names are original existing categories
#' in the \code{from} and values are new categories in the new variable.
#' @return When \code{newTo = NULL}, a factor object of the new variable.
#' Otherwise, the input object with variable \code{newTo} updated in
#' \code{cellMeta(object)}.
#' @export
#' @examples
#' pbmc <- mapCellMeta(pbmc, from = "dataset", newTo = "modal",
#'                     ctrl = "rna", stim = "rna")
mapCellMeta <- function(
        object,
        from,
        newTo = NULL,
        ...
) {
    object <- recordCommand(object, ...)
    from <- cellMeta(object, from)
    if (!is.factor(from))
        cli::cli_abort("{.code from} must be a {.cls factor}.")
    mapping <- list(...)
    fromCats <- names(mapping)
    notFound <- fromCats[!fromCats %in% levels(from)]
    if (length(notFound) > 0) {
        cli::cli_abort("{length(notFound)} categor{?y is/ies are} requested but not found: {.val {notFound}}")
    }

    toCats <- unlist(mapping)
    unchangedCats <- levels(from)
    unchangedCats <- unchangedCats[!unchangedCats %in% fromCats]
    names(unchangedCats) <- unchangedCats
    if (length(unchangedCats) > 0) toCats <- c(toCats, unchangedCats)
    to <- toCats[as.character(from)]
    to <- factor(unname(to), levels = unique(toCats))
    if (is.null(newTo)) return(to)
    cellMeta(object, newTo) <- to
    return(object)
}
