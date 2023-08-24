#' SNN Graph Based Community Detection
#' @description After quantile normalization, users can additionally run the
#' Leiden or Louvain algorithm for community detection, which is widely used in
#' single-cell analysis and excels at merging small clusters into broad cell
#' classes.
#'
#' While using quantile normalized factor loadings (result from
#' \code{\link{quantileNorm}}) is recommended, this function looks for
#' unnormalized factor loadings (result from \code{\link{optimizeALS}} or
#' \code{\link{online_iNMF}}) when the former is not available.
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
#' Default \code{100}.
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
#' @param partitionType Choose from \code{"ModularityVertexPartition",
#' "RBConfigurationVertexPartition", "RBERVertexPartition",
#' "SignificanceVertexPartition", "CPMVertexPartition",
#' "SurpriseVertexPartition"}. See
#' \code{\link[leidenbase]{leiden_find_partition}} for detail. Default
#' \code{"ModularityVertexPartition"}.
#' @param ... Additional arguments passed to
#' \code{\link[leidenbase]{leiden_find_partition}}, including
#' \code{initial_membership} and \code{node_sizes}.
#' @return \code{object} with refined cluster assignment updated in
#' \code{clusterName} variable in \code{cellMeta} slot. Can be fetched
#' with \code{object[[clusterName]]}
#' @rdname clustering
#' @export
#' @examples
#' pbmcPlot <- runLeidenCluster(pbmcPlot)
#' head(pbmcPlot$leiden_cluster)
#' pbmcPlot <- runLouvainCluster(pbmcPlot)
#' head(pbmcPlot$louvain_cluster)
runLeidenCluster <- function(
        object,
        resolution = 1.0,
        nNeighbors = 20,
        partitionType = c("ModularityVertexPartition",
                          "RBConfigurationVertexPartition",
                          "RBERVertexPartition",
                          "SignificanceVertexPartition",
                          "CPMVertexPartition",
                          "SurpriseVertexPartition"),
        prune = 1 / 15,
        eps = 0.1,
        nRandomStarts = 10,
        nIterations = 100,
        useDims = NULL,
        groupSingletons = TRUE,
        clusterName = "leiden_cluster",
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    partitionType <- match.arg(partitionType)
    object <- recordCommand(object,
                            dependencies = c("RANN", "leidenbase", "igraph"))
    H.norm <- getMatrix(object, "H.norm")
    if (is.null(H.norm)) {
        type <- " unnormalized "
        H.norm <- Reduce(cbind, getMatrix(object, "H"))
    } else type <- " quantile normalized "
    if (is.null(H.norm))
        stop("No factor loading ('H.norm' or 'H') found in `object`.")
    if (type == " unnormalized ") H.norm <- t(H.norm)

    if (!is.null(useDims)) H.norm <- H.norm[, useDims]

    if (isTRUE(verbose))
        .log("Leiden clustering on", type, "cell factor loadings...")
    knn <- RANN::nn2(H.norm, k = nNeighbors, eps = eps)
    snn <- ComputeSNN(knn$nn.idx, prune = prune)
    snnSummary <- summary(snn)
    edges <- as.vector(t(snnSummary[,c(1, 2)]))
    g <- igraph::graph(edges = edges, n = nrow(snn),directed = FALSE)
    igraph::E(g)$weight <- snnSummary[, 3]

    set.seed(seed)
    maxQuality <- -1
    if (isTRUE(verbose))
        pb <- utils::txtProgressBar(0, nRandomStarts, style = 3)
    for (i in seq(nRandomStarts)) {
        seed <- sample(1000, 1)
        part <- leidenbase::leiden_find_partition(
            igraph = g,
            partition_type = partitionType,
            resolution_parameter = resolution,
            edge_weights = igraph::E(g)$weight,
            num_iter = nIterations,
            verbose = FALSE,
            seed = seed,
            ...
        )
        #if (is.null(part$quality)) {
        #    clusts <- part$membership
        #    if (isTRUE(verbose)) utils::setTxtProgressBar(pb, nRandomStarts)
        #    break
        #}
        if (part$quality > maxQuality) {
            clusts <- part$membership
            maxQuality <- part$quality
        }
        if (isTRUE(verbose)) utils::setTxtProgressBar(pb, i)
    }
    if (isTRUE(verbose)) cat("\n")
    names(clusts) <- colnames(object)
    rownames(snn) <- colnames(object)
    colnames(snn) <- colnames(object)

    clusts <- groupSingletons(ids = clusts, SNN = snn,
                              groupSingletons = groupSingletons,
                              verbose = verbose)
    if (is.numeric(clusts)) clusts <- clusts - 1
    clusts <- factor(clusts)
    cellMeta(object, clusterName, check = FALSE) <- clusts
    return(object)
}

#' @rdname clustering
#' @export
runLouvainCluster <- function(
        object,
        resolution = 1.0,
        nNeighbors = 20,
        prune = 1 / 15,
        eps = 0.1,
        nRandomStarts = 10,
        nIterations = 100,
        useDims = NULL,
        groupSingletons = TRUE,
        clusterName = "louvain_cluster",
        seed = 1,
        verbose = getOption("ligerVerbose")
) {
    object <- recordCommand(object, dependencies = "RANN")
    H.norm <- getMatrix(object, "H.norm")
    if (is.null(H.norm)) {
        type <- " unnormalized "
        H.norm <- Reduce(cbind, getMatrix(object, "H"))
    } else type <- " quantile normalized "
    if (is.null(H.norm))
        stop("No factor loading ('H.norm' or 'H') found in `object`.")
    # Not transposing when cbind'ing becausing `t(NULL)` causes error
    if (type == " unnormalized ") H.norm <- t(H.norm)

    edgeOutPath <- paste0("edge_", sub("\\s", "_", Sys.time()), '.txt')
    edgeOutPath <- gsub("-", "", edgeOutPath)
    edgeOutPath <- gsub(":", "", edgeOutPath)

    if (!is.null(useDims)) H.norm <- H.norm[, useDims]

    if (isTRUE(verbose))
        .log("Louvain clustering on", type, "cell factor loadings...")
    knn <- RANN::nn2(H.norm, k = nNeighbors, eps = eps)
    snn <- ComputeSNN(knn$nn.idx, prune = prune)
    WriteEdgeFile(snn, edgeOutPath, display_progress = FALSE)
    clusts <- RunModularityClusteringCpp(
        snn,
        modularityFunction = 1,
        resolution = resolution,
        nRandomStarts = nRandomStarts,
        nIterations = nIterations,
        algorithm = 1,
        randomSeed = seed,
        printOutput = TRUE,
        edgefilename = edgeOutPath
    )
    names(clusts) <- colnames(object)
    rownames(snn) <- colnames(object)
    colnames(snn) <- colnames(object)
    # clusts must not be a factor at this point
    clusts <- groupSingletons(ids = clusts, SNN = snn,
                              groupSingletons = groupSingletons,
                              verbose = verbose)
    cellMeta(object, columns = clusterName, check = FALSE) <- factor(clusts)
    unlink(edgeOutPath)
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
#' For \code{louvainCluster}, use \code{\link{runLouvainCluster}} as the
#' replacement, while \code{\link{runLeidenCluster}} is more recommended.
#' @export
louvainCluster <- function(
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
    lifecycle::deprecate_warn("1.99.0", "louvainCluster()",
                              "runLouvainCluster()")
    runLouvainCluster(
        object, resolution = resolution, nNeighbors = k, prune = prune,
        eps = eps, nRandomStarts = nRandomStarts, nIterations = nIterations,
        useDims = dims.use, groupSingletons = TRUE,
        clusterName = "louvain_cluster", seed = random.seed, verbose = verbose
    )
}

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
    singletons <- names(which(table(ids) == 1))
    singletons <- intersect(unique(ids), singletons)
    if (!isTRUE(groupSingletons)) {
        ids[ids %in% singletons] <- "singleton"
        return(ids)
    }
    # calculate connectivity of singletons to other clusters, add singleton
    # to cluster it is most connected to
    clusterNames <- unique(ids)
    clusterNames <- setdiff(clusterNames, singletons)
    connectivity <- vector(mode = "numeric", length = length(clusterNames))
    names(connectivity) <- clusterNames
    for (i in singletons) {
        i.cells <- names(which(ids == i))
        for (j in clusterNames) {
            j.cells <- names(which(ids == j))
            subSNN <- SNN[i.cells, j.cells]
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
    if (length(singletons) > 0 && isTRUE(verbose))
        .log(length(singletons), " singletons identified. ",
             length(unique(ids)), " final clusters.")
    return(ids)
}
