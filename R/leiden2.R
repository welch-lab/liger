#' @export
#' @rdname clustering
runLeidenCluster2 <- function(
        object,
        resolution = 1.0,
        nNeighbors = 20,
        partitionType = c("ModularityVertexPartition", "CPMVertexPartition",
                          "RBConfigurationVertexPartition", "RBERVertexPartition",
                          "SignificanceVertexPartition", "SurpriseVertexPartition"),
        prune = 1 / 15,
        eps = 0.1,
        nRandomStarts = 10,
        nIterations = 5,
        useDims = NULL,
        groupSingletons = TRUE,
        clusterName = "leiden_cluster",
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    partitionType <- match.arg(partitionType)
    object <- recordCommand(object,
                            dependencies = c("RANN", "igraph"))
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
    clusts <- leidenalg_find_partition(
        snn, nrow(H.norm), resolution = resolution, nStart = nRandomStarts,
        nIter = nIterations, seed = seed, verbose = verbose
    )
    names(clusts) <- colnames(object)
    rownames(snn) <- colnames(object)
    colnames(snn) <- colnames(object)
    clusts <- groupSingletons(ids = clusts, SNN = snn,
                              groupSingletons = groupSingletons,
                              verbose = verbose)
    if (is.numeric(clusts)) clusts <- clusts - 1
    clusts <- .labelClustBySize(clusts)
    names(clusts) <- colnames(object)
    # clusts <- factor(clusts)
    cellMeta(object, clusterName, check = FALSE) <- clusts
    return(object)
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


leidenalg_find_partition <- function(snn, nCells, resolution=1.0, nStart=10, nIter = 10, seed = 123, verbose = TRUE) {
    snnSummary <- summary(snn)
    edgelist <- as.vector(t(snnSummary[, c(2, 1)])) - 1
    edgelist_length <- length(edgelist)
    edge_weights <- snnSummary[,3]
    res <- leidenalg_find_partition_rcpp(
        edgelist = edgelist, edgelist_length = edgelist_length,
        num_vertices = nCells, direction = FALSE,
        edge_weights = edge_weights, resolution = resolution, nStart = nStart,
        nIter = nIter, seed = seed, verbose = verbose)
    return(res)
}
