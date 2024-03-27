#' Perform UMAP Dimensionality Reduction
#' @description
#' Run UMAP on the quantile normalized cell factors (result from
#' \code{\link{quantileNorm}}), or unnormalized cell factors (result from
#' \code{\link{runIntegration}})) to generate a 2D embedding for visualization
#' (or general dimensionality reduction). Has option to run on subset of
#' factors. It is generally recommended to use this method for dimensionality
#' reduction with extremely large datasets. The underlying UMAP calculation
#' imports uwot \code{\link[uwot]{umap}}.
#' @details For \code{nNeighbors}, larger values will result in more global
#' structure being preserved at the loss of detailed local structure. In general
#' this parameter should often be in the range 5 to 50, with a choice of 10 to
#' 15 being a sensible default.
#'
#' For \code{minDist}, larger values ensure embedded points are more evenly
#' distributed, while smaller values allow the algorithm to optimize more
#' accurately with regard to local structure. Sensible values are in the range
#' 0.001 to 0.5, with 0.1 being a reasonable default.
#' @param object \linkS4class{liger} object with factorization results.
#' @param useRaw Whether to use un-aligned cell factor loadings (\eqn{H}
#' matrices). Default \code{NULL} search for quantile-normalized loadings first
#' and un-aligned loadings then.
#' @param useDims Index of factors to use for computing UMAP embedding. Default
#' \code{NULL} uses all factors.
#' @param nDims Number of dimensions to reduce to. Default \code{2}.
#' @param distance Character. Metric used to measure distance in the input
#' space. Default \code{"cosine"}, alternative options include:
#' \code{"euclidean"}, \code{"manhattan"} and \code{"hamming"}.
#' @param nNeighbors Number of neighboring points used in local approximations
#' of manifold structure. Default \code{10}.
#' @param minDist Numeric. Controls how tightly the embedding is allowed
#' compress points together. Default \code{0.1}.
#' @param dimredName Name of the variable in \code{cellMeta} slot to store the
#' result matrix. Default \code{"UMAP"}.
#' @param seed Random seed for reproducibility. Default \code{42}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param k,use.raw,dims.use,n_neighbors,min_dist,rand.seed \bold{Deprecated}.
#' See Usage section for replacement.
#' @return The \code{object} where a \code{"UMAP"} variable is updated in the
#' \code{cellMeta} slot with the whole 2D embedding matrix.
#' @seealso \code{\link{runTSNE}}
#' @export
#' @examples
#' pbmc <- runUMAP(pbmcPlot)
runUMAP <- function(
        object,
        useRaw = NULL,
        useDims = NULL,
        nDims = 2,
        distance = c("cosine", "euclidean", "manhattan", "hamming"),
        nNeighbors = 20,
        minDist = 0.1,
        dimredName = "UMAP",
        seed = 42,
        verbose = getOption("ligerVerbose", TRUE),
        # Deprecated coding style
        k = nDims,
        use.raw = useRaw,
        dims.use = useDims,
        n_neighbors = nNeighbors,
        min_dist = minDist,
        rand.seed = seed
) {
    .deprecateArgs(list(k = "nDims", use.raw = "useRaw", dims.use = "useDims",
                        n_neighbors = "nNeighbors", min_dist = "minDist",
                        rand.seed = "seed"))
    distance <- match.arg(distance)
    object <- recordCommand(object, dependencies = "uwot")
    set.seed(seed)
    Hsearch <- searchH(object, useRaw)
    H <- Hsearch$H
    useRaw <- Hsearch$useRaw
    type <- ifelse(useRaw, "unnormalized", "quantile normalized")
    if (isTRUE(verbose))
        cli::cli_process_start("Generating UMAP on {type} cell factor loadings")
    if (!is.null(useDims)) H <- H[, useDims, drop = FALSE]
    umap <- uwot::umap(H,
                       n_components = as.integer(nDims),
                       metric = distance,
                       n_neighbors = as.integer(nNeighbors),
                       min_dist = minDist)
    if (isTRUE(verbose)) cli::cli_process_done()
    dimRed(object, dimredName) <- umap
    if (isTRUE(verbose))
        cli::cli_alert_info("{.field DimRed} {.val {dimredName}} is now set as default.")
    return(object)
}

#' Perform t-SNE dimensionality reduction
#' @description
#' Runs t-SNE on the quantile normalized cell factors (result from
#' \code{\link{quantileNorm}}), or unnormalized cell factors (result from
#' \code{\link{runIntegration}})) to generate a 2D embedding for visualization.
#' By default \code{\link[Rtsne]{Rtsne}} (Barnes-Hut implementation of t-SNE)
#' method is invoked, while alternative "fftRtsne" method (FFT-accelerated
#' Interpolation-based t-SNE, using Kluger Lab implementation) is also
#' supported. For very large datasets, it is recommended to use
#' \code{method = "fftRtsne"} due to its efficiency and scalability.
#'
#' Extra external installation steps are required for using "fftRtsne" method.
#' Please consult
#' \href{https://welch-lab.github.io/liger/articles/installation.html}{detailed guide}.
#' @param object \linkS4class{liger} object with factorization results.
#' @param useRaw Whether to use un-aligned cell factor loadings (\eqn{H}
#' matrices). Default \code{NULL} search for quantile-normalized loadings first
#' and un-aligned loadings then.
#' @param useDims Index of factors to use for computing UMAP embedding. Default
#' \code{NULL} uses all factors.
#' @param nDims Number of dimensions to reduce to. Default \code{2}.
#' @param usePCA Whether to perform initial PCA step for Rtsne. Default
#' \code{FALSE}.
#' @param perplexity Numeric parameter to pass to Rtsne (expected number of
#' neighbors). Default \code{30}.
#' @param theta Speed/accuracy trade-off (increase for less accuracy), set to
#' \code{0.0} for exact TSNE. Default \code{0.5}.
#' @param method Choose from \code{"Rtsne"} or \code{"fftRtsne"}. See
#' Description. Default \code{"Rtsne"}.
#' @param dimredName Name of the variable in \code{cellMeta} slot to store the
#' result matrix. Default \code{"TSNE"}.
#' @param fitsnePath Path to the cloned FIt-SNE directory (i.e.
#' \code{'/path/to/dir/FIt-SNE'}). Required only when first time using
#' \code{runTSNE} with \code{method = "fftRtsne"}. Default \code{NULL}.
#' @param seed Random seed for reproducibility. Default \code{42}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param use.raw,dims.use,k,use.pca,fitsne.path,rand.seed \bold{Deprecated}.
#' See Usage section for replacement.
#' @return The \code{object} where a \code{"TSNE"} variable is updated in the
#' \code{cellMeta} slot with the whole 2D embedding matrix.
#' @seealso \code{\link{runUMAP}}
#' @export
#' @examples
#' pbmc <- runTSNE(pbmcPlot)
runTSNE <- function(
        object,
        useRaw = NULL,
        useDims = NULL,
        nDims = 2,
        usePCA = FALSE,
        perplexity = 30,
        theta = 0.5,
        method = c("Rtsne", "fftRtsne"),
        dimredName = "TSNE",
        fitsnePath = NULL,
        seed = 42,
        verbose = getOption("ligerVerbose", TRUE),
        # Deprecated coding styles
        k = nDims,
        use.raw = useRaw,
        dims.use = useDims,
        use.pca = usePCA,
        fitsne.path = fitsnePath,
        rand.seed = seed
) {
    .deprecateArgs(list(k = "nDims", use.raw = "useRaw", dims.use = "useDims",
                        use.pca = "usePCA", fitsne.path = "fitsnePath",
                        rand.seed = "seed"))
    method <- match.arg(method)
    object <- recordCommand(object, dependencies = "Rtsne")
    Hsearch <- searchH(object, useRaw)
    H <- Hsearch$H
    useRaw <- Hsearch$useRaw
    type <- ifelse(useRaw, "unnormalized", "quantile normalized")
    if (isTRUE(verbose))
        cli::cli_process_start("Generating TSNE ({method}) on {type} cell factor loadings")
    if (!is.null(useDims)) H <- H[, useDims, drop = FALSE]
    if (method == "Rtsne") {
        set.seed(seed)
        tsne <- Rtsne::Rtsne(H,
                             dims = nDims,
                             pca = usePCA,
                             check_duplicates = FALSE,
                             theta = theta,
                             perplexity = perplexity)
        tsne <- tsne$Y
    } else if (method == "fftRtsne") {
        tsne <- .fftRtsne(H,
                          dims = nDims,
                          rand_seed = seed,
                          fast_tsne_path = fitsnePath,
                          theta = theta,
                          perplexity = perplexity)
    }
    if (isTRUE(verbose)) cli::cli_process_done()
    dimRed(object, dimredName) <- tsne
    object@uns$TSNE <- list(method = method)
    if (isTRUE(verbose))
        cli::cli_alert_info("{.field DimRed} {.val {dimredName}} is now set as default.")
    return(object)
}

# FIt-SNE helper function for calling fast_tsne from R
#
# Codes from Seurat (https://github.com/satijalab/seurat)
# Originally Based on Kluger Lab FIt-SNE v1.1.0 code on
# https://github.com/KlugerLab/FIt-SNE/blob/master/fast_tsne.R
# commit d2cf403 on Feb 8, 2019
.fftRtsne <- function(
        X,
        dims = 2,
        perplexity = 30,
        theta = 0.5,
        check_duplicates = TRUE,
        max_iter = 1000,
        fft_not_bh = TRUE,
        ann_not_vptree = TRUE,
        stop_early_exag_iter = 250,
        exaggeration_factor = 12.0,
        no_momentum_during_exag = FALSE,
        start_late_exag_iter = -1.0,
        late_exag_coeff = 1.0,
        mom_switch_iter = 250,
        momentum = 0.5,
        final_momentum = 0.8,
        learning_rate = 200,
        n_trees = 50,
        search_k = -1,
        rand_seed = -1,
        nterms = 3,
        intervals_per_integer = 1,
        min_num_intervals = 50,
        K = -1,
        sigma = -30,
        initialization = NULL,
        data_path = NULL,
        result_path = NULL,
        load_affinities = NULL,
        fast_tsne_path = NULL,
        nthreads = getOption("mc.cores", default = 1),
        perplexity_list = NULL,
        get_costs = FALSE,
        df = 1.0,
        ...
) {
    if (is.null(fast_tsne_path)) {
        stop("Please pass in path to FIt-SNE directory as fitsne.path.")
    } else {
        if (.Platform$OS.type == "unix") {
            fast_tsne_path <- file.path(fast_tsne_path, "bin", "fast_tsne")
        } else {
            fast_tsne_path <- file.path(fast_tsne_path, "bin", "FItSNE.exe")
        }
    }

    if (is.null(data_path)) {
        data_path <- tempfile(pattern = "fftRtsne_data_", fileext = ".dat")
    }
    if (is.null(result_path)) {
        result_path <-
            tempfile(pattern = "fftRtsne_result_", fileext = ".dat")
    }

    fast_tsne_path <- normalizePath(path = fast_tsne_path)
    if (!utils::file_test(op = "-x", x = fast_tsne_path)) {
        stop("fast_tsne_path '", fast_tsne_path,
             "' does not exist or is not executable")
    }
    # check fast_tsne version
    ft.out <- suppressWarnings(system2(fast_tsne_path, stdout = TRUE))
    if (grepl(pattern = "= t-SNE v1.1", x = ft.out[1])) {
        version_number <- "1.1.0"
    } else if (grepl(pattern = "= t-SNE v1.0", x = ft.out[1])) {
        version_number <- "1.0"
    } else {
        message("First line of fast_tsne output is")
        message(ft.out[1])
        stop("Our FIt-SNE wrapper requires FIt-SNE v1.X.X, ",
             "please install the appropriate version from ",
             "github.com/KlugerLab/FIt-SNE and have fast_tsne_path ",
             "point to it if it's not in your path")
    }
    is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) {
        return(abs(x = x - round(x = x)) < tol)
    }
    if (version_number == "1.0" && df != 1.0) {
        stop("This version of FIt-SNE does not support df!=1. ",
             "Please install the appropriate version from ",
             "github.com/KlugerLab/FIt-SNE")
    }
    if (!is.numeric(theta) || (theta < 0.0) || (theta > 1.0)) {
        stop("Incorrect theta.")
    }
    if (nrow(X) - 1 < 3 * perplexity) {
        stop("Perplexity is too large.")
    }
    if (!is.matrix(X)) {
        stop("Input X is not a matrix")
    }
    if (!(max_iter > 0)) {
        stop("Incorrect number of iterations.")
    }
    if (!is.wholenumber(x = stop_early_exag_iter) ||
        stop_early_exag_iter < 0) {
        stop("stop_early_exag_iter should be a positive integer")
    }
    if (!is.numeric(x = exaggeration_factor)) {
        stop("exaggeration_factor should be numeric")
    }
    if (!is.wholenumber(x = dims) || dims <= 0) {
        stop("Incorrect dimensionality.")
    }
    if (search_k == -1) {
        if (perplexity > 0) {
            search_k <- n_trees * perplexity * 3
        } else if (perplexity == 0) {
            search_k <- n_trees * max(perplexity_list) * 3
        } else {
            search_k <- n_trees * K * 3
        }
    }
    nbody_algo <- ifelse(test = fft_not_bh, yes = 2, no = 1)
    if (is.null(load_affinities)) {
        load_affinities <- 0
    } else {
        if (load_affinities == "load") {
            load_affinities <- 1
        } else if (load_affinities == "save") {
            load_affinities <- 2
        } else {
            load_affinities <- 0
        }
    }
    knn_algo <- ifelse(test = ann_not_vptree, yes = 1, no = 2)
    f <- file(description = data_path, open = "wb")
    n <- nrow(X)
    D <- ncol(X)
    writeBin(object = as.integer(n), con = f, size = 4)
    writeBin(object = as.integer(D), con = f, size = 4)
    # theta
    writeBin(object = as.numeric(theta), con = f, size = 8)
    # theta
    writeBin(object = as.numeric(perplexity), con = f, size = 8)
    if (perplexity == 0) {
        writeBin(object = as.integer(length(perplexity_list)),
                 con = f, size = 4)
        writeBin(object = perplexity_list, con = f)
    }
    # theta
    writeBin(object = as.integer(dims), con = f, size = 4)
    writeBin(object = as.integer(max_iter), con = f, size = 4)
    writeBin(object = as.integer(stop_early_exag_iter), con = f, size = 4)
    writeBin(object = as.integer(mom_switch_iter), con = f, size = 4)
    writeBin(object = as.numeric(momentum), con = f, size = 8)
    writeBin(object = as.numeric(final_momentum), con = f, size = 8)
    writeBin(object = as.numeric(learning_rate), con = f, size = 8)
    # K
    writeBin(object = as.integer(K), con = f, size = 4)
    # sigma
    writeBin(object = as.numeric(sigma), con = f, size = 8)
    # not barnes hut
    writeBin(object = as.integer(nbody_algo), con = f, size = 4)
    writeBin(object = as.integer(knn_algo), con = f, size = 4)
    # compexag
    writeBin(object = as.numeric(exaggeration_factor), con = f, size = 8)
    writeBin(object = as.integer(no_momentum_during_exag), con = f, size = 4)
    writeBin(object = as.integer(n_trees), con = f, size = 4)
    writeBin(object = as.integer(search_k), con = f, size = 4)
    writeBin(object = as.integer(start_late_exag_iter), con = f, size = 4)
    writeBin(object = as.numeric(late_exag_coeff), con = f, size = 8)
    writeBin(object = as.integer(nterms), con = f, size = 4)
    writeBin(object = as.numeric(intervals_per_integer), con = f, size = 8)
    writeBin(object = as.integer(min_num_intervals), con = f, size = 4)
    tX <- c(t(X))
    writeBin(object = tX, con = f)
    writeBin(object = as.integer(rand_seed), con = f, size = 4)
    if (version_number != "1.0")
        writeBin(object = as.numeric(df), con = f, size = 8)
    writeBin(object = as.integer(load_affinities),
             con = f,
             size = 4)
    if (!is.null(x = initialization))
        writeBin(object = c(t(x = initialization)), con = f)
    close(con = f)
    if (version_number == "1.0")
        flag <- system2(fast_tsne_path, c(data_path, result_path, nthreads))
    else
        flag <- system2(fast_tsne_path, c(version_number, data_path,
                                          result_path, nthreads))
    if (flag != 0)
        stop("tsne call failed")
    f <- file(description = result_path, open = "rb")
    n <- readBin(con = f, what = integer(), n = 1, size = 4)
    d <- readBin(con = f, what = integer(), n = 1, size = 4 )
    Y <- readBin(con = f, what = numeric(), n = n * d)
    Y <- t(x = matrix(Y, nrow = d))
    if (get_costs) {
        tmp <- readBin(con = f, what = integer(), n = 1, size = 4)
        costs <- readBin(con = f, what = numeric(), n = max_iter, size = 8)
        Yout <- list(Y = Y, costs = costs)
    } else {
        Yout <- Y
    }
    close(con = f)
    file.remove(data_path)
    file.remove(result_path)
    return(Yout)
}
