#' Perform UMAP dimensionality reduction
#' @description
#' Run UMAP on the normalized cell factors (or raw cell factors) to generate a 2D embedding for
#' visualization (or general dimensionality reduction). Has option to run on subset of factors.
#' Note that running multiple times will overwrite tsne.coords values. It is generally
#' recommended to use this method for dimensionality reduction with extremely large datasets.
#'
#' Note that this method requires that the package uwot is installed. It does not depend
#' on reticulate or python umap-learn.
#' @param object \code{liger} object. Should run quantile_norm before calling with defaults.
#' @param use.raw Whether to use un-aligned cell factor loadings (H matrices) (default FALSE).
#' @param dims.use Factors to use for computing tSNE embedding (default 1:ncol(H.norm)).
#' @param k Number of dimensions to reduce to (default 2).
#' @param distance Mtric used to measure distance in the input space. A wide variety of metrics are
#'   already coded, and a user defined function can be passed as long as it has been JITd by numba.
#'   (default "euclidean", alternatives: "cosine", "manhattan", "hamming")
#' @param n_neighbors Number of neighboring points used in local approximations of manifold
#'   structure. Larger values will result in more global structure being preserved at the loss of
#'   detailed local structure. In general this parameter should often be in the range 5 to 50, with
#'   a choice of 10 to 15 being a sensible default. (default 10)
#' @param min_dist Controls how tightly the embedding is allowed compress points together. Larger
#'   values ensure embedded points are more evenly distributed, while smaller values allow the
#'   algorithm to optimise more accurately with regard to local structure. Sensible values are in
#'   the range 0.001 to 0.5, with 0.1 being a reasonable default. (default 0.1)
#' @param rand.seed Random seed for reproducibility (default 42).
#' @return \code{liger} object with tsne.coords slot set.
#' @export
runUMAP <- function(object,
                    use.raw = FALSE,
                    dims.use = NULL,
                    k = 2,
                    distance = "euclidean",
                    n_neighbors = 10,
                    min_dist = 0.1,
                    rand.seed = 42) {
    set.seed(rand.seed)
    if (isTRUE(use.raw)) {
        data.use <- t(Reduce(cbind, getMatrix(object, "H")))
    } else {
        data.use <- getMatrix(object, "H.norm")
    }
    if (is.null(dims.use)) dims.use <- seq(ncol(data.use))
    umap <- uwot::umap(data.use[, dims.use],
                       n_components = as.integer(k),
                       metric = distance,
                       n_neighbors = as.integer(n_neighbors),
                       min_dist = min_dist)
    rownames(umap) <- colnames(object)
    colnames(umap) <- paste0("UMAP_", seq(k))
    cell.meta(object) <- cbind(cell.meta(object), umap)
    object@uns$UMAP <- list(variableName = colnames(umap))
    return(object)
}

#' Perform t-SNE dimensionality reduction
#'
#' Runs t-SNE on the normalized cell factors (or raw cell factors) to generate a 2D embedding for
#' visualization. Has option to run on subset of factors. Note that running multiple times will
#' reset tsne.coords values.
#'
#' In order to run fftRtsne (recommended for large datasets), you must first install FIt-SNE as
#' detailed \href{https://github.com/KlugerLab/FIt-SNE}{here}. Include the path to the cloned
#' FIt-SNE directory as the fitsne.path parameter, though this is only necessary for the first call
#' to runTSNE. For more detailed FIt-SNE installation instructions, see the liger repo README.
#'
#' @param object \code{liger} object. Should run quantile_norm before calling with defaults.
#' @param use.raw Whether to use un-aligned cell factor loadings (H matrices) (default FALSE).
#' @param dims.use Factors to use for computing tSNE embedding (default 1:ncol(H.norm)).
#' @param use.pca Whether to perform initial PCA step for Rtsne (default FALSE).
#' @param perplexity Parameter to pass to Rtsne (expected number of neighbors) (default 30).
#' @param theta Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE
#'   (default 0.5).
#' @param method Supports two methods for estimating tSNE values: Rtsne (Barnes-Hut implementation
#'   of t-SNE) and fftRtsne (FFT-accelerated Interpolation-based t-SNE) (using Kluger Lab
#'   implementation). (default Rtsne)
#' @param fitsne.path Path to the cloned FIt-SNE directory (ie. '/path/to/dir/FIt-SNE') (required
#'   for using fftRtsne -- only first time runTSNE is called) (default NULL).
#' @param rand.seed Random seed for reproducibility (default 42).
#' @return \code{liger} object with tsne.coords slot set.
#' @export
runTSNE <- function(object,
                    use.raw = FALSE,
                    dims.use = NULL,
                    use.pca = FALSE,
                    perplexity = 30,
                    theta = 0.5,
                    method = c("Rtsne", "fftRtsne"),
                    fitsne.path = NULL,
                    rand.seed = 42) {
    method <- match.arg(method)
    if (isTRUE(use.raw)) {
        data.use <- t(Reduce(cbind, getMatrix(object, "H")))
    } else {
        data.use <- getMatrix(object, "H.norm")
    }
    if (is.null(dims.use)) dims.use <- seq(ncol(data.use))
    if (method == "Rtsne") {
        set.seed(rand.seed)
        tsne <- Rtsne::Rtsne(data.use[, dims.use],
                             pca = use.pca,
                             check_duplicates = FALSE,
                             theta = theta,
                             perplexity = perplexity)
        tsne <- tsne$Y
    } else if (method == "fftRtsne") {
        tsne <- .fftRtsne(X = data.use[, dims.use],
            rand_seed = rand.seed,
            fast_tsne_path = fitsne.path,
            theta = theta,
            perplexity = perplexity)
    }
    rownames(tsne) <- colnames(object)
    colnames(tsne) <- paste0("TSNE_", seq(ncol(tsne)))
    cell.meta(object) <- cbind(cell.meta(object), tsne)
    object@uns$TSNE <- list(variableName = colnames(tsne), method = method)
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
    }
    else {
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
