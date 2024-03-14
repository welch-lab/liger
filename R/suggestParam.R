#' Visually suggest appropiate k value
#' @description
#' This can be used to select appropriate value of k for factorization of
#' particular dataset. Plots median (across cells in all datasets) K-L
#' divergence from uniform for cell factor loadings as a function of k. This
#' should increase as k increases but is expected to level off above
#' sufficiently high number of factors (k). This is because cells should have
#' factor loadings which are not uniformly distributed when an appropriate
#' number of factors is reached.
#' @param object A \linkS4class{liger} object with scaled data.
#' @param kTest Numeric vector of "k" values to test. Default \code{seq(5, 50,
#' 5)}.
#' @param lambda Lambda to use for all factorization runs. Default \code{5}.
#' @param thresh Convergence threshold. Convergence occurs when
#' \eqn{|obj_0-obj|/(mean(obj_0,obj)) < thresh}. Default \code{1e-4}.
#' @param maxIter Maximum number of block coordinate descent iterations to
#' perform. Default \code{100}.
#' @param nCores Number of cores to use for optimizing factorization in
#' parallel. Default \code{1}.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param genNew Logical. Whether to use \code{\link{optimizeNewK}} in
#' factorization. Results in slower factorization. Default \code{FALSE}.
#' @param nrep See \code{\link{optimizeALS}} for detail. Number restarts to
#' perform at each k value tested. Increase to produce smoother curve if results
#' unclear. Default \code{1}.
#' @param plotLog2 Logical. Whether to plot log2 curve for reference on K-L
#' plot. log2 is upper bound and can sometimes help in identifying "elbow" of
#' plot. Default \code{TRUE}.
#' @param returnRaw Logical. Whether to return raw data or data.frame. Raw data
#' is a list of matrices of K-L divergences (length(k.test) by n_cells). Length
#' of list corresponds to \code{nrep}. (default FALSE)
#' @param return.data \bold{Defuncted}. Will always return data now. Use
#' \code{\link{plotSuggestK}} to show the plot.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param k.test,max.iters,num.cores,rand.seed,gen.new,plot.log2,return.raw
#' \bold{Deprecated}. See Usage section for replacement.
#' @return Matrix of results
#' @export
suggestK <- function(
        object,
        kTest = seq(5, 50, 5),
        lambda = 5,
        thresh = 1e-4,
        maxIter = 100,
        nCores = 1,
        seed = 1,
        genNew = FALSE,
        plotLog2 = TRUE,
        nrep = 1,
        returnRaw = FALSE,
        verbose = getOption("ligerVerbose"),
        # Deprecated coding style
        k.test = kTest,
        max.iters = maxIter,
        num.cores = nCores,
        rand.seed = seed,
        gen.new = genNew,
        plot.log2 = plotLog2,
        return.raw = returnRaw,
        # Defunct
        return.data = TRUE
) {
    .deprecateArgs(list(k.test = "kTest", max.iters = "maxIter",
                        num.cores = "nCores", rand.seed = "seed",
                        gen.new = "genNew", plot.log2 = FALSE,
                        return.raw = "returnRaw"),
                   defunct = "return.data")
    if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("Package \"doParallel\" needed for this function to work. ",
             "Please install it by command:\n",
             "install.packages('doParallel')",
             call. = FALSE)
    }
    if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Package \"parallel\" needed for this function to work. ",
             "Please install it by command:\n",
             "install.packages('parallel')",
             call. = FALSE)
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
        stop("Package \"foreach\" needed for this function to work. ",
             "Please install it by command:\n",
             "install.packages('foreach')",
             call. = FALSE)
    }
    `%dopar%` <- foreach::`%dopar%`
    # namespaceImportFrom(ns = "foreach", self = rlang::current_env(),
    #                     vars = "%dopar%")
    scaleData <- getMatrix(object, "scaleData")
    notScaled <- sapply(scaleData, is.null)
    if (any(notScaled))
        stop("Dataset not scaled (",
             paste(names(object)[notScaled], collapse = ", "),
             "), please run `scaleNotCenter` in advance")

    time_start <- Sys.time()
    if (isTRUE(verbose))
        .log("Testing k values. This may take several minutes depending on ",
             "the number of values being tested...")

    # optimize largest k value first to take advantage of efficient updating
    rep_data <- list()
    kTest <- sort(kTest, decreasing = TRUE)
    for (r in seq(nrep)) {
        if (isTRUE(verbose))
            .log("Preprocessing for rep ", r,
                 ": optimizing initial factorization with largest test k=",
                 kTest[1])
        object <- optimizeALS(
            object,
            k = kTest[1],
            lambda = lambda,
            thresh = thresh,
            maxIter = maxIter,
            nrep = 1,
            seed = (seed + r - 1),
            verbose = verbose
        )

        if (isTRUE(verbose)) {
            .log("Testing different choices of k")
        }
        cl <- parallel::makeCluster(nCores)
        doParallel::registerDoParallel(cl)
        #pb <- txtProgressBar(min = 0, max = length(kTest), style = 3,
        #                     initial = 1, file = "")
        # define progress bar function
        #progress <- function(n) setTxtProgressBar(pb, n)
        #opts <- list(progress = progress)
        i <- 0
        data_matrix <- foreach::foreach(
            i = seq_along(kTest),
            .combine = "rbind",
            .packages = "rliger2"
        ) %dopar% {
            if (i != length(kTest)) {
                if (isTRUE(genNew)) {
                    ob.test <- optimizeALS(
                        object,
                        k = kTest[i],
                        lambda = lambda,
                        thresh = thresh,
                        maxIter = maxIter,
                        seed = (seed + r - 1),
                        verbose = verbose
                    )
                } else {
                    ob.test <- optimizeNewK(
                        object,
                        k.new = kTest[i],
                        lambda = lambda,
                        thresh = thresh,
                        max.iters = maxIter,
                        rand.seed = (seed + r - 1),
                        verbose = verbose
                    )
                }
            } else {
                ob.test <- object
            }
            dataset_split <- kl_divergence_uniform(ob.test)
            unlist(dataset_split)
        }
        parallel::stopCluster(cl)
        #data_matrix <- data_matrix[1:nrow(data_matrix), ]
        rep_data[[r]] <- data_matrix
    }

    medians <- Reduce(cbind, lapply(rep_data, function(x) {
        apply(x, 1, stats::median)
    }))
    if (is.null(dim(medians))) {
        medians <- matrix(medians, ncol = 1)
    }
    mean_kls <- apply(medians, 1, mean)

    time_elapsed <- difftime(Sys.time(), time_start, units = "auto")

    if (isTRUE(verbose))
        .log("Completed in: ", as.double(time_elapsed), units(time_elapsed))

    # make dataframe
    df_kl <- data.frame(
        median_kl = c(mean_kls, log2(kTest)),
        k = c(kTest, kTest),
        calc = factor(c(rep("KL_div", length(kTest)),
                        rep("log2(k)", length(kTest))))
    )


    if (isFALSE(plotLog2)) df_kl <- df_kl[df_kl$calc == "KL_div", ]

    if (isTRUE(returnRaw)) {
        rep_data <- lapply(rep_data, function(x) {
            rownames(x) <- kTest
            x
        })
        return(rep_data)
    }
    return(df_kl)
}

#' Visualize suggestK result statistics
#' @param result Returned data.frame from \code{\link{suggestK}}.
#' @param xlab,ylab The title on x-/y-axis of the plot. See Usage for default.
#' @param legendColorTitle The title of legend. Default \code{""}.
#' @param ... Additional plot theme setting passed to
#' \code{\link{.ggplotLigerTheme}}.
#' @return ggplot
plotSuggestK <- function(
        result,
        xlab = "k",
        ylab = "Median KL divergence (across all cells)",
        legendColorTitle = "",
        ...) {
    p <- ggplot2::ggplot(
        result,
        ggplot2::aes(x = .data[["k"]],
                     y = .data[["median_kl"]],
                     color = .data[["calc"]],
                     label = round(.data[["median_kl"]], 2))) +
        ggplot2::geom_line(size = 1) +
        ggplot2::geom_point() +
        ggplot2::geom_text(color = "black", nudge_y = 0.1) +
        ggplot2::scale_x_continuous(breaks = unique(result$k))
    .ggplotLigerTheme(p, xlab = xlab, ylab = ylab,
                      legendColorTitle = legendColorTitle, ...)
}

# helper function for calculating KL divergence from uniform distribution
# (related to Shannon entropy) for factorization
kl_divergence_uniform <- function(object, Hs = NULL) {
    # Input Hs should be k x c
    if (is.null(Hs)) {
        Hs <- getMatrix(object, "H", returnList = TRUE)
    }

    Hs <- lapply(Hs, t)
    n_factors <- ncol(Hs[[1]])
    dataset_list <- list()
    for (i in seq_along(object)) {
        scaled = scale(Hs[[i]], center = FALSE, scale = TRUE)

        inflated = t(apply(scaled, 1, function(x) {
            replace(x, x == 0, 1e-20)
        }))
        inflated = inflated / rowSums(inflated)
        divs = apply(inflated, 1, function(x) {
            log2(n_factors) + sum(log2(x) * x)
        })
        dataset_list[[i]] <- divs
    }
    return(dataset_list)
}

#' Visually suggest appropriate lambda value
#'
#' Can be used to select appropriate value of lambda for factorization of
#' particular dataset. Plot alignment and agreement for various test values of
#' lambda. Most appropriate lambda is likely around the "elbow" of the alignment
#' plot (when alignment stops increasing). This will likely also correspond to
#' slower decrease in agreement. Depending on number of cores used, this process
#' can take 10-20 minutes.
#'
#'
#' @param object A \linkS4class{liger} object with scaled data.
#' @param k Number of factors to use in the test.
#' @param lambdaTest Numeric vector of "lambda" values to test. Default
#' \code{NULL} will test the following sequence: 0.25, 0.5, 0.75; 1, 2 to 10;
#' and 15, 20, ... to 60, totally 23 numbers.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param nCores Number of cores to use for optimizing factorization in
#' parallel. Default \code{1}.
#' @param thresh Convergence threshold. Convergence occurs when
#' \eqn{|obj_0-obj|/(mean(obj_0,obj)) < thresh}. Default \code{1e-4}.
#' @param maxIter Maximum number of block coordinate descent iterations to
#' perform. Default \code{100}.
#'
#' @param nNeighbors Number of nearest neighbors for within-dataset KNN in
#' \code{\link{quantileNorm}}. Default \code{20}.
#' @param reference Reference dataset for \code{\link{quantileNorm}}. See link
#' for detail. Default \code{NULL}.
#'
#' @param genNew Logical. Whether to use \code{\link{optimizeNewLambda}} in
#' factorization. Recommended \code{TRUE} when looking at only a small range of
#' lambdas (i.e. 1:7). Default \code{FALSE}.
#' @param nrep See \code{\link{optimizeALS}} for detail. Number restarts to
#' perform at each k value tested. Increase to produce smoother curve if results
#' unclear. Default \code{1}.
#' @param returnRaw Logical. Whether to return raw data or data.frame. Raw data
#' is matrix of alignment values for each lambda value tested (each column
#' represents a different rep for nrep). Default \code{FALSE}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param lambda.test,rand.seed,num.cores,max.iters,knn_k,ref_dataset,gen.new,return.raw
#' \bold{Deprecated}. See Usage section for replacement.
#' @param return.data \bold{Defuncted}. Will always return data now. Use
#' \code{\link{plotSuggestLambda}} to show the plot.
#' @param k2,resolution \bold{Defuncted}. Parameter passed to
#' \code{\link{quantileAlignSNF}} which is deprecated.
#' @return Matrix of results
#' @export
suggestLambda <- function(
        object,
        k,
        lambdaTest = NULL,
        seed = 1,
        nCores = 1,
        thresh = 1e-4,
        maxIter = 100,
        nNeighbors = 20,
        reference = NULL,
        genNew = FALSE,
        nrep = 1,
        returnRaw = FALSE,
        verbose = getOption("ligerVerbose"),
        # Deprecated coding style
        lambda.test = lambdaTest,
        rand.seed = seed,
        num.cores = nCores,
        max.iters = maxIter,
        knn_k = nNeighbors,
        ref_dataset = reference,
        gen.new = genNew,
        return.raw = returnRaw,
        # Defunct
        return.data = TRUE,
        k2 = 500,
        resolution = 1
) {
    .deprecateArgs(list(lambda.test = "lambdaTest", rand.seed = "seed",
                        num.cores = "nCores", max.iters = "maxIter",
                        knn_k = "nNeighbors", ref_dataset = "reference",
                        gen.new = "genNew", return.raw = "returnRaw"),
                   defunct = "return.data")
    if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("Package \"doParallel\" needed for this function to work. ",
             "Please install it by command:\n",
             "install.packages('doParallel')",
             call. = FALSE)
    }
    if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Package \"parallel\" needed for this function to work. ",
             "Please install it by command:\n",
             "install.packages('parallel')",
             call. = FALSE)
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
        stop("Package \"foreach\" needed for this function to work. ",
             "Please install it by command:\n",
             "install.packages('foreach')",
             call. = FALSE)
    }
    `%dopar%` <- foreach::`%dopar%`
    if (is.null(lambdaTest)) {
        lambdaTest <- c(seq(0.25, 1, 0.25), seq(2, 10, 1), seq(15, 60, 5))
    }
    lambdaTest <- sort(lambdaTest)
    time_start <- Sys.time()
    # optimize smallest lambda value first to take advantage of efficient
    # updating
    if (isTRUE(verbose))
        .log("Testing lambda values. This may take several minutes depending ",
             "on the number of values being tested...")

    rep_data <- list()
    for (r in seq(nrep)) {
        if (isTRUE(verbose))
            .log("Preprocessing for rep ", r,
                 ": optimizing initial factorization with smallest test ",
                 "lambda=", lambdaTest[1])
        object <- optimizeALS(
            object,
            k = k,
            thresh = thresh,
            lambda = lambdaTest[1],
            maxIter = maxIter,
            nrep = 1,
            seed = (seed + r - 1),
            verbose = verbose
        )
        if (isTRUE(verbose)) {
            .log("Testing different choices of lambda")
        }
        cl <- parallel::makeCluster(nCores)
        doParallel::registerDoParallel(cl)
        i <- 0
        data_matrix <- foreach::foreach(
            i = seq_along(lambdaTest),
            .combine = "rbind",
            .packages = "rliger2"
        ) %dopar% {
            if (i != 1) {
                if (isTRUE(genNew)) {
                    ob.test <- optimizeALS(
                        object,
                        k = k,
                        lambda = lambdaTest[i],
                        thresh = thresh,
                        maxIter = maxIter,
                        seed = (seed + r - 1)
                    )
                } else {
                    ob.test <- optimizeNewLambda(
                        object,
                        new.lambda = lambdaTest[i],
                        thresh = thresh,
                        max.iters = maxIter,
                        rand.seed = (seed + r - 1)
                    )
                }
            } else {
                ob.test <- object
            }
            ob.test <- quantileNorm(
                object = object, reference = reference, nNeighbors = nNeighbors
            )
            #ob.test <-
            #    quantileAlignSNF(
            #        ob.test,
            #        knn_k = nNeighbors,
            #        k2 = k2,
            #        resolution = resolution,
            #        ref_dataset = reference,
            #        id.number = i
            #    )
            calcAlignment(ob.test)
        }

        #close(pb)
        parallel::stopCluster(cl)
        rep_data[[r]] <- data_matrix
    }

    aligns <- Reduce(cbind, rep_data)
    if (is.null(dim(aligns))) {
        aligns <- matrix(aligns, ncol = 1)
    }
    mean_aligns <- apply(aligns, 1, mean)

    time_elapsed <- difftime(Sys.time(), time_start, units = "auto")
    if (isTRUE(verbose))
        .log("Completed in: ", as.double(time_elapsed), units(time_elapsed))

    # make dataframe
    df_al <- data.frame(align = mean_aligns, lambda = lambdaTest)

    if (returnRaw) {
        rownames(aligns) <- lambdaTest
        return(aligns)
    }
    return(df_al)
}

#' Visualize suggestLambda result statistics
#' @param result Returned data.frame from \code{\link{suggestLambda}}.
#' @param xlab,ylab The title on x-/y-axis of the plot. See Usage for default.
#' @param legendColorTitle The title of legend. Default \code{""}.
#' @param ... Additional plot theme setting passed to
#' \code{\link{.ggplotLigerTheme}}.
#' @return ggplot
plotSuggestLambda <- function(
        result,
        xlab = "Lambda",
        ylab = "Alignment",
        legendColorTitle = "",
        ...
) {
    p <- ggplot2::ggplot(
        result,
        ggplot2::aes(x = .data[['lambda']],
                     y = .data[['mean_aligns']],
                     label = round(.data[["mean_aligns"]], 2))) +
        ggplot2::geom_line(size = 1) +
        ggplot2::geom_point() +
        ggplot2::geom_text(color = "black", nudge_y = 0.1)
    .ggplotLigerTheme(p, xlab = xlab, ylab = ylab,
                      legendColorTitle = legendColorTitle, ...)
}
