# selectGeneGlobalRank <- function(
#         object,
#         n = 4000,
#         alpha = 0.99,
#         useDatasets = NULL,
#         unsharedDatasets = NULL,
#         chunk = 1000,
#         verbose = getOption("ligerVerbose")
# ) {
#     .checkObjVersion(object)
#     # A bunch of input checks at first ####
#     useDatasets <- .checkUseDatasets(object, useDatasets)
#     object <- recordCommand(object, dependencies = "hdf5r")
#     if (!is.null(unsharedDatasets))
#         unsharedDatasets <- .checkUseDatasets(object, unsharedDatasets)
#     involved <- unique(c(useDatasets, unsharedDatasets))
#
#     shared.features <- Reduce(intersect, lapply(datasets(object)[involved],
#                                                 rownames))
#     perDatasetSelect <- list()
#     for (d in involved) {
#         ld <- dataset(object, d)
#         if (is.null(normData(ld))) {
#             warning("Dataset \"", d, "\" is not normalized, skipped")
#             next
#         }
#         ## Make sure that all required feature meta values exist ####
#         if (isH5Liger(ld)) {
#             ld <- calcGeneVars.H5(ld, chunkSize = chunk,
#                                   verbose = verbose)
#         } else {
#             featureMeta(ld, check = FALSE)$geneMeans <-
#                 Matrix::rowMeans(normData(ld))
#             featureMeta(ld, check = FALSE)$geneVars <-
#                 rowVars_sparse_rcpp(normData(ld), featureMeta(ld)$geneMeans)
#         }
#         datasets(object, check = FALSE)[[d]] <- ld
#         ## The real calculation starts here ####
#         geneMeans <- featureMeta(ld)$geneMeans
#         geneVars <- featureMeta(ld)$geneVars
#         trx_per_cell <- cellMeta(object, "nUMI", cellIdx = object$dataset == d)
#         nolan_constant <- mean((1 / trx_per_cell))
#         alphathresh.corrected <- alpha / nrow(ld)
#         geneMeanUpper <- geneMeans +
#             stats::qnorm(1 - alphathresh.corrected / 2) *
#             sqrt(geneMeans * nolan_constant / ncol(ld))
#         basegenelower <- log10(geneMeans * nolan_constant)
#         pass.upper <- geneVars / nolan_constant > geneMeanUpper
#         pass.lower <- log10(geneVars) > basegenelower
#         preselected <- data.frame(
#             gene = rownames(ld),
#             dataset = d,
#             shared = d %in% useDatasets,
#             unshared = d %in% unsharedDatasets,
#             varianceDiff = geneVars - basegenelower
#         )
#         perDatasetSelect[[d]] <- preselected[pass.upper & pass.lower,]
#     }
#     perDatasetSelect <- Reduce(rbind, perDatasetSelect)
#     # For shared
#     shareTable <- perDatasetSelect[perDatasetSelect$gene %in% shared.features,]
#     rank <- order(shareTable$varianceDiff, decreasing = TRUE)
#     shareTable <- shareTable[rank,]
#     line <- 1
#     while (length(unique(shareTable$gene[seq(line)])) < n) line <- line + 1
#     return(unique(shareTable$gene[seq(line)]))
# }
#
#
# .scaleH5Matrix <- function(ld, featureIdx, resultH5Path, chunk, verbose) {
#     features <- rownames(ld)[featureIdx]
#     geneSumSq <- featureMeta(ld)$geneSumSq[featureIdx]
#     nCells <- ncol(ld)
#     geneRootMeanSumSq = sqrt(geneSumSq / (nCells - 1))
#     h5file <- getH5File(ld)
#     safeH5Create(
#         ld,
#         dataPath = resultH5Path,
#         dims = c(length(features), nCells),
#         dtype = "double",
#         chunkSize = c(length(features), chunk)
#     )
#     H5Apply(
#         ld,
#         useData = "normData",
#         chunkSize = chunk,
#         verbose = verbose,
#         FUN = function(chunk, sparseXIdx, cellIdx, values) {
#             chunk <- chunk[featureIdx, , drop = FALSE]
#             chunk = as.matrix(chunk)
#             chunk = sweep(chunk, 1, geneRootMeanSumSq, "/")
#             rownames(chunk) <- features
#             chunk[is.na(chunk)] = 0
#             chunk[chunk == Inf] = 0
#             h5file[[resultH5Path]][seq_along(features),
#                                    cellIdx] <- chunk
#         }
#     )
#     h5fileInfo(ld, "scaleData", check = FALSE) <- resultH5Path
#     safeH5Create(
#         ld,
#         dataPath = paste0(resultH5Path, ".featureIdx"),
#         dims = length(features),
#         dtype = "int"
#     )
#     h5file[[paste0(resultH5Path, ".featureIdx")]][1:length(featureIdx)] <-
#         featureIdx
#     return(ld)
# }


# `r lifecycle::badge("experimental")` Perform iterative cross-validation to
# suggest the best k for the given data
# @description
# For a series of k values, we perform n-fold cross validation on the given
# data. We use the train data to get factorized low-rank matrices (the \eqn{W}
# and \eqn{V} matrices for gene loading in metagenes) and predict the \eqn{H}
# matrices (for metagene loading in cells) for the test data. The objective
# error (see \code{\link{runINMF}}) is calculated for the test data and shown
# as the mean squared error (MSE) for each k. The k with the lowest MSE is
# suggested as the best k for the given data.
# @noRd

# suggestK <- function(
        #         object,
#         kTest = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
#         nFold = 5,
#         lambda = 5,
#         nIteration = 30,
#         verbose = getOption("ligerVerbose", TRUE)
# ) {
#     if ((nFold %% 1) != 0) {
#         cli::cli_abort("{.field nFold} must be an integer")
#     }
#     if (nFold < 2) {
#         cli::cli_abort("{.field nFold} must be at least 2")
#     }
#     resultDF <- data.frame(
#         # Should be like 5, 5, 5, 10, 10, 10, ..
#         k = rep(kTest, each = nFold),
#         objErr = rep(NA, length(kTest*nFold)),
#         fold = rep(seq(nFold), length(kTest))
#     )
#     on.exit({
#         return(list(result = resultDF, figure = .plotSuggestK(resultDF)))
#     })
#     nCells <- ncol(object)
#     cliID <- cli::cli_progress_bar(total = length(kTest) * nFold)
#     if (isTRUE(verbose)) {
#         cli::cli_alert_info("The progress might take long. Current result will still be returned if interrupted.")
#     }
#     for (k in kTest) {
#         permuteIdx <- sample(nCells)
#         foldAssign <- (seq(nCells) %% nFold) + 1
#         # Check if Each fold has all batchces covered
#         datasetVar <- object$dataset
#         if (table(datasetVar[permuteIdx], foldAssign) %>% `==`(0) %>% any()) {
#             cli::cli_abort(c(
#                 x = "Each fold must have all dataset covered.",
#                 i = "Please consider using smaller {.field nFold}."
#             ))
#         }
#         for (fold in seq(nFold)) {
#             testIdx <- permuteIdx[foldAssign == fold]
#             testIdx <- sort(testIdx)
#             testObj <- subsetLiger(object, cellIdx = testIdx, useSlot = "scaleData", verbose = FALSE)
#             trainIdx <- permuteIdx[foldAssign != fold]
#             trainIdx <- sort(trainIdx)
#             trainObj <- subsetLiger(object, cellIdx = trainIdx, useSlot = "scaleData", verbose = FALSE)
#             cli::cli_progress_update(
#                 inc = 0,
#                 id = cliID,
#                 status = sprintf("k = %d, fold = %d, training", k, fold)
#             )
#             trainObj <- runINMF(trainObj, k = k, lambda = lambda, nIteration = nIteration, verbose = FALSE)
#             cli::cli_progress_update(
#                 inc = 0,
#                 id = cliID,
#                 status = sprintf("k = %d, fold = %d, validating", k, fold)
#             )
#             W <- getMatrix(trainObj, "W")
#             Vs <- getMatrix(trainObj, "V")
#             H_pred <- lapply(
#                 X = seq_along(Vs),
#                 FUN = .solveH_i,
#                 # Params to `solveH_i`
#                 object = testObj,
#                 W = W,
#                 Vs = Vs,
#                 lambda = lambda
#             )
#             objErr <- 0
#             for (i in seq_along(H_pred)) {
#                 objErr <- objErr +
#                     objErr_i(
#                         H = H_pred[[i]],
#                         W = W,
#                         V = Vs[[i]],
#                         E = scaleData(testObj, i),
#                         lambda = 0
#                     )
#             }
#             resultDF[resultDF$k == k & resultDF$fold == fold, "objErr"] <- objErr
#             cli::cli_progress_update(
#                 inc = 1,
#                 id = cliID,
#                 status = sprintf("k = %d, fold = %d, done", k, fold)
#             )
#         }
#     }
# }
#
# .solveH_i <- function(i, object, W, Vs, lambda) {
#     Xi <- scaleData(object, i)
#     Vi <- Vs[[i]]
#     WV <- W + Vi
#     CtC <- t(WV) %*% WV + lambda * t(Vi) %*% Vi
#     CtB <- t(WV) %*% Xi
#     RcppPlanc::bppnnls_prod(CtC, as.matrix(CtB))
# }
#
# .plotSuggestK <- function(stats) {
#     stats.mean <- stats %>%
#         dplyr::group_by(.data[["k"]]) %>%
#         dplyr::summarize(MSE = sqrt(mean(.data[["objErr"]])))
#     stats.error <- stats %>%
#         dplyr::group_by(.data[["k"]]) %>%
#         dplyr::summarize(ymax = sqrt(max(.data[["objErr"]])), ymin = sqrt(min(.data[["objErr"]])))
#     ggplot2::ggplot() +
#         ggplot2::geom_point(
#             data = stats,
#             mapping = ggplot2::aes(.data[["k"]], sqrt(.data[["objErr"]]), group = .data[["fold"]]),
#             position = ggplot2::position_jitter(width = 0.4, height = 0),
#             size = .8
#         ) +
#         ggplot2::geom_errorbar(
#             data = stats.error,
#             mapping = ggplot2::aes(x = .data[["k"]], ymin = .data[["ymin"]], ymax = .data[["ymax"]]),
#             inherit.aes = FALSE,
#             width = 1.5, linewidth = 0.5
#         ) +
#         ggplot2::geom_line(
#             data = stats.mean,
#             mapping = ggplot2::aes(.data[["k"]], .data[["MSE"]]),
#             inherit.aes = FALSE,
#             color = "#54B0E4",
#         ) +
#         ggplot2::labs(
#             x = "k",
#             y = "Squared Error"
#         ) +
#         ggplot2::theme(
#             panel.border = ggplot2::element_rect(fill = NA, color = "black"),
#             panel.background = ggplot2::element_blank(),
#             panel.grid = ggplot2::element_line(color = "grey", linetype = "dashed")
#         )
# }
