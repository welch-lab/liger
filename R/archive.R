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

#' #' Perform Wilcoxon rank-sum test
#' #' @description Perform Wilcoxon rank-sum tests on specified dataset using
#' #' given method.
#' #' @param object A \linkS4class{liger} object with cluster assignment available.
#' #' @param useDatasets A character vector of the names, a numeric or logical
#' #' vector of the index of the datasets to be normalized. Default
#' #' \code{NULL} uses all datasets.
#' #' @param method Choose from \code{"clusters"} or \code{"datasets"}. Default
#' #' \code{"clusters"} compares between clusters across all datasets, while
#' #' \code{"datasets"} compares between datasets within each cluster.
#' #' @param useCluster The name of the column in \code{cellMeta} slot storing the
#' #' cluster assignment variable. Default \code{"leiden_cluster"}
#' #' @param usePeak Logical, whether to test peak counts instead of gene
#' #' expression. Requires presence of ATAC modility datasets. Default
#' #' \code{FALSE}.
#' #' @param verbose Logical. Whether to show information of the progress. Default
#' #' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' #' @param data.use,compare.method \bold{Deprecated}. See Usage section for
#' #' replacement.
#' #' @return A 10-columns data.frame with test results.
#' #' @export
#' #' @examples
#' #' library(dplyr)
#' #' result <- runWilcoxon(pbmcPlot)
#' #' result %>% group_by(group) %>% top_n(2, logFC)
#' runWilcoxon <- function(
#'         object,
#'         useDatasets = NULL,
#'         method = c("clusters", "datasets"),
#'         useCluster = "leiden_cluster",
#'         usePeak = FALSE,
#'         verbose = getOption("ligerVerbose"),
#'         # Deprecated coding style,
#'         data.use = useDatasets,
#'         compare.method = method
#' ) {
#'     .deprecateArgs(list(data.use = "useDatasets", compare.method = "method"))
#'     method <- match.arg(method)
#'     # Input checks ####
#'     useDatasets <- .checkUseDatasets(object, useDatasets,
#'                                      modal = ifelse(usePeak, "atac", "default"))
#'     if (!isTRUE(usePeak)) {
#'         if (method == "datasets" & length(useDatasets) < 2)
#'             stop("Should have at least 2 datasets as input ",
#'                  "when compare between datasets")
#'         if (isH5Liger(object, useDatasets)) {
#'             stop("HDF5 based datasets detected but is not supported. \n",
#'                  "Try `object.sub <- downsample(object, useSlot = ",
#'                  "'normData')` to create ANOTHER object with in memory data.")
#'         }
#'         allNormed <- all(sapply(datasets(object),
#'                                 function(ld) !is.null(normData(ld))))
#'         if (!allNormed)
#'             stop("All datasets being involved has to be normalized")
#'
#'         ## get all shared genes of datasets involved
#'         normDataList <- getMatrix(object, "normData", dataset = useDatasets,
#'                                   returnList = TRUE)
#'         features <- Reduce(intersect, lapply(normDataList, rownames))
#'         normDataList <- lapply(normDataList, function(x) x[features,])
#'         featureMatrix <- Reduce(cbind, normDataList)
#'     } else {
#'         if (method == "datasets" || length(useDatasets) != 1)
#'             stop("For wilcoxon test on peak counts, can only use ",
#'                  "\"cluster\" method on one dataset.")
#'         normPeakList <- lapply(useDatasets, function(d) normPeak(object, d))
#'         features <- Reduce(intersect, lapply(normPeakList, rownames))
#'         featureMatrix <- Reduce(cbind, normPeakList)
#'         #featureMatrix <- normPeak(object, useDatasets)
#'         if (is.null(featureMatrix))
#'             stop("Peak counts of specified dataset has to be normalized. ",
#'                  "Please try `normalizePeak(object, useDatasets = '",
#'                  useDatasets, "')`.")
#'         #features <- rownames(featureMatrix)
#'     }
#'
#'     ## Subset metadata involved
#'     cellIdx <- object$dataset %in% useDatasets
#'     cellSource <- object$dataset[cellIdx]
#'     clusters <- .fetchCellMetaVar(object, useCluster, cellIdx = cellIdx,
#'                                   checkCategorical = TRUE)
#'
#'     if (isTRUE(verbose))
#'         .log("Performing Wilcoxon test on ", length(useDatasets), " datasets: ",
#'              paste(useDatasets, collapse = ", "))
#'     # perform wilcoxon test ####
#'     if (method == "clusters") {
#'         # compare between clusters across datasets
#'         nfeatures <- length(features)
#'         if (nfeatures > 100000) {
#'             if (isTRUE(verbose)) .log("Calculating Large-scale Input...")
#'             results <- Reduce(rbind, lapply(
#'                 suppressWarnings(split(seq(nfeatures),
#'                                        seq(nfeatures / 100000))),
#'                 function(index) {
#'                     fm <- log1p(1e10*featureMatrix[index, ])
#'                     wilcoxauc(fm, clusters)
#'                 }))
#'         } else {
#'             # TODO: If we add log-transformation to normalization method in the
#'             # future, remember to have conditions here.
#'             fm <- log1p(1e10*featureMatrix)
#'             results <- wilcoxauc(fm, clusters)
#'         }
#'     } else {
#'         # compare between datasets within each cluster
#'         results <- Reduce(rbind, lapply(levels(clusters), function(cluster) {
#'             clusterIdx <- clusters == cluster
#'             subLabel <- paste0(cluster, "-", cellSource[clusterIdx])
#'             if (length(unique(subLabel)) == 1) {
#'                 # if cluster has only 1 data source
#'                 warning("Skipped Cluster ", cluster,
#'                         " since it has only one dataset source.")
#'                 return()
#'             } else {
#'                 subMatrix <- log1p(1e10*featureMatrix[, clusterIdx])
#'                 return(wilcoxauc(subMatrix, subLabel))
#'             }
#'         }))
#'     }
#'     return(results)
#' }
