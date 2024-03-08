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
