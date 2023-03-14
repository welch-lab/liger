#' Create liger object
#' @description This function allows creating \linkS4class{liger} object from
#' multiple datasets of various forms (See \code{rawData}).
#' @param rawData Named list of datasets. Required. Elements allowed include a
#' matrix, a \code{Seurat} object, a \code{SingleCellExperiment} object, an
#' \code{AnnData} object, a \linkS4class{ligerDataset} object or a filename to
#' an HDF5 file. See detail for HDF5 reading.
#' @param modal Character vector for modality setting. Currently options of
#' \code{"default"}, \code{"rna"}, and \code{"atac"} are supported.
#' @param cellMeta data.frame of metadata at single-cell level. Default
#' \code{NULL}.
#' @param removeMissing Logical. Whether to remove cells that do not have any
#' counts and features not expressed in any cells from each dataset. Default
#' \code{TRUE}. H5 based dataset with less than 8000 cells will be subset into
#' memory.
#' @param formatType Select preset of H5 file structure. Current available
#' options are \code{"10X"} and \code{"AnnData"}. Can be either a single
#' specification for all datasets or a character vector that match with each
#' dataset.
#' @param dataName,indicesName,indptrName The path in a H5 file for the raw
#' sparse matrix data. These three types of data stands for the \code{x},
#' \code{i}, and \code{p} slots of a \code{\link[Matrix]{dgCMatrix-class}}
#' object. Default \code{NULL} uses \code{formatType} preset.
#' @param genesName,barcodesName The path in a H5 file for the gene names and
#' cell barcodes. Default \code{NULL} uses \code{formatType} preset.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Additional slot values that should be directly placed in object.
#' @param remove.missing,format.type,data.name,indices.name,indptr.name,genes.name,barcodes.name
#' \bold{Deprecated.} See Usage section for replacement.
#' @export
#' @seealso \code{\link{createLigerDataset}}, \code{\link{createH5LigerDataset}}
createLiger <- function(
        rawData,
        modal = NULL,
        cellMeta = NULL,
        removeMissing = TRUE,
        formatType = "10X",
        dataName = NULL,
        indicesName = NULL,
        indptrName = NULL,
        genesName = NULL,
        barcodesName = NULL,
        verbose = getOption("ligerVerbose"),
        ...,
        # Deprecated coding style
        remove.missing = removeMissing,
        format.type = formatType,
        data.name = dataName,
        indices.name = indicesName,
        indptr.name = indptrName,
        genes.name = genesName,
        barcodes.name = barcodesName
) {
    .deprecateArgs(list(remove.missing = "removeMissing",
                        format.type = "formatType", data.name = "dataName",
                        indices.name = "indicesName",
                        indptr.name = "indptrName", genes.name = "genesName",
                        barcodes.name = "barcodesName"))
    if (!is.list(rawData)) stop("`rawData` has to be a named list.")

    nData <- length(rawData)
    if (missing(modal) || is.null(modal)) modal <- "default"
    modal <- tolower(modal)
    if (length(modal) == 1) modal <- rep(modal, nData)
    else if (length(modal) != nData)
        stop("Wrong length of `modal`. ",
             "Specify only 1 or match the length of `datasets`. ",
             "See ?createLiger for valid options.")
    # TODO handle h5 specific argument for hybrid of H5 and in memory stuff.
    datasets <- list()
    barcodesOrig <- NULL
    cellID <- NULL
    for (i in seq_along(rawData)) {
        dname <- names(rawData)[i]
        data <- rawData[[i]]
        if (is.character(data)) {
            # Assuming character input is a filename
            datasets[[dname]] <- createH5LigerDataset(
                h5file = data,
                formatType = formatType,
                rawData = dataName,
                barcodesName = barcodesName,
                genesName = genesName,
                indicesName = indicesName,
                indptrName = indptrName,
                modal = modal[i]
            )
            barcodesOrig <- c(barcodesOrig, colnames(datasets[[dname]]))
            cellID <- paste0(dname, "_", colnames(datasets[[dname]]))
            colnames(datasets[[dname]]) <- cellID
        } else if (inherits(data, c("matrixLike"))) {
            datasets[[dname]] <- as.ligerDataset(data, modal = modal[i])
            barcodesOrig <- c(barcodesOrig, colnames(datasets[[dname]]))
            cellID <- paste0(dname, "_", colnames(datasets[[dname]]))
            colnames(datasets[[dname]]) <- cellID
        } else {
            datasets[[dname]] <- as.ligerDataset(data, modal = modal[i])
            barcodesOrig <- c(barcodesOrig, colnames(datasets[[dname]]))
            cellID <- colnames(datasets[[dname]])
        }
    }

    #barcodesOrig <- unlist(lapply(datasets, colnames), use.names = FALSE)
    datasets <- .dedupLigerDatasets(datasets)
    cellID <- unlist(lapply(datasets, colnames), use.names = FALSE)
    if (is.null(cellMeta)) {
        cellMeta <- S4Vectors::DataFrame(
            dataset = factor(rep(names(datasets), lapply(datasets, ncol)),
                             levels = names(datasets)),
            barcode = barcodesOrig,
            row.names = cellID)
    } else {
        cellMeta <- S4Vectors::DataFrame(cellMeta)
        # TODO Whether the given cell metadata dataframe should have original
        # barcodes or processed cellID?
        cellMeta <- cellMeta[barcodesOrig, , drop = FALSE]
        rownames(cellMeta) <- cellID
        # Force writing `dataset` variable as named by @datasets
        cellMeta$dataset <- factor(rep(names(datasets),
                                       lapply(datasets, ncol)))
    }
    obj <- methods::new("liger",
                        datasets = datasets,
                        cellMeta = cellMeta,
                        ...)
    obj <- runGeneralQC(obj, verbose = verbose)
    if (isTRUE(removeMissing)) {
        obj <- removeMissing(obj, "both", filenameSuffix = "qc",
                             verbose = verbose)
    }

    return(obj)
}

#' Deduplicate barcodes from all datasets
#' @param  datasets list of ligerDataset object
#' @noRd
.dedupLigerDatasets <- function(datasets) {
    #for (d in names(datasets)) {
    #    colnames(datasets[[d]]) <- paste0(d, "_", colnames(datasets[[d]]))
    #}

    allBarcodes <- unlist(lapply(datasets, colnames), use.names = FALSE)
    if (any(duplicated(allBarcodes))) {
        warning("Duplicated barcodes detected within a single dataset.")
        datasetVar <- rep(names(datasets), lapply(datasets, ncol))
        # Main deduplication process. Wondering if make stand-alone dot-func
        dups <- unique(allBarcodes[duplicated(allBarcodes)])
        for (bc in dups) {
            idx <- which(allBarcodes == bc)
            allBarcodes[idx] <- paste0(bc, "-", seq_along(idx))
        }
        # Done. Then assign them back to ligerDataset objects
        for (d in names(datasets)) {
            colnames(datasets[[d]]) <- allBarcodes[datasetVar == d]
        }
    }
    return(datasets)
}

#' Create in-memory ligerDataset object
#' @param rawData,normData,scaleData A \code{\link[Matrix]{dgCMatrix-class}}
#' object for the raw or normalized expression count or a dense matrix of scaled
#' variable gene expression, respectively. Default \code{NULL} for all three but
#' at lease one has to be specified.
#' @param modal Name of modality for this dataset. Currently options of
#' \code{"default"}, \code{"rna"}, and \code{"atac"} are supported. Default
#' \code{"default"}.
#' @param featureMeta Data frame of feature metadata. Default \code{NULL}.
#' @param ... Additional slot data. See \linkS4class{ligerDataset} for detail.
#' Given values will be directly placed at corresponding slots.
#' @seealso \linkS4class{ligerDataset}, \linkS4class{ligerATACDataset},
#' \code{\link{createH5LigerDataset}}
#' @export
createLigerDataset <- function(
        rawData = NULL,
        modal = c("default", "rna", "atac"),
        normData = NULL,
        scaleData = NULL,
        featureMeta = NULL,
        ...
) {
    modal <- match.arg(modal)
    args <- as.list(environment())
    additional <- list(...)
    # Necessary initialization of slots
    if (is.null(rawData) && is.null(normData) && is.null(scaleData)) {
        stop("At least one type of expression data (rawData, normData or ",
             "scaleData) has to be provided")
    }
    # Look for proper colnames and rownames
    cn <- NULL
    rn <- NULL
    for (i in c("rawData", "normData", "scaleData")) {
        cn <- colnames(args[[i]])
        if (!is.null(cn)) break
    }
    if (!is.null(rawData)) {
        rn <- rownames(rawData)
        if (!inherits(rawData, "dgCMatrix"))
            rawData <- methods::as(rawData, "CsparseMatrix")
    }
    if (!is.null(normData)) {
        if (is.null(rn)) rn <- rownames(normData)
        if (!inherits(normData, "dgCMatrix"))
            rawData <- methods::as(normData, "CsparseMatrix")
    }
    if (!is.null(scaleData)) {
        if (is.null(rn)) rn <- rownames(scaleData)
    }
    if (is.null(h5fileInfo)) h5fileInfo <- list()
    if (is.null(featureMeta))
        featureMeta <- S4Vectors::DataFrame(row.names = rn)
    else if (!inherits(featureMeta, "DFrame"))
        featureMeta <- S4Vectors::DataFrame(featureMeta)
    # Create ligerDataset
    allData <- list(.modalClassDict[[modal]],
                    rawData = rawData, normData = normData,
                    scaleData = scaleData, featureMeta = featureMeta,
                    colnames = cn, rownames = rn)
    allData <- c(allData, additional)

    x <- do.call("new", allData)
    return(x)
}

#' Create on-disk ligerDataset Object
#' @param h5file Filename of a H5 file
#' @param formatType Select preset of H5 file structure. Current available
#' options are \code{"10X"} and \code{"AnnData"}.
#' @param rawData,indicesName,indptrName The path in a H5 file for the raw
#' sparse matrix data. These three types of data stands for the \code{x},
#' \code{i}, and \code{p} slots of a \code{\link[Matrix]{dgCMatrix-class}}
#' object. Default \code{NULL} uses \code{formatType} preset.
#' @param normData The path in a H5 file for the "x" vector of the normalized
#' sparse matrix. Default \code{NULL}.
#' @param scaleData The path in a H5 file for the dense 2D scaled matrix.
#' Default \code{NULL}.
#' @param genesName,barcodesName The path in a H5 file for the gene names and
#' cell barcodes. Default \code{NULL} uses \code{formatType} preset.
#' @param modal Name of modality for this dataset. Currently options of
#' \code{"default"}, \code{"rna"}, and \code{"atac"} are supported. Default
#' \code{"default"}.
#' @param featureMeta Data frame for feature metadata. Default \code{NULL}.
#' @param ... Additional slot data. See \linkS4class{ligerDataset} for detail.
#' Given values will be directly placed at corresponding slots.
createH5LigerDataset <- function(
        h5file,
        formatType = NULL,
        rawData = NULL,
        normData = NULL,
        scaleData = NULL,
        barcodesName = NULL,
        genesName = NULL,
        indicesName = NULL,
        indptrName = NULL,
        modal = c("default", "rna", "atac"),
        featureMeta = NULL,
        ...
) {
    if (!hdf5r::is_hdf5(h5file)) {
        stop("Please specify an HDF5 filename to argument `h5file`.")
    }
    modal <- match.arg(modal)
    additional <- list(...)
    h5file <- hdf5r::H5File$new(h5file, mode = "r+")
    if (!is.null(formatType) &&
        formatType %in% c("10X", "AnnData")) {
        if (formatType == "10X") {
            barcodesName <- "matrix/barcodes"
            barcodes <- h5file[[barcodesName]][]
            rawData <- "matrix/data"
            indicesName <- "matrix/indices"
            indptrName <- "matrix/indptr"
            genesName <- "matrix/features/name"
            genes <- h5file[[genesName]][]
        } else if (formatType == "AnnData") {
            barcodesName <- "obs"
            barcodes <- h5file[[barcodesName]][]$cell
            rawData <- "raw.X/data"
            indicesName <- "raw.X/indices"
            indptrName <- "raw.X/indptr"
            genesName <- "raw.var"
            genes <- h5file[[genesName]][]
        } else {
            stop("Specified `formatType` '", formatType,
                 "' is not supported for now.")
        }
    } else {
        barcodes <- h5file[[barcodesName]][]
        genes <- h5file[[genesName]][]
    }
    # The order of list elements matters. Put "paths" together so easier for
    # checking link existence.
    h5.meta <- list(
        H5File = h5file,
        filename = h5file$filename,
        formatType = formatType,
        indicesName = indicesName,
        indptrName = indptrName,
        barcodesName = barcodesName,
        genesName = genesName,
        rawData = rawData,
        normData = normData,
        scaleData = scaleData
    )
    if (!is.null(rawData)) rawData <- h5file[[rawData]]
    if (!is.null(normData)) normData <- h5file[[normData]]
    if (!is.null(scaleData)) scaleData <- h5file[[scaleData]]
    if (is.null(featureMeta))
        featureMeta <- S4Vectors::DataFrame(row.names = genes)
    else if (!inherits(featureMeta, "DFrame"))
        featureMeta <- S4Vectors::DataFrame(featureMeta)
    allData <- list(.modalClassDict[[modal]],
                    rawData = rawData,
                    normData = normData,
                    scaleData = scaleData,
                    #H = H, V = V, A = A, B = B, U = U,
                    h5fileInfo = h5.meta, featureMeta = featureMeta,
                    colnames = barcodes, rownames = genes)
    allData <- c(allData, additional)
    x <- do.call("new", allData)
    x
}



#' Read 10X alignment data (including V3)
#'
#' This function generates a sparse matrix (genes x cells) from the data generated by 10X's
#' cellranger count pipeline. It can process V2 and V3 data together, producing either a single
#' merged matrix or list of matrices. Also handles multiple data types produced by 10X V3 (Gene
#' Expression, Antibody Capture, CRISPR, CUSTOM).
#'
#' @param sample.dirs List of directories containing either matrix.mtx(.gz) file along with genes.tsv,
#'   (features.tsv), and barcodes.tsv, or outer level 10X output directory (containing outs directory).
#' @param sample.names Vector of names to use for samples (corresponding to sample.dirs)
#' @param merge Whether to merge all matrices of the same data type across samples or leave as list
#'   of matrices (default TRUE).
#' @param num.cells Optional limit on number of cells returned for each sample (only for Gene
#'   Expression data). Retains the cells with the highest numbers of transcripts (default NULL).
#' @param min.umis Minimum UMI threshold for cells (default 0).
#' @param use.filtered Whether to use 10X's filtered data (as opposed to raw). Only relevant for
#'   sample.dirs containing 10X outs directory (default FALSE).
#' @param reference For 10X V<3, specify which reference directory to use if sample.dir is outer
#'   level 10X directory (only necessary if more than one reference used for sequencing).
#'   (default NULL)
#' @param data.type Indicates the protocol of the input data. If not specified, input data will be
#' considered scRNA-seq data (default 'rna', alternatives: 'atac').
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#'
#' @return List of merged matrices across data types (returns sparse matrix if only one data type
#'   detected), or nested list of matrices organized by sample if merge=F.
#'
#' @export
#' @examples
#' \dontrun{
#' # 10X output directory V2 -- contains outs/raw_gene_bc_matrices/<reference>/...
#' sample.dir1 <- "path/to/outer/dir1"
#' # 10X output directory V3 -- for two data types, Gene Expression and CUSTOM
#' sample.dir2 <- "path/to/outer/dir2"
#' dges1 <- read10X(list(sample.dir1, sample.dir2), c("sample1", "sample2"), min.umis = 50)
#' ligerex <- createLiger(expr = dges1[["Gene Expression"]], custom = dges1[["CUSTOM"]])
#' }

read10X <-
    function(sample.dirs,
             sample.names,
             merge = TRUE,
             num.cells = NULL,
             min.umis = 0,
             use.filtered = FALSE,
             reference = NULL,
             data.type = "rna",
             verbose = getOption("ligerVerbose")) {
        datalist <- list()
        datatypes <- c("Gene Expression")

        if (length(num.cells) == 1) {
            num.cells <- rep(num.cells, length(sample.dirs))
        }
        for (i in seq_along(sample.dirs)) {
            print(paste0("Processing sample ", sample.names[i]))
            sample.dir <- sample.dirs[[i]]
            inner1 <- paste0(sample.dir, "/outs")
            if (dir.exists(inner1)) {
                sample.dir <- inner1
                is_v3 <-
                    dir.exists(paste0(sample.dir, "/filtered_feature_bc_matrix"))
                matrix.prefix <- ifelse(use.filtered, "filtered", "raw")
                if (is_v3) {
                    sample.dir <-
                        paste0(sample.dir,
                               "/",
                               matrix.prefix,
                               "_feature_bc_matrix")
                } else {
                    if (is.null(reference)) {
                        references <- list.dirs(
                            paste0(sample.dir, "/raw_gene_bc_matrices"),
                            full.names = FALSE,
                            recursive = FALSE
                        )
                        if (length(references) > 1) {
                            stop(
                                "Multiple reference genomes found. Please specify a single one."
                            )
                        } else {
                            reference <- references[1]
                        }
                    }
                    sample.dir <-
                        paste0(sample.dir,
                               "/",
                               matrix.prefix,
                               "_gene_bc_matrices/",
                               reference)
                }
            } else {
                is_v3 <- file.exists(paste0(sample.dir, "/features.tsv.gz"))
            }
            suffix <- ifelse(is_v3, ".gz", "")
            if (data.type == "rna") {
                features.file <-
                    ifelse(
                        is_v3,
                        paste0(sample.dir, "/features.tsv.gz"),
                        paste0(sample.dir, "/genes.tsv")
                    )
            } else if (data.type == "atac") {
                features.file <- ifelse(
                    is_v3,
                    paste0(sample.dir, "/peaks.bed.gz"),
                    paste0(sample.dir, "/peaks.bed")
                )
            }
            matrix.file <- paste0(sample.dir, "/matrix.mtx", suffix)
            barcodes.file <- paste0(sample.dir, "/barcodes.tsv", suffix)

            rawdata <- Matrix::readMM(matrix.file)
            # convert to dgc matrix
            if (class(rawdata)[1] == "dgTMatrix") {
                rawdata <- methods::as(rawdata, "CsparseMatrix")
            }

            # filter for UMIs first to increase speed
            umi.pass <- which(colSums(rawdata) > min.umis)
            if (length(umi.pass) == 0) {
                message("No cells pass UMI cutoff. Please lower it.")
            }
            rawdata <- rawdata[, umi.pass, drop = FALSE]

            barcodes <- readLines(barcodes.file)[umi.pass]
            # Remove -1 tag from barcodes
            if (all(grepl(barcodes, pattern = "\\-1$"))) {
                barcodes <- as.vector(sapply(barcodes, function(x) {
                    strsplit(x, "-")[[1]][1]
                }))
            }
            if (data.type == "rna") {
                features <-
                    utils::read.delim(features.file,
                               header = FALSE,
                               stringsAsFactors = FALSE)
                rownames(rawdata) <- make.unique(features[, 2])
            } else if (data.type == "atac") {
                features <- utils::read.table(features.file, header = FALSE)
                features <-
                    paste0(features[, 1], ":", features[, 2], "-", features[, 3])
                rownames(rawdata) <- features
            }
            # since some genes are only differentiated by ENSMBL
            colnames(rawdata) <- barcodes

            # split based on 10X datatype -- V3 has Gene Expression, Antibody Capture, CRISPR, CUSTOM
            # V2 has only Gene Expression by default and just two columns
            if (is.null(ncol(features))) {
                samplelist <- list(rawdata)
                names(samplelist) <- c("Chromatin Accessibility")
            } else if (ncol(features) < 3) {
                samplelist <- list(rawdata)
                names(samplelist) <- c("Gene Expression")
            } else {
                sam.datatypes <- features[, 3]
                sam.datatypes.unique <- unique(sam.datatypes)
                # keep track of all unique datatypes
                datatypes <- union(datatypes, unique(sam.datatypes))
                samplelist <- lapply(sam.datatypes.unique, function(x) {
                    rawdata[which(sam.datatypes == x),]
                })
                names(samplelist) <- sam.datatypes.unique
            }

            # num.cells filter only for gene expression data
            if (!is.null(num.cells)) {
                if (names(samplelist) == "Gene Expression" |
                    names(samplelist) == "Chromatin Accessibility") {
                    data_label <- names(samplelist)
                    cs <- colSums(samplelist[[data_label]])
                    limit <- ncol(samplelist[[data_label]])
                    if (num.cells[i] > limit) {
                        if (verbose) {
                            message(
                                "You selected more cells than are in matrix ",
                                i,
                                ". Returning all ",
                                limit,
                                " cells."
                            )
                        }
                        num.cells[i] <- limit
                    }
                    samplelist[[data_label]] <-
                        samplelist[[data_label]][, order(cs, decreasing = TRUE)
                                                 [1:num.cells[i]]]
                }

                # cs <- colSums(samplelist[["Gene Expression"]])
                # limit <- ncol(samplelist[["Gene Expression"]])
                # if (num.cells[i] > limit) {
                #   print(paste0(
                #     "You selected more cells than are in matrix ", i,
                #     ". Returning all ", limit, " cells."
                #   ))
                #   num.cells[i] <- limit
                # }
                # samplelist[["Gene Expression"]] <- samplelist[["Gene Expression"]][, order(cs, decreasing = TRUE)
                #                                                                    [1:num.cells[i]]]
            }

            datalist[[i]] <- samplelist
        }
        if (merge) {
            if (verbose) {
                message("Merging samples")
            }
            return_dges <- lapply(datatypes, function(x) {
                mergelist <- lapply(datalist, function(d) {
                    d[[x]]
                })
                mergelist <- mergelist[!sapply(mergelist, is.null)]
                sample.names.x <-
                    sample.names[!sapply(mergelist, is.null)]
                mergeSparseAll(mergelist, sample.names)
            })
            names(return_dges) <- datatypes

            # if only one type of data present
            if (length(return_dges) == 1) {
                if (verbose) {
                    message("Returning ", datatypes, " data matrix")
                }
                return(return_dges[[1]])
            }
            return(return_dges)
        } else {
            names(datalist) <- sample.names
            return(datalist)
        }
    }

#' Read liger object from RDS file
#' @param filename Path to an RDS file of a \code{liger} object of old versions.
#' @param dimredName The name of variable in \code{cellMeta} slot to store the
#' dimensionality reduction matrix, which originally located in
#' \code{tsne.coords} slot. Default \code{"tsne.coords"}.
#' @param clusterName The name of variable in \code{cellMeta} slot to store the
#' clustering assignment, which originally located in \code{clusters} slot.
#' Default \code{"clusters"}.
#' @param update Logical, whether to update an old (<=1.0.0) \code{liger} object
#' to the currect version of structure. Default \code{TRUE}.
#' @return New version of \linkS4class{liger} object
#' @export
readLiger <- function(
        filename,
        dimredName = "tsne_coords",
        clusterName = "clusters",
        update = TRUE) {
    oldObj <- readRDS(filename)
    if (!inherits(oldObj, "liger"))
        stop("Object is not of class \"liger\".")
    oldVer <- oldObj@version
    if (oldVer >= package_version("1.99.0")) return(oldObj)
    .log("Older version (", oldVer, ") of liger object detected.")
    if (isTRUE(update)) {
        .log("Updating the object structure to make it compatible ",
             "with current version (", utils::packageVersion("rliger"), ")")
        return(convertOldLiger(oldObj, dimredName = dimredName,
                               clusterName = clusterName))
    } else {
        .log("`update = FALSE` specified. Returning the original object.")
        return(oldObj)
    }
}
