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
#' @param addPrefix Logical. Whether to add "<dataset name>_" as a prefix of
#' cell identifiers (e.g. barcodes) to avoid duplicates in multiple libraries (
#' common with 10X data). Default \code{TRUE}
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
        addPrefix = TRUE,
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
        } else {
            datasets[[dname]] <- as.ligerDataset(data, modal = modal[i])
        }
        barcodesOrig <- c(barcodesOrig, colnames(datasets[[dname]]))
        if (isTRUE(addPrefix)) {
            cellID <- paste0(dname, "_", colnames(datasets[[dname]]))
            colnames(datasets[[dname]]) <- cellID
        }
    }

    #barcodesOrig <- unlist(lapply(datasets, colnames), use.names = FALSE)
    datasets <- lapply(datasets, function(ld) {
        colnames(ld) <- make.unique(colnames(ld))
        return(ld)
    })
    cellID <- unlist(lapply(datasets, colnames), use.names = FALSE)
    if (is.null(cellMeta)) {
        cellMeta <- S4Vectors::DataFrame(
            dataset = factor(rep(names(datasets), sapply(datasets, ncol)),
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
            normData <- methods::as(normData, "CsparseMatrix")
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
#' @param h5file Filename of an H5 file
#' @param formatType Select preset of H5 file structure. Default \code{"10X"}.
#' Current available options are \code{"10X"} and \code{"AnnData"}.
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
        formatType = "10X",
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
    if (!is.null(formatType)) {
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


#' Read liger object from RDS file
#' @param filename Path to an RDS file of a \code{liger} object of old versions.
#' @param dimredName The name of variable in \code{cellMeta} slot to store the
#' dimensionality reduction matrix, which originally located in
#' \code{tsne.coords} slot. Default \code{"tsne.coords"}.
#' @param clusterName The name of variable in \code{cellMeta} slot to store the
#' clustering assignment, which originally located in \code{clusters} slot.
#' Default \code{"clusters"}.
#' @param h5FilePath Named list, to specify the path to the H5 file of each
#' dataset if location has been changed. Default \code{NULL} looks at the file
#' paths stored in object.
#' @param update Logical, whether to update an old (<=1.0.0) \code{liger} object
#' to the currect version of structure. Default \code{TRUE}.
#' @return New version of \linkS4class{liger} object
#' @export
readLiger <- function(
        filename,
        dimredName = "tsne_coords",
        clusterName = "clusters",
        h5FilePath = NULL,
        update = TRUE) {
    oldObj <- readRDS(filename)
    if (!inherits(oldObj, "liger"))
        stop("Object is not of class \"liger\".")
    oldVer <- oldObj@version
    if (oldVer >= package_version("1.99.0")) return(oldObj)
    .log("Older version (", oldVer, ") of liger object detected.")
    if (isTRUE(update)) {
        .log("Updating the object structure to make it compatible ",
             "with current version (", utils::packageVersion("rliger2"), ")")
        return(convertOldLiger(oldObj, dimredName = dimredName,
                               clusterName = clusterName,
                               h5FilePath = h5FilePath))
    } else {
        .log("`update = FALSE` specified. Returning the original object.")
        return(oldObj)
    }
}

#' Import prepared dataset publically available
#' @param dataset Name of dataset, see available options with
#' \code{names(.manifest)}.
#' @param overwrite Logical, if a file exists at corresponding download
#' location, whether to re-download or directly use this file. Default
#' \code{FALSE}.
#' @param dir Path to download datasets. Default \code{getwd()}.
#' @param method \code{method} argument directly passed to
#' \code{\link[utils]{download.file}}. Using \code{"libcurl"} while other
#' options might not work depending on platform.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Additional arguments passed to \code{\link{download.file}}
#' @return \code{\linkS4class{liger}} object of specified dataset.
#' @export
#' @examples
#' if (FALSE) {
#'     pbmc <- importVignetteData("pbmc")
#' }
importVignetteData <- function(
        dataset,
        overwrite = FALSE,
        dir = getwd(),
        method = "libcurl",
        verbose = getOption("ligerVerbose"),
        ...
) {
    if (!dataset %in% names(.manifest)) {
        stop("Requested dataset \"", dataset, "\" not supported yet. \n",
             "Available options: ", paste(names(.manifest), collapse = ", "))
    }
    fsep <- ifelse(Sys.info()["sysname"] == "Windows", "\\", "/")

    info <- .manifest[[dataset]]
    url <- sapply(info, function(x) x$url)
    dataNames <- names(info)
    filenames <- sapply(info, function(x) x$filename)
    modal <- sapply(info, function(x) x$modal)

    # ATAC assay specific processing
    # If non of them, peakURL and peakFilename should be empty lists
    atacIdx <- modal == "atac"
    peakURL <- sapply(info[atacIdx], function(x) x$peak$url)

    peakFilename <- sapply(info[atacIdx], function(x) x$peak$filename)
    peakFilename <- file.path(dir, peakFilename, fsep = fsep)
    names(peakFilename) <- names(info[atacIdx])

    filenames <- file.path(dir, filenames, fsep = fsep)

    allURLs <- c(url, peakURL)
    allFiles <- c(filenames, peakFilename)

    doDownload <- rep(TRUE, length(allURLs))
    for (i in seq_along(allURLs)) {
        f <- allFiles[i]
        if (file.exists(f) && isFALSE(overwrite)) {
            warning("File already exists, skipped. set `overwrite = TRUE` ",
                    "to force downloading: ", f)
            doDownload[i] <- FALSE
            next
        }
        if (isTRUE(verbose)) .log("Downloading from ", allURLs[i], " to ", f)
    }
    if (sum(doDownload) > 0) {
        allURLs <- allURLs[doDownload]
        allURLs <- unlist(allURLs)
        utils::download.file(allURLs,
                             destfile = allFiles[doDownload],
                             mode = "wb", quiet = !verbose, method = method,
                             ...)
    }

    rawList <- lapply(filenames, readRDS)
    names(rawList) <- dataNames
    object <- createLiger(rawList, modal = modal, verbose = verbose)
    for (i in which(atacIdx)) {
        name <- names(object)[i]
        rawPeak(object, dataset = name) <- readRDS(peakFilename[name])
    }
    return(object)
}

