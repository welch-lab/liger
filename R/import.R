#' Create liger object
#' @description This function allows creating \linkS4class{liger} object from
#' multiple datasets of various forms (See \code{rawData}).
#' @param rawData Named list of datasets. Required. Elements allowed include a
#' matrix, a \code{Seurat} object, a \code{SingleCellExperiment} object, an
#' \code{AnnData} object, a \linkS4class{ligerDataset} object or a filename to
#' an HDF5 file. See detail for HDF5 reading.
#' @param modal Character vector for modality setting. Currently options of
#' \code{"default"}, \code{"rna"}, \code{"atac"} and \code{"spatial"} are
#' supported.
#' @param cellMeta data.frame of metadata at single-cell level. Default
#' \code{NULL}.
#' @param removeMissing Logical. Whether to remove cells that do not have any
#' counts and features not expressed in any cells from each dataset. Default
#' \code{TRUE}.
#' @param addPrefix Logical. Whether to add "<dataset name>_" as a prefix of
#' cell identifiers (e.g. barcodes) to avoid duplicates in multiple libraries (
#' common with 10X data). Default \code{"auto"} detects if matrix columns
#' already has the exact prefix or not. Logical value forces the action.
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
#' @param newH5 When using HDF5 based data and subsets created after removing
#' missing cells/features, whether to create new HDF5 files for the subset.
#' Default \code{TRUE}. If \code{FALSE}, data will be subset into memory and
#' can be dangerous for large scale analysis.
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
        addPrefix = "auto",
        formatType = "10X",
        dataName = NULL,
        indicesName = NULL,
        indptrName = NULL,
        genesName = NULL,
        barcodesName = NULL,
        newH5 = TRUE,
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
        colnameOrig <- colnames(datasets[[dname]])
        prefix <- paste0(dname, "_")
        .addPrefix <- FALSE
        if (addPrefix == "auto") {
            # If all colnames starts with the prefix wanted, don't add it again
            .addPrefix <- !all(startsWith(colnameOrig, prefix))
        }
        barcodesOrig <- c(barcodesOrig, colnameOrig)
        if (.addPrefix) {
            colnames(datasets[[dname]]) <- paste0(prefix, colnameOrig)
        }
    }

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
                             verbose = verbose, newH5 = newH5)
    }

    return(obj)
}

#' Create in-memory ligerDataset object
#' @param rawData,normData,scaleData A \code{\link[Matrix]{dgCMatrix-class}}
#' object for the raw or normalized expression count or a dense matrix of scaled
#' variable gene expression, respectively. Default \code{NULL} for all three but
#' at lease one has to be specified.
#' @param modal Name of modality for this dataset. Currently options of
#' \code{"default"}, \code{"rna"}, \code{"atac"} and \code{"spatial"} are
#' supported. Default \code{"default"}.
#' @param featureMeta Data frame of feature metadata. Default \code{NULL}.
#' @param ... Additional slot data. See \linkS4class{ligerDataset} for detail.
#' Given values will be directly placed at corresponding slots.
#' @seealso \linkS4class{ligerDataset}, \linkS4class{ligerATACDataset},
#' \code{\link{createH5LigerDataset}}
#' @export
createLigerDataset <- function(
        rawData = NULL,
        modal = c("default", "rna", "atac", "spatial"),
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
        modal = c("default", "rna", "atac", "spatial"),
        featureMeta = NULL,
        ...
) {
    # Currently commented this check because it fails for unknown reason on
    # Windows platform when subsetting H5 to H5, even if the following part
    # works without passing this check
    # if (!hdf5r::is.h5file(h5file)) {
    #     stop("`h5file`: Invalid HDF5 file or file path: ", h5file)
    # }
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
    obj <- readRDS(filename)
    if (!inherits(obj, "liger"))
        stop("Object is not of class \"liger\".")
    ver <- obj@version
    if (ver >= package_version("1.99.0")) {
        if (isH5Liger(obj)) obj <- restoreH5Liger(obj)
        return(obj)
    }
    .log("Older version (", ver, ") of liger object detected.")
    if (isTRUE(update)) {
        .log("Updating the object structure to make it compatible ",
             "with current version (", utils::packageVersion("rliger2"), ")")
        return(convertOldLiger(obj, dimredName = dimredName,
                               clusterName = clusterName,
                               h5FilePath = h5FilePath))
    } else {
        .log("`update = FALSE` specified. Returning the original object.")
        return(obj)
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






#' Load in data from 10X
#' @description
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#'
#' \code{read10X} works generally for 10X cellranger pipelines including:
#' CellRanger < 3.0 & >= 3.0 and CellRanger-ARC.
#'
#' \code{read10XRNA} invokes \code{read10X} and takes the "Gene Expression" out,
#' so that the result can directly be used to construct a \linkS4class{liger}
#' object. See Examples for demonstration.
#'
#' \code{read10XATAC} works for both cellRanger-ARC and cellRanger-ATAC
#' pipelines but needs user arguments for correct recognition. Similarly, the
#' returned value can directly be used for constructing a \linkS4class{liger}
#' object.
#' @param path [A.] A Directory containing the matrix.mtx, genes.tsv (or
#' features.tsv), and barcodes.tsv files provided by 10X. A vector, a named
#' vector, a list or a named list can be given in order to load several data
#' directories. [B.] The 10X root directory where subdirectories of per-sample
#' output folders can be found. Sample names will by default take the name of
#' the vector, list or subfolders.
#' @param sampleNames A vector of names to override the detected or set sample
#' names for what is given to \code{path}. Default \code{NULL}. If no name
#' detected at all and multiple samples are given, will name them by numbers.
#' @param useFiltered Logical, if \code{path} is given as case B, whether to use
#' the filtered feature barcode matrix instead of raw (unfiltered). Default
#' \code{TRUE}.
#' @param reference In case of specifying a CellRanger<3 root folder to
#' \code{path}, import the matrix from the output using which reference. Only
#' needed when multiple references present. Default \code{NULL}.
#' @param geneCol Specify which column of genes.tsv or features.tsv to use for
#' gene names. Default \code{2}.
#' @param cellCol Specify which column of barcodes.tsv to use for cell names.
#' Default \code{1}.
#' @param returnList Logical, whether to still return a structured list instead
#' of a single matrix object, in the case where only one sample and only one
#' feature type can be found. Otherwise will always return a list. Default
#' \code{FALSE}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return
#' \itemize{
#'  \item{When only one sample is given or detected, and only one feature type
#'  is detected or using CellRanger < 3.0, and \code{returnList = FALSE}, a
#'  sparse matrix object (dgCMatrix class) will be returned.}
#'  \item{When using \code{read10XRNA} or \code{read10XATAC}, which are modality
#'  specific, returns a list named by samples, and each element is the
#'  corresponding sparse matrix object (dgCMatrix class).}
#'  \item{\code{read10X} generally returns a list named by samples. Each sample
#'  element will be another list named by feature types even if only one feature
#'  type is detected (or using CellRanger < 3.0) for data structure consistency.
#'  The feature type "Gene Expression" always comes as the first type if
#'  available.}
#' }
#' @export
#' @rdname read10X
#' @examples
#' \dontrun{
#' # For output from CellRanger < 3.0
#' dir <- 'path/to/data/directory'
#' list.files(dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' mat <- read10X(dir)
#' class(mat) # Should show dgCMatrix
#'
#' # For root directory from CellRanger < 3.0
#' dir <- 'path/to/root'
#' list.dirs(dir) # Should show sample names
#' matList <- read10X(dir)
#' names(matList) # Should show the sample names
#' class(matList[[1]][["Gene Expression"]]) # Should show dgCMatrix
#'
#' # For output from CellRanger >= 3.0 with multiple data types
#' dir <- 'path/to/data/directory'
#' list.files(dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
#' matList <- read10X(dir, sampleNames = "tissue1")
#' names(matList) # Shoud show "tissue1"
#' names(matList$tissue1) # Should show feature types, e.g. "Gene Expression" and etc.
#'
#' # For root directory from CellRanger >= 3.0 with multiple data types
#' dir <- 'path/to/root'
#' list.dirs(dir) # Should show sample names, e.g. "rep1", "rep2", "rep3"
#' matList <- read10X(dir)
#' names(matList) # Should show the sample names: "rep1", "rep2", "rep3"
#' names(matList$rep1) # Should show the avalable feature types for rep1
#' }
read10X <- function(
        path,
        sampleNames = NULL,
        useFiltered = NULL,
        reference = NULL,
        geneCol = 2,
        cellCol = 1,
        returnList = FALSE,
        verbose = getOption("ligerVerbose"),
        # Renamed
        sample.dirs = path,
        sample.names = sampleNames,
        use.filtered = useFiltered,
        # Defunct
        data.type = NULL,
        merge = NULL,
        num.cells = NULL,
        min.umis = NULL
) {
    .deprecateArgs(replace = list(sample.dirs = "path",
                                  sample.names = "sampleNames",
                                  use.filtered = "useFiltered"),
                   defunct = c("data.type", "merge", "num.cells", "min.umis"))
    # Adopted from Seurat V4.4.0 repo, with minor modification
    # Whether a folder with samples inside, or a list/vector of mtxDirs
    if (length(path) == 1) {
        dirSampleNames <- list.dirs(path, full.names = FALSE, recursive = FALSE)
        outsPaths <- file.path(path, dirSampleNames, "outs")
        if (any(dir.exists(outsPaths))) {
            # Case B, make the final `path` to case A
            dirSampleNames <- dirSampleNames[dir.exists(outsPaths)]
            outsPaths <- outsPaths[dir.exists(outsPaths)]
            useFiltered <- useFiltered %||% TRUE
            prefix <- ifelse(useFiltered, "filtered", "raw")
            subdir <- paste0(prefix, "_feature_bc_matrix")
            if (dir.exists(file.path(outsPaths[1], subdir))) {
                isV3 <- TRUE
            } else {
                isV3 <- FALSE
                subdir <- paste0(prefix, "_gene_bc_matrix")
            }
            # Now paths are sample/outs/*_*_bc_matrix/
            path <- file.path(outsPaths, subdir)
            if (!isV3) {
                # Account the inner reference folder for V2
                refsExist <- list.dirs(path[1], full.names = FALSE,
                                       recursive = FALSE)
                if (is.null(reference)) {
                    if (length(refsExist) == 1) {
                        reference <- refsExist
                        .log("Using reference: ", reference)
                    } else {
                        stop("Multiple references found, please select one ",
                             "from: ", paste0(refsExist, collapse = ", "))
                    }
                } else if (length(reference) == 1) {
                    if (!reference %in% refsExist) {
                        stop("Specified reference not found, please select ",
                             "one from: ", paste0(refsExist, collapse = ", "))
                    }
                } else {
                    stop("Multiple reference specified but only one allowed.")
                }
                path <- file.path(path, reference)
            }
            names(path) <- dirSampleNames
            .log("Found the following sample folders with possible sub-folder ",
                 "structure: \n", paste0(dirSampleNames, collapse = ", "))
        } # else mtxDirs
    } # else mtxDirs

    allData <- list()
    sampleNames <- .checkArgLen(sampleNames, length(path), repN = FALSE)
    if (is.null(sampleNames) && !is.null(names(path))) {
        sampleNames <- names(path)
    } else {
        if (any(duplicated(sampleNames))) {
            stop("Cannot set duplicated sample names.")
        }
    }

    for (i in seq_along(path)) {
        if (isTRUE(verbose)) {
            name <- sampleNames[i]
            if (is.null(name)) name <- paste0("sample ", i)
            .log("Reading from ", name, "...")
        }
        if (is.list(path)) run <- path[[i]]
        else run <- path[i]

        if (!dir.exists(run)) {
            stop("Directory provided does not exist: ",
                 normalizePath(run, mustWork = FALSE))
        }
        barcode.loc <- file.path(run, 'barcodes.tsv')
        gene.loc <- file.path(run, 'genes.tsv')
        features.loc <- file.path(run, 'features.tsv.gz')
        matrix.loc <- file.path(run, 'matrix.mtx')
        # Flag to indicate if this data is from CellRanger >= 3.0
        isOldVer <- file.exists(gene.loc)
        if (!isOldVer) {
            addgz <- function(s) paste0(s, ".gz")
            barcode.loc <- addgz(barcode.loc)
            matrix.loc <- addgz(matrix.loc)
        }
        if (!file.exists(barcode.loc)) {
            stop("Barcode file missing. Expecting ", basename(barcode.loc))
        }
        if (!isOldVer && !file.exists(features.loc) ) {
            stop("Gene name or features file missing. Expecting ",
                 basename(features.loc))
        }
        if (!file.exists(matrix.loc)) {
            stop("Expression matrix file missing. Expecting ",
                 basename(matrix.loc))
        }
        data <- read10XFiles(matrixPath = matrix.loc, barcodesPath = barcode.loc,
                             featuresPath = ifelse(isOldVer, gene.loc, features.loc),
                             sampleName = sampleNames[i], geneCol = geneCol,
                             cellCol = cellCol)
        if (isOldVer) names(data) <- "Gene Expression"
        allData[[i]] <- data
    }
    if (!is.null(sampleNames)) names(allData) <- sampleNames
    if (isTRUE(returnList)) return(allData)
    # By default, if only one sample and only one feature type,
    # directly return the matrix, otherwise still return the structured list
    if (length(allData) == 1 && length(allData[[1]]) == 1) {
        return(allData[[1]][[1]])
    } else {
        return(allData)
    }
}

#' @rdname read10X
#' @export
#' @param ... Arguments passed to \code{read10X}
#' @examples
#' \dontrun{
#' # For creating LIGER object from root directory of CellRanger >= 3.0
#' dir <- 'path/to/root'
#' list.dirs(dir) # Should show sample names, e.g. "rep1", "rep2", "rep3"
#' matList <- read10XRNA(dir)
#' names(matList) # Should show the sample names: "rep1", "rep2", "rep3"
#' sapply(matList, class) # Should show matrix class all are "dgCMatrix"
#' lig <- createLigerObject(matList)
#' }
read10XRNA <- function(
        path,
        sampleNames = NULL,
        useFiltered = NULL,
        reference = NULL,
        returnList = FALSE,
        ...
) {
    dataList <- read10X(path, sampleNames = sampleNames,
                        useFiltered = useFiltered, reference = reference,
                        returnList = TRUE, ...)
    dataList <- lapply(dataList, `[[`, i = "Gene Expression")
    if (length(dataList) == 1 && isFALSE(returnList)) return(dataList[[1]])
    else return(dataList)
}

#' @rdname read10X
#' @export
#' @param pipeline Which cellRanger pipeline type to find the ATAC data. Choose
#' \code{"atac"} to read the peak matrix from cellranger-atac pipeline output
#' folder(s), or \code{"arc"} to split the ATAC feature subset out from the
#' multiomic cellranger-arc pipeline output folder(s). Default \code{"atac"}.
#' @param arcFeatureType When \code{pipeline = "arc"}, which feature type is
#' for the ATAC data of interests. Default \code{"Peaks"}. Other possible
#' feature types can be \code{"Chromatin Accessibility"}. Error message will
#' show available options if argument specification cannot be found.
read10XATAC <- function(
        path,
        sampleNames = NULL,
        useFiltered = NULL,
        pipeline = c("atac", "arc"),
        arcFeatureType = "Peaks",
        returnList = FALSE,
        geneCol = 2,
        cellCol = 1,
        verbose = getOption("ligerVerbose")
) {
    pipeline <- match.arg(pipeline)
    if (length(path) == 1) {
        dirSampleNames <- list.dirs(path, full.names = FALSE, recursive = FALSE)
        outsPaths <- file.path(path, dirSampleNames, "outs")
        if (any(dir.exists(outsPaths))) {
            # Case B, make the final `path` to case A
            dirSampleNames <- dirSampleNames[dir.exists(outsPaths)]
            outsPaths <- outsPaths[dir.exists(outsPaths)]
            useFiltered <- useFiltered %||% TRUE
            prefix <- ifelse(useFiltered, "filtered", "raw")
            subdir <- paste0(prefix,
                             ifelse(pipeline == "atac",
                                    "_peak_bc_matrix", "_feature_bc_matrix"))
            # Now paths are sample/outs/*_peak_bc_matrix/
            path <- file.path(outsPaths, subdir)
            if (!dir.exists(path)) {
                stop("Cannot find folder '", path, "', not standard ",
                     "`cellranger-", pipeline, "` output. ",
                     "Please try with the other `pipeline`.")
            }
            names(path) <- dirSampleNames
            .log("Found the following sample folders with possible sub-folder ",
                 "structure: \n", paste0(dirSampleNames, collapse = ", "))
        } # else mtxDirs
    } # else mtxDirs

    allData <- list()
    sampleNames <- .checkArgLen(sampleNames, length(path), repN = FALSE)
    if (is.null(sampleNames) && !is.null(names(path))) {
        sampleNames <- names(path)
    } else {
        if (any(duplicated(sampleNames))) {
            stop("Cannot set duplicated sample names.")
        }
    }

    for (i in seq_along(path)) {
        if (isTRUE(verbose)) {
            name <- sampleNames[i]
            if (is.null(name)) name <- paste0("sample ", i)
            .log("Reading from ", name, "...")
        }
        if (is.list(path)) run <- path[[i]]
        else run <- path[i]

        if (!dir.exists(run)) {
            stop("Directory provided does not exist: ",
                 normalizePath(run, mustWork = FALSE))
        }
        barcode.loc <- switch(pipeline,
                              arc = "barcodes.tsv.gz",
                              atac = "barcodes.tsv")
        feature.loc <- switch(pipeline,
                              arc = "features.tsv.gz",
                              atac = "peaks.bed"
        )
        matrix.loc <- switch(pipeline,
                             arc = "matrix.mtx.gz",
                             atac = "matrix.mtx"
        )

        if (!file.exists(barcode.loc)) {
            stop("Barcode file missing. Expecting ", basename(barcode.loc))
        }
        if (!file.exists(feature.loc) ) {
            stop("Peak or feature file missing. Expecting ",
                 basename(feature.loc))
        }
        if (!file.exists(matrix.loc)) {
            stop("Expression matrix file missing. Expecting ",
                 basename(matrix.loc))
        }
        data <- read10XFiles(matrixPath = matrix.loc,
                             barcodesPath = barcode.loc,
                             featuresPath = feature.loc,
                             sampleName = sampleNames[i],
                             geneCol = geneCol, cellCol = cellCol,
                             isATAC = pipeline == "atac")
        if (pipeline == "arc" && !arcFeatureType %in% names(data)) {
            stop("No ATAC data retrieved from cellranger-arc pipeline. ",
                 "Please see if the following available feature types match ",
                 "with need and select one for `arcFeatureType`: ",
                 paste0(names(data), collapse = ", "))
        }
        data <- switch(pipeline,
                       arc = data[[arcFeatureType]],
                       atac = data[[1]]
        )
        allData[[i]] <- data
    }
    if (!is.null(sampleNames)) names(allData) <- sampleNames
    if (isTRUE(returnList)) return(allData)
    # By default, if only one sample and only one feature type,
    # directly return the matrix, otherwise still return the structured list
    if (length(allData) == 1 && length(allData[[1]]) == 1) {
        return(allData[[1]][[1]])
    } else {
        return(allData)
    }
}

read10XFiles <- function(
        matrixPath,
        barcodesPath,
        featuresPath,
        sampleName = NULL,
        geneCol = 2,
        cellCol = 1,
        isATAC = FALSE
) {
    # Matrix can be easily read
    data <- methods::as(Matrix::readMM(matrixPath), "CsparseMatrix")

    # Processing barcodes
    cell.barcodes <- read.table(barcodesPath, header = FALSE,
                                sep = '\t', row.names = NULL)
    if (ncol(cell.barcodes) > 1) {
        cell.names <- cell.barcodes[, cellCol]
    } else {
        cell.names <- readLines(barcodesPath)
    }

    # If sampleNames given, prefix with it; otherwise looks at list/vector
    # names; finally number prefix if more than one sample.
    # `sampleNames` is preprocessed and checked at top of the function
    if (is.null(sampleName)) {
        colnames(data) <- cell.names
    } else {
        colnames(data) <- paste0(sampleName, "_", cell.names)
    }

    # Processing features
    feature.names <- utils::read.delim(
        file = featuresPath,
        header = FALSE,
        stringsAsFactors = FALSE
    )
    if (isTRUE(isATAC)) {
        features <- paste0(feature.names[, 1], ":", feature.names[, 2],
                           "-", feature.names[, 3])
    } else {
        if (ncol(feature.names) < geneCol) {
            stop("`geneCol` was set to ", geneCol, " but feature.tsv.gz ",
                 "(or genes.tsv) only has ", ncol(feature.names), " columns.",
                 " Try setting `geneCol` to a value <= ",
                 ncol(feature.names), ".")
        }
        if (any(is.na(feature.names[, geneCol]))) {
            warning(
                "Some features names are NA. Replacing NA names with ID from the ",
                "opposite column requested",
                call. = FALSE, immediate. = TRUE
            )
            na.features <- which(is.na(feature.names[, geneCol]))
            replacement.column <- ifelse(geneCol == 2, 1, 2)
            feature.names[na.features, geneCol] <-
                feature.names[na.features, replacement.column]
        }
        features <- feature.names[, geneCol]
    }
    rownames(data) <- make.unique(feature.names[, geneCol])
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(feature.names) > 2) {
        data_types <- factor(feature.names$V3)
        lvls <- levels(data_types)
        if (length(lvls) > 1) {
            .log("10X data contains more than one type and is being ",
                 "returned as a list containing matrices of each type.")
        }
        expr_name <- "Gene Expression"
        # Return Gene Expression first
        if (expr_name %in% lvls) {
            lvls <- c(expr_name, lvls[-which(lvls == expr_name)])
        }
        data <- lapply(lvls, function(l) {
            data[data_types == l, , drop = FALSE]
        })
        names(data) <- lvls
    } else{
        data <- list(data)
    }
    return(data)
}
