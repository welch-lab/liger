#' Create liger object
#' @description This function allows creating \linkS4class{liger} object from
#' multiple datasets of various forms (See \code{rawData}).
#'
#' \bold{DO} make a copy of the H5AD files because rliger functions write to
#' the files and they will not be able to be read back to Python. This will be
#' fixed in the future.
#' @param rawData Named list of datasets. Required. Elements allowed include a
#' matrix, a \code{Seurat} object, a \code{SingleCellExperiment} object, an
#' \code{AnnData} object, a \linkS4class{ligerDataset} object or a filename to
#' an HDF5 file. See detail for HDF5 reading.
#' @param modal Character vector for modality setting. Use one string for all
#' datasets, or the same number of strings as the number of datasets. Currently
#' options of \code{"default"}, \code{"rna"}, \code{"atac"}, \code{"spatial"}
#' and \code{"meth"} are supported.
#' @param organism Character vector for setting organism for identifying mito,
#' ribo and hemo genes for expression percentage calculation. Use one string for
#' all datasets, or the same number of strings as the number of datasets.
#' Currently options of \code{"mouse"}, \code{"human"}, \code{"zebrafish"},
#' \code{"rat"}, and \code{"drosophila"} are supported.
#' @param cellMeta data.frame of metadata at single-cell level. Default
#' \code{NULL}.
#' @param removeMissing Logical. Whether to remove cells that do not have any
#' counts from each dataset. Default \code{TRUE}.
#' @param addPrefix Logical. Whether to add "datasetName_" as a prefix of
#' cell identifiers (e.g. barcodes) to avoid duplicates in multiple libraries (
#' common with 10X data). Default \code{"auto"} detects if matrix columns
#' already has the exact prefix or not. Logical value forces the action.
#' @param formatType Select preset of H5 file structure. Current available
#' options are \code{"10x"} and \code{"anndata"}. Can be either a single
#' specification for all datasets or a character vector that match with each
#' dataset.
#' @param anndataX The HDF5 path to the raw count data in an H5AD file. See
#' \code{\link{createH5LigerDataset}} Details. Default \code{"X"}.
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
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param ... Additional slot values that should be directly placed in object.
#' @param raw.data,remove.missing,format.type,data.name,indices.name,indptr.name,genes.name,barcodes.name
#' `r lifecycle::badge("superseded")` See Usage section for replacement.
#' @param take.gene.union `r lifecycle::badge("defunct")` Will be ignored.
#' @export
#' @seealso \code{\link{createLigerDataset}}, \code{\link{createH5LigerDataset}}
#' @examples
#' # Create from raw count matrices
#' ctrl.raw <- rawData(pbmc, "ctrl")
#' stim.raw <- rawData(pbmc, "stim")
#' pbmc1 <- createLiger(list(ctrl = ctrl.raw, stim = stim.raw))
#'
#' # Create from H5 files
#' h5Path <- system.file("extdata/ctrl.h5", package = "rliger")
#' tempPath <- tempfile(fileext = ".h5")
#' file.copy(from = h5Path, to = tempPath)
#' lig <- createLiger(list(ctrl = tempPath))
#'
#' # Create from other container object
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'     ctrl.seu <- SeuratObject::CreateSeuratObject(ctrl.raw)
#'     stim.seu <- SeuratObject::CreateSeuratObject(stim.raw)
#'     pbmc2 <- createLiger(list(ctrl = ctrl.seu, stim = stim.seu))
#' }
createLiger <- function(
        rawData,
        modal = NULL,
        organism = "human",
        cellMeta = NULL,
        removeMissing = TRUE,
        addPrefix = "auto",
        formatType = "10X",
        anndataX = "X",
        dataName = NULL,
        indicesName = NULL,
        indptrName = NULL,
        genesName = NULL,
        barcodesName = NULL,
        newH5 = TRUE,
        verbose = getOption("ligerVerbose", TRUE),
        ...,
        # Deprecated coding style
        raw.data = rawData,
        take.gene.union = NULL,
        remove.missing = removeMissing,
        format.type = formatType,
        data.name = dataName,
        indices.name = indicesName,
        indptr.name = indptrName,
        genes.name = genesName,
        barcodes.name = barcodesName
) {
    .deprecateArgs(list(raw.data = "rawData", remove.missing = "removeMissing",
                        format.type = "formatType", data.name = "dataName",
                        indices.name = "indicesName",
                        indptr.name = "indptrName", genes.name = "genesName",
                        barcodes.name = "barcodesName"),
                   defunct = "take.gene.union")
    if (!is.list(rawData) ||
        is.null(names(rawData)) ||
        any(nchar(names(rawData)) == 0)) {
        cli::cli_abort("{.var rawData} has to be a named list.")
    }

    nData <- length(rawData)
    if (missing(modal) || is.null(modal)) modal <- "default"
    modal <- tolower(modal)
    modal <- .checkArgLen(modal, nData, repN = TRUE, class = "character")

    organism <- tolower(organism)
    organism <- .checkArgLen(organism, nData, repN = TRUE, class = "character")
    if (!all(organism %in% c("mouse", "human", "zebrafish", "rat", "drosophila"))) {
        cli::cli_abort("Invalid {.var organism} value. Only support: mouse, human, zebrafish, rat, drosophila.")
    }

    # TODO handle h5 specific argument for hybrid of H5 and in memory stuff.
    datasets <- list()
    barcodesOrig <- NULL
    cellID <- NULL
    for (i in seq_along(rawData)) {
        dname <- names(rawData)[i]
        data <- rawData[[i]]
        if (is.character(data)) {
            # Assuming character input is a filename
             ld <- createH5LigerDataset(
                h5file = data,
                formatType = formatType,
                anndataX = anndataX,
                rawData = dataName,
                barcodesName = barcodesName,
                genesName = genesName,
                indicesName = indicesName,
                indptrName = indptrName,
                modal = modal[i]
            )
        } else {
            ld <- as.ligerDataset(data, modal = modal[i])
        }
        if (is.null(ld)) next
        colnameOrig <- colnames(ld)
        prefix <- paste0(dname, "_")
        .addPrefix <- FALSE
        if (addPrefix == "auto") {
            # If all colnames starts with the prefix wanted, don't add it again
            .addPrefix <- !all(startsWith(colnameOrig, prefix))
        }
        barcodesOrig <- c(barcodesOrig, colnameOrig)
        if (.addPrefix) {
            colnames(ld) <- paste0(prefix, colnameOrig)
        }
        datasets[[dname]] <- ld
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
    for (i in seq_along(datasets)) {
        obj <- runGeneralQC(
            object = obj,
            organism = organism[i],
            useDatasets = i,
            overwrite = TRUE,
            verbose = verbose
        )
    }

    if (isTRUE(removeMissing)) {
        obj <- removeMissing(obj, "cell", filenameSuffix = "qc",
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
#' \code{"default"}, \code{"rna"}, \code{"atac"}, \code{"spatial"} and
#' \code{"meth"} are supported. Default \code{"default"}.
#' @param featureMeta Data frame of feature metadata. Default \code{NULL}.
#' @param ... Additional slot data. See \linkS4class{ligerDataset} for detail.
#' Given values will be directly placed at corresponding slots.
#' @seealso \linkS4class{ligerDataset}, \linkS4class{ligerATACDataset},
#' \linkS4class{ligerSpatialDataset}, \linkS4class{ligerMethDataset}
#' @export
#' @examples
#' ctrl.raw <- rawData(pbmc, "ctrl")
#' ctrl.ld <- createLigerDataset(ctrl.raw)
createLigerDataset <- function(
        rawData = NULL,
        modal = c("default", "rna", "atac", "spatial", "meth"),
        normData = NULL,
        scaleData = NULL,
        featureMeta = NULL,
        ...
) {
    modal <- match.arg(modal)
    args <- as.list(environment())
    additional <- list(...)
    # Necessary initialization of slots
    # if (is.null(rawData) && is.null(normData)) {
    #     cli::cli_abort("At least one of {.field rawData} or {.field normData} has to be provided.")
    # }
    # Look for proper colnames and rownames
    cn <- NULL
    rn <- NULL
    for (i in c("rawData", "normData", "scaleData")) {
        cn <- colnames(args[[i]])
        if (!is.null(cn)) break
    }
    if (!is.null(rawData)) {
        rn <- rownames(rawData)
        if (!inherits(rawData, "dgCMatrix") && !inherits(rawData, "DelayedArray"))
            rawData <- methods::as(rawData, "CsparseMatrix")
    }
    if (!is.null(normData)) {
        if (is.null(rn)) rn <- rownames(normData)
        if (!inherits(normData, "dgCMatrix") && !inherits(normData, "DelayedArray"))
            normData <- methods::as(normData, "CsparseMatrix")
    }
    if (!is.null(scaleData)) {
        if (is.null(rn)) rn <- rownames(scaleData)
        if (!inherits(scaleData, "dgCMatrix")  && !inherits(scaleData, "DelayedArray"))
            scaleData <- methods::as(scaleData, "CsparseMatrix")
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
#' @description
#' For convenience, the default \code{formatType = "10x"} directly fits the
#' structure of cellranger output. \code{formatType = "anndata"} works for
#' current AnnData H5AD file specification (see Details). If a customized H5
#' file structure is presented, any of the \code{rawData},
#' \code{indicesName}, \code{indptrName}, \code{genesName}, \code{barcodesName}
#' should be specified accordingly to override the \code{formatType} preset.
#'
#' \bold{DO} make a copy of the H5AD files because rliger functions write to
#' the files and they will not be able to be read back to Python. This will be
#' fixed in the future.
#' @details
#' For H5AD file written from an AnnData object, we allow using
#' \code{formatType = "anndata"} for the function to infer the proper structure.
#' However, while a typical AnnData-based analysis tends to in-place update the
#' \code{adata.X} attribute and there is no standard/forced convention for where
#' the raw count data, as needed from LIGER, is stored. Therefore, we expose
#' argument \code{anndataX} for specifying this information. The default value
#' \code{"X"} looks for \code{adata.X}. If the raw data is stored in a layer,
#' e.g. \code{adata.layers['count']}, then \code{anndataX = "layers/count"}.
#' If it is stored to \code{adata.raw.X}, then \code{anndataX = "raw/X"}. If
#' your AnnData object does not have the raw count retained, you will have to
#' go back to the Python work flow to have it inserted at desired object space
#' and re-write the H5AD file, or just go from upstream source files with which
#' the AnnData was originally created.
#' @param h5file Filename of an H5 file
#' @param formatType Select preset of H5 file structure. Default \code{"10X"}.
#' Alternatively, we also support \code{"anndata"} for H5AD files.
#' @param rawData,indicesName,indptrName The path in a H5 file for the raw
#' sparse matrix data. These three types of data stands for the \code{x},
#' \code{i}, and \code{p} slots of a \code{\link[Matrix]{dgCMatrix-class}}
#' object. Default \code{NULL} uses \code{formatType} preset.
#' @param normData The path in a H5 file for the "x" vector of the normalized
#' sparse matrix. Default \code{NULL}.
#' @param scaleData The path in a H5 file for the Group that contains the sparse
#' matrix constructing information for the scaled data. Default \code{NULL}.
#' @param genesName,barcodesName The path in a H5 file for the gene names and
#' cell barcodes. Default \code{NULL} uses \code{formatType} preset.
#' @param anndataX The HDF5 path to the raw count data in an H5AD file. See
#' Details. Default \code{"X"}.
#' @param modal Name of modality for this dataset. Currently options of
#' \code{"default"}, \code{"rna"}, \code{"atac"}, \code{"spatial"} and
#' \code{"meth"} are supported. Default \code{"default"}.
#' @param featureMeta Data frame for feature metadata. Default \code{NULL}.
#' @param ... Additional slot data. See \linkS4class{ligerDataset} for detail.
#' Given values will be directly placed at corresponding slots.
#' @export
#' @return H5-based \linkS4class{ligerDataset} object
#' @examples
#' h5Path <- system.file("extdata/ctrl.h5", package = "rliger")
#' tempPath <- tempfile(fileext = ".h5")
#' file.copy(from = h5Path, to = tempPath)
#' ld <- createH5LigerDataset(tempPath)
createH5LigerDataset <- function(
        h5file,
        formatType = "10x",
        rawData = NULL,
        normData = NULL,
        scaleData = NULL,
        barcodesName = NULL,
        genesName = NULL,
        indicesName = NULL,
        indptrName = NULL,
        anndataX = "X",
        modal = c("default", "rna", "atac", "spatial", "meth"),
        featureMeta = NULL,
        ...
) {
    replaceFunc <- ifelse(endsWith(h5file, 'h5ad'), 'readH5AD', 'read10XH5')
    lifecycle::deprecate_warn(
        when = "2.2.0",
        what = I('Building liger object with specifying paths to HDF5 files'),
        with = sprintf('%s()', replaceFunc),
        details = sprintf("It is now possible to load HDF5 backed data into memory with %s() by default, or featuring the on-disk scalability with %s(inMemory = FALSE)",
                          replaceFunc, replaceFunc)
    )
    # Currently commented this check because it fails for unknown reason on
    # Windows platform when subsetting H5 to H5, even if the following part
    # works without passing this check
    # if (!hdf5r::is.h5file(h5file)) {
    #     stop("`h5file`: Invalid HDF5 file or file path: ", h5file)
    # }
    modal <- match.arg(modal)
    additional <- list(...)
    h5file <- hdf5r::H5File$new(h5file, mode = "r+")
    tryCatch(
        {
            if (!is.null(formatType)) {
                formatType <- tolower(formatType)
                if (formatType == "10x") {
                    barcodesName <- barcodesName %||% "matrix/barcodes"
                    barcodes <- h5file[[barcodesName]][]
                    rawData <- rawData %||% "matrix/data"
                    indicesName <- indicesName %||% "matrix/indices"
                    indptrName <- indptrName %||% "matrix/indptr"
                    genesName <- genesName %||% "matrix/features/name"
                    genes <- h5file[[genesName]][]
                } else if (formatType == "anndata") {
                    # AnnData specs changed around 0.7.0
                    if (inherits(h5file[["obs"]], "H5Group")) {
                        barcodesName <- paste0("obs/", hdf5r::h5attr(h5file[["obs"]], "_index"))
                        barcodes <- h5file[[barcodesName]][]
                        genesName <- paste0("var/", hdf5r::h5attr(h5file[["var"]], "_index"))
                        genes <- h5file[[genesName]][]
                    } else if (inherits(h5file[["obs"]], "H5D")) {
                        barcodesName <- "obs"
                        barcodes <- h5file[["obs"]][]$cell
                        genesName <- genesName %||% "raw.var"
                        genes <- h5file[[genesName]][]
                    }
                    rawData <- rawData %||% paste0(anndataX, "/data")
                    indicesName <- indicesName %||% paste0(anndataX, "/indices")
                    indptrName <- indptrName %||% paste0(anndataX, "/indptr")
                    if (!inherits(h5file[[rawData]]$key_info$type, "H5T_INTEGER")) {
                        cli::cli_alert_warning(
                            "Warning: H5AD matrix data detected is not of integer dtype while {.pkg rliger} requires raw counts input."
                        )
                        cli::cli_alert_warning(
                            "Use caution when using this dataset if the data really does not have integer values."
                        )
                        cli::cli_alert_info(
                            "A different matrix can be selected using argument {.var anndataX}."
                        )
                        cli::cli_alert_info(
                            "See {.code ?createH5LigerDataset} {.emph Details} for more information."
                        )
                    }
                } else {
                    cli::cli_abort("Specified {.var formatType} ({.val {formatType}}) is not supported for now.")
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
            allData <- list(Class = .modalClassDict[[modal]],
                            rawData = rawData,
                            normData = normData,
                            scaleData = scaleData,
                            #H = H, V = V, A = A, B = B, U = U,
                            h5fileInfo = h5.meta, featureMeta = featureMeta,
                            colnames = barcodes, rownames = genes)
            allData <- c(allData, additional)
            x <- do.call("new", allData)
            return(x)
        },
        error = function(e) {
            message(format(e))
            cli::cli_alert_warning("Closing all linking to H5 file: {.file {h5file$filename}}")
            h5file$close_all()
            return(invisible(NULL))
        }
    )
}


#' Read liger object from RDS file
#' @description
#' This file reads a liger object stored in RDS files under all kinds of types.
#' 1. A \linkS4class{liger} object with in-memory data created from package
#' version since 1.99.
#' 2. A liger object with on-disk H5 data associated, where the link to H5 files
#' will be automatically restored.
#' 3. A liger object created with older package version, and can be updated to
#' the latest data structure by default.
#' @param filename Path to an RDS file of a \code{liger} object of old versions.
#' @param dimredName The name of variable in \code{cellMeta} slot to store the
#' dimensionality reduction matrix, which originally located in
#' \code{tsne.coords} slot. Default \code{"tsne.coords"}.
#' @param clusterName The name of variable in \code{cellMeta} slot to store the
#' clustering assignment, which originally located in \code{clusters} slot.
#' Default \code{"clusters"}.
#' @param h5FilePath Named character vector for all H5 file paths. Not required
#' for object run with in-memory analysis. For object containing H5-based
#' analysis (e.g. online iNMF), this must be supplied if the H5 file location is
#' different from that at creation time.
#' @param update Logical, whether to update an old (<=1.99.0) \code{liger} object
#' to the currect version of structure. Default \code{TRUE}.
#' @return New version of \linkS4class{liger} object
#' @export
#' @examples
#' # Save and read regular current-version liger object
#' tempPath <- tempfile(fileext = ".rds")
#' saveRDS(pbmc, tempPath)
#'
#' pbmc <- readLiger(tempPath, dimredName = NULL)
#'
#' # Save and read H5-based liger object
#' h5Path <- system.file("extdata/ctrl.h5", package = "rliger")
#' h5tempPath <- tempfile(fileext = ".h5")
#' file.copy(from = h5Path, to = h5tempPath)
#' lig <- createLiger(list(ctrl = h5tempPath))
#' tempPath <- tempfile(fileext = ".rds")
#' saveRDS(lig, tempPath)
#'
#' lig <- readLiger(tempPath, h5FilePath = list(ctrl = h5tempPath))
#'
#' \dontrun{
#' # Read a old liger object <= 1.0.1
#' # Assume the dimensionality reduction method applied was UMAP
#' # Assume the clustering was derived with Louvain method
#' lig <- readLiger(
#'     filename = "path/to/oldLiger.rds",
#'     dimredName = "UMAP",
#'     clusterName = "louvain"
#' )
#' }
readLiger <- function(
        filename,
        dimredName,
        clusterName = "clusters",
        h5FilePath = NULL,
        update = TRUE) {
    object <- readRDS(filename)
    if (isTRUE(update)) {
        object <- updateLigerObject(object, dimredName, clusterName, h5FilePath)
    }
    return(object)
}

# nocov start
#' Import prepared dataset publically available
#' @description
#' These are functions to download example datasets that are subset from public
#' data.
#' \itemize{
#' \item{\bold{PBMC} - Downsampled from GSE96583, Kang et al, Nature
#' Biotechnology, 2018. Contains two scRNAseq datasets.}
#' \item{\bold{BMMC} - Downsampled from GSE139369, Granja et al, Nature
#' Biotechnology, 2019. Contains two scRNAseq datasets and one scATAC data.}
#' \item{\bold{CGE} - Downsampled from GSE97179, Luo et al, Science, 2017.
#' Contains one scRNAseq dataset and one DNA methylation data.}
#' }
#'
#' @rdname importVignetteData
#' @param overwrite Logical, if a file exists at corresponding download
#' location, whether to re-download or directly use this file. Default
#' \code{FALSE}.
#' @param dir Path to download datasets. Default current working directory
#' \code{getwd()}.
#' @param method \code{method} argument directly passed to
#' \code{\link[utils]{download.file}}. Using \code{"libcurl"} while other
#' options might not work depending on platform.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param ... Additional arguments passed to \code{\link{download.file}}
#' @export
#' @return Constructed \linkS4class{liger} object with QC performed and missing
#' data removed.
#' @examplesIf interactive()
#' \donttest{
#' pbmc <- importPBMC()
#' bmmc <- importBMMC()
#' cge <- importCGE()
#' }
importPBMC <- function(
        dir = getwd(),
        overwrite = FALSE,
        method = "libcurl",
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    fsep <- ifelse(Sys.info()["sysname"] == "Windows", "\\", "/")
    info <- data.frame(
        name = c("ctrl", "stim"),
        url = c(
            "https://figshare.com/ndownloader/files/40054042",
            "https://figshare.com/ndownloader/files/40054048"
        ),
        filename = file.path(dir,
                             c("liger_PBMCs_ctrl.rds", "liger_PBMCs_stim.rds"),
                             fsep = fsep),
        modal = "default"
    )
    doDownload <- rep(TRUE, nrow(info))
    for (i in seq(nrow(info))) {
        f <- info$filename[i]
        if (file.exists(f) && isFALSE(overwrite)) {
            cli::cli_alert_warning(
                "Skipping file already exists at: {.file {f}}. "
            )
            cli::cli_alert_info("Set {.code overwrite = TRUE} to forcing download.")
            doDownload[i] <- FALSE
            next
        }
        if (isTRUE(verbose))
            cli::cli_alert_info("Downloading from {.url {info$url[i]}} to {.file {f}}")
    }
    if (sum(doDownload) > 0) {
        utils::download.file(info$url[doDownload],
                             destfile = info$filename[doDownload],
                             mode = "wb", quiet = !verbose, method = method,
                             ...)
    }
    rawList <- lapply(info$filename, readRDS)
    names(rawList) <- info$name
    object <- createLiger(rawList, modal = info$modal, verbose = verbose)
    return(object)
}

#' @rdname importVignetteData
#' @export
importBMMC <- function(
        dir = getwd(),
        overwrite = FALSE,
        method = "libcurl",
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    fsep <- ifelse(Sys.info()["sysname"] == "Windows", "\\", "/")
    info <- data.frame(
        name = c("rna_D1T1", "rna_D1T2", "atac_D5T1", NA),
        url = c(
            "https://figshare.com/ndownloader/files/40054858",
            "https://figshare.com/ndownloader/files/40054861",
            "https://figshare.com/ndownloader/files/40054891",
            "https://figshare.com/ndownloader/files/40054864"
        ),
        filename = file.path(dir,
                             c("liger_BMMC_rna_D1T1.rds",
                               "liger_BMMC_rna_D1T2.rds",
                               "liger_BMMC_atac_D5T1.rds",
                               "liger_BMMC_atac_D5T1_peak.rds"),
                             fsep = fsep),
        modal = c("default", "default", "atac", NA)
    )
    doDownload <- rep(TRUE, nrow(info))
    for (i in seq(nrow(info))) {
        f <- info$filename[i]
        if (file.exists(f) && isFALSE(overwrite)) {
            cli::cli_alert_warning(
                "Skipping file already exists at: {.file {f}}. "
            )
            cli::cli_alert_info("Set {.code overwrite = TRUE} to forcing download.")
            doDownload[i] <- FALSE
            next
        }
        if (isTRUE(verbose))
            cli::cli_alert_info("Downloading from {.url {info$url[i]}} to {.file {f}}")
    }
    if (sum(doDownload) > 0) {
        utils::download.file(info$url[doDownload],
                             destfile = info$filename[doDownload],
                             mode = "wb", quiet = !verbose, method = method,
                             ...)
    }
    rawList <- lapply(info$filename, readRDS)
    names(rawList) <- info$name
    # BMMC raw data has prefix added to barcodes, so don't add it
    object <- createLiger(rawList[1:3], modal = info$modal[1:3], verbose = verbose,
                          addPrefix = FALSE)
    rawPeak(object, info$name[3]) <- rawList[[4]]

    return(object)
}

#' @rdname importVignetteData
#' @export
importCGE <- function(
        dir = getwd(),
        overwrite = FALSE,
        method = "libcurl",
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    fsep <- ifelse(Sys.info()["sysname"] == "Windows", "\\", "/")
    info <- data.frame(
        name = c("rna", "met"),
        url = c(
            "https://figshare.com/ndownloader/files/40222699",
            "https://figshare.com/ndownloader/files/40222702"
        ),
        filename = file.path(dir,
                             c("liger_CGE_rna.rds", "liger_CGE_met.rds"),
                             fsep = fsep),
        modal = c("default", "meth")
    )
    doDownload <- rep(TRUE, nrow(info))
    for (i in seq(nrow(info))) {
        f <- info$filename[i]
        if (file.exists(f) && isFALSE(overwrite)) {
            cli::cli_alert_warning(
                "Skipping file already exists at: {.file {f}}. "
            )
            cli::cli_alert_info("Set {.code overwrite = TRUE} to forcing download.")
            doDownload[i] <- FALSE
            next
        }
        if (isTRUE(verbose))
            cli::cli_alert_info("Downloading from {.url {info$url[i]}} to {.file {f}}")
    }
    if (sum(doDownload) > 0) {
        utils::download.file(info$url[doDownload],
                             destfile = info$filename[doDownload],
                             mode = "wb", quiet = !verbose, method = method,
                             ...)
    }
    rawList <- lapply(info$filename, readRDS)
    names(rawList) <- info$name
    object <- createLiger(rawList, modal = info$modal, verbose = verbose)
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
#' @param path (A.) A Directory containing the matrix.mtx, genes.tsv (or
#' features.tsv), and barcodes.tsv files provided by 10X. A vector, a named
#' vector, a list or a named list can be given in order to load several data
#' directories. (B.) The 10X root directory where subdirectories of per-sample
#' output folders can be found. Sample names will by default take the name of
#' the vector, list or subfolders.
#' @param sampleNames A vector of names to override the detected or set sample
#' names for what is given to \code{path}. Default \code{NULL}. If no name
#' detected at all and multiple samples are given, will name them by numbers.
#' @param addPrefix Logical, whether to add sample names as a prefix to the
#' barcodes. Default \code{FALSE}.
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
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param sample.dirs,sample.names,use.filtered These arguments are renamed and
#' will be deprecated in the future. Please see usage for corresponding
#' arguments.
#' @param data.type,merge,num.cells,min.umis These arguments are defuncted
#' because the functionality can/should be fulfilled with other functions.
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
        addPrefix = FALSE,
        useFiltered = NULL,
        reference = NULL,
        geneCol = 2,
        cellCol = 1,
        returnList = FALSE,
        verbose = getOption("ligerVerbose", TRUE),
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
                        cli::cli_alert_info("Using referece {.val {reference}}")
                    } else {
                        cli::cli_abort(
                            "Multiple references found, please select one from: {.val {refsExist}}"
                        )
                    }
                } else if (length(reference) == 1) {
                    if (!reference %in% refsExist) {
                        cli::cli_abort(
                            "Specified reference not found, please select one from: {.val {refsExist}}"
                        )
                    }
                } else {
                    cli::cli_abort("Multiple reference specified but only one allowed.")
                }
                path <- file.path(path, reference)
            }
            names(path) <- dirSampleNames
            cli::cli_alert_info(
                c("Found the following sample folders with possible sub-folder structure: ",
                  "{.val {dirSampleNames}}")
            )
        } # else mtxDirs
    } # else mtxDirs

    allData <- list()
    sampleNames <- .checkArgLen(sampleNames, length(path), repN = FALSE, class = "character")
    if (is.null(sampleNames) && !is.null(names(path))) {
        sampleNames <- names(path)
    } else {
        if (any(duplicated(sampleNames))) {
            cli::cli_abort("Cannot set duplicated sample names.")
        }
    }

    for (i in seq_along(path)) {
        if (isTRUE(verbose)) {
            name <- sampleNames[i]
            if (is.null(name)) name <- paste0("sample ", i)
            cliID <- cli::cli_process_start("Reading from {.val {name}}")
        }
        if (is.list(path)) run <- path[[i]]
        else run <- path[i]

        if (!dir.exists(run)) {
            cli::cli_abort("Directory provided does not exist: {.file {normalizePath(run, mustWork = FALSE)}}")
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
            cli::cli_abort("Barcode file is missing. Expecting {.file {barcode.loc}}")
        }
        if (!isOldVer && !file.exists(features.loc) ) {
            cli::cli_abort("Gene name or features file is missing. Expecting {.file {features.loc}}")
        }
        if (!file.exists(matrix.loc)) {
            cli::cli_abort("Expression matrix file is missing. Expecting {.file {matrix.loc}}")
        }
        data <- read10XFiles(matrixPath = matrix.loc, barcodesPath = barcode.loc,
                             featuresPath = ifelse(isOldVer, gene.loc, features.loc),
                             sampleName = if (isTRUE(addPrefix)) sampleNames[i] else NULL,
                             geneCol = geneCol, cellCol = cellCol, returnList = TRUE)
        if (isOldVer) names(data) <- "Gene Expression"
        allData[[i]] <- data
        if (isTRUE(verbose)) cli::cli_process_done(id = cliID)
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
        addPrefix = FALSE,
        useFiltered = NULL,
        reference = NULL,
        returnList = FALSE,
        ...
) {
    dataList <- read10X(path, sampleNames = sampleNames, addPrefix = addPrefix,
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
        addPrefix = FALSE,
        useFiltered = NULL,
        pipeline = c("atac", "arc"),
        arcFeatureType = "Peaks",
        returnList = FALSE,
        geneCol = 2,
        cellCol = 1,
        verbose = getOption("ligerVerbose", TRUE)
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
                cli::cli_abort(
                    c("Cannot find folder {.file {path}}, not standard {.code cellranger-{pipeline}} output. ",
                      "i" = "Please try with the other {.code pipeline}.")
                )
            }
            names(path) <- dirSampleNames
            cli::cli_alert_info(
                "Found the following sample folders with possible sub-folder structure: {.val {dirSampleNames}}"
            )
        } # else mtxDirs
    } # else mtxDirs

    allData <- list()
    sampleNames <- .checkArgLen(sampleNames, length(path), repN = FALSE, class = "character")
    if (is.null(sampleNames) && !is.null(names(path))) {
        sampleNames <- names(path)
    } else {
        if (any(duplicated(sampleNames))) {
            cli::cli_abort("Cannot set duplicated sample names.")
        }
    }

    for (i in seq_along(path)) {
        if (isTRUE(verbose)) {
            name <- sampleNames[i]
            if (is.null(name)) name <- paste0("sample ", i)
            cliID <- cli::cli_process_start("Reading from {.val {name}}")
        }
        if (is.list(path)) run <- path[[i]]
        else run <- path[i]

        if (!dir.exists(run)) {
            cli::cli_abort("Directory provided does not exist: {.file {normalizePath(run, mustWork = FALSE)}}")
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
            cli::cli_abort("Barcode file is missing. Expecting {.file {barcode.loc}}")
        }
        if (!file.exists(feature.loc) ) {
            cli::cli_abort("Peak or feature file is missing. Expecting {.file {feature.loc}}")
        }
        if (!file.exists(matrix.loc)) {
            cli::cli_abort("Expression matrix file is missing. Expecting {.file {matrix.loc}}")
        }
        data <- read10XFiles(matrixPath = matrix.loc,
                             barcodesPath = barcode.loc,
                             featuresPath = feature.loc,
                             sampleName = if (isTRUE(addPrefix)) sampleNames[i] else NULL,
                             geneCol = geneCol, cellCol = cellCol,
                             isATAC = pipeline == "atac",
                             returnList = TRUE)
        if (pipeline == "arc" && !arcFeatureType %in% names(data)) {
            cli::cli_abort(
                c("No ATAC data retrieved from cellranger-arc pipeline. ",
                "Please see if the following available feature types match ",
                "with need and select one for `arcFeatureType`: {.val {names(data)}}")
            )
        }
        data <- switch(pipeline,
                       arc = data[[arcFeatureType]],
                       atac = data[[1]]
        )
        allData[[i]] <- data
        cli::cli_process_done(id = cliID)
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


#' Read 10X cellranger files (matrix, barcodes and features) into R session
#' @description
#' This function works for loading a single sample with specifying the paths to
#' the matrix.mtx, barcodes.tsv, and features.tsv files. This function is
#' internally used by \code{\link{read10X}} functions for loading individual
#' samples from cellranger output directory, while it can also be convenient
#' when out-of-standard files are presented (e.g. data downloaded from GEO).
#' @param matrixPath Character string, path to the matrix MTX file. Can be
#' gzipped.
#' @param barcodesPath Character string, path to the barcodes TSV file. Can be
#' gzipped.
#' @param featuresPath Character string, path to the features TSV file. Can be
#' gzipped.
#' @param sampleName Character string attached as a prefix to the cell barcodes
#' loaded from the barcodes file. Default \code{NULL} does not add any prefix.
#' Useful when users plan to merge multiple samples into one matrix and need
#' to avoid duplicated cell barcodes from different batches.
#' @param geneCol An integer indicating which column in the features file to
#' extract as the feature identifiers. Default \code{2}.
#' @param cellCol An integer indicating which column in the barcodes file to
#' extract as the cell identifiers. Default \code{1}.
#' @param isATAC Logical, whether the data is for ATAC-seq. Default
#' \code{FALSE}. If \code{TRUE}, feature identifiers will be generated by
#' combining the first three columns of the features file in the format of
#' "chr:start-end".
#' @param returnList Logical, used internally by wrapper functions. Whether to
#' force putting the loaded matrix in a list even if there's only one matrix.
#' Default \code{FALSE}.
#' @export
#' @return For a single-modal sample, a dgCMatrix object, or a list of one
#' dgCMatrix when \code{returnList = TRUE}. A list of multiple dgCMatrix objects
#' when multiple feature types are detected.
#' @examples
#' \dontrun{
#' matrix <- read10XFiles(
#'     matrixPath = "path/to/matrix.mtx.gz",
#'     barcodesPath = "path/to/barcodes.tsv.gz",
#'     featuresPath = "path/to/features.tsv.gz"
#' )
#' }
read10XFiles <- function(
        matrixPath,
        barcodesPath,
        featuresPath,
        sampleName = NULL,
        geneCol = 2,
        cellCol = 1,
        isATAC = FALSE,
        returnList = FALSE
) {
    # Matrix can be easily read
    data <- methods::as(Matrix::readMM(matrixPath), "CsparseMatrix")

    # Processing barcodes
    cell.barcodes <- utils::read.table(barcodesPath, header = FALSE,
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
            cli::cli_abort(
                c("{.var geneCol} was set to {.val {geneCol}} but {.file feature.tsv.gz} (or {.file genes.tsv}) only has {ncol(fetures.names)} columns.",
                  "i" = "Try setting {.var geneCol} to a value <= {ncol(feature.names)}.")
            )
        }
        if (any(is.na(feature.names[, geneCol]))) {
            cli::cli_alert_warning(
                "Some feature names are NA. Replacing NA names with ID from the opposite column requested"
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
            cli::cli_alert_warning(
                "10X data contains more than one type and is being returned as a list containing matrices of each type."
            )
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
        if (isTRUE(returnList)) data <- list(data)
    }
    return(data)
}

#' Read 10X HDF5 file
#' @rdname read10XH5
#' @export
#' @description
#' Read count matrix from 10X CellRanger HDF5 file. By default, \code{read10XH5}
#' load scRNA, scATAC or multimodal data into memory (\code{inMemory = TRUE}).
#' To use LIGER in delayed mode for handling large datasets, set
#' \code{inMemory = FALSE} to load the data as a \code{DelayedArray} object. The
#' delayed mode only supports scRNA data for now.
#' @param filename Character string, path to the HDF5 file.
#' @param inMemory Logical, whether to load the data into memory. Default
#' \code{TRUE}. \code{FALSE} loads the data as a \code{DelayedArray} object.
#' @param useNames Logical, whether to use gene names as row names. Default
#' \code{TRUE}. \code{FALSE} uses gene IDs instead.
#' @param featureMakeUniq Logical, whether to make gene names unique. Default
#' \code{TRUE}.
#' @return A sparse matrix when only using older CellRanger output HDF5 file or
#' when only one genome and one modality is detected. When multiple genomes are
#' available, will return a list for each genome. When using multimodal data,
#' each genome will be a list of matrices for each modality. The matrix will be
#' of dgCMatrix class when in memory, or a TENxMatrix object when in delayed
#' mode.
#' @examples
#' matrix <- read10XH5(
#'     filename = system.file("extdata/ctrl.h5", package = "rliger"),
#'     inMemory = TRUE
#' )
#' class(matrix) # Should show dgCMatrix
#' if (requireNamespace("HDF5Array", quietly = TRUE)) {
#'    matrix <- read10XH5(
#'       filename = system.file("extdata/ctrl.h5", package = "rliger"),
#'       inMemory = FALSE
#'    )
#'    print(class(matrix)) # Should show TENxMatrix
#' }
read10XH5 <- function(
        filename,
        inMemory = TRUE,
        useNames = TRUE,
        featureMakeUniq = TRUE
) {
    if (isTRUE(inMemory)) {
        read10XH5Mem(filename, useNames, featureMakeUniq)
    } else {
        read10XH5Delay(filename, useNames, featureMakeUniq)
    }
}

#' @rdname read10XH5
#' @export
read10XH5Mem <- function(
        filename,
        useNames = TRUE,
        featureMakeUniq = TRUE
) {
    # Adopted from Seurat. with minimum modification to adapt to LIGER coding style
    if (!is.logical(useNames) || length(useNames) != 1) cli::cli_abort("{.field useNames} must be a logical.")
    if (!file.exists(filename)) {
        cli::cli_abort("File not found at {.file {filename}}")
    }
    infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
    genomes <- names(infile)
    output <- list()
    if (hdf5r::existsGroup(infile, 'matrix')) {
        # cellranger version 3
        featureSlot <- ifelse(useNames, 'features/name', 'features/id')
    } else {
        featureSlot <- ifelse(useNames, 'gene_names', 'genes')
    }
    for (genome in genomes) {
        # Modern cellranger output H5 really got genomes as the root group
        # Older ones just got "matrix" which is the "genome" here
        counts <- infile[[paste0(genome, '/data')]]
        indices <- infile[[paste0(genome, '/indices')]]
        indptr <- infile[[paste0(genome, '/indptr')]]
        shp <- infile[[paste0(genome, '/shape')]]
        features <- infile[[paste0(genome, '/', featureSlot)]][]
        barcodes <- infile[[paste0(genome, '/barcodes')]]
        mat <- Matrix::sparseMatrix(
            i = indices[] + 1,
            p = indptr[],
            x = as.numeric(counts[]),
            dims = shp[],
            repr = "T"
        )
        if (isTRUE(featureMakeUniq)) {
            features <- make.unique(features)
        }
        rownames(mat) <- features
        colnames(mat) <- barcodes[]
        mat <- methods::as(mat, "CsparseMatrix")
        # Split v3 multimodal
        if (infile$exists(paste0(genome, '/features/feature_type'))) {
            types <- infile[[paste0(genome, '/features/feature_type')]][]
            typesUniq <- unique(types)
            if (length(typesUniq) > 1) {
                cli::cli_alert_info(
                    "Genome {genome} has multiple modalities, returning a list of matrices for this genome"
                )
                mat <- sapply(
                    X = typesUniq,
                    FUN = function(x) {
                        return(mat[which(types == x), ])
                    },
                    simplify = FALSE,
                    USE.NAMES = TRUE
                )
            }
        }
        output[[genome]] <- mat
    }
    infile$close_all()
    if (length(output) == 1) {
        return(output[[genome]])
    } else{
        return(output)
    }
}

#' @rdname read10XH5
#' @export
read10XH5Delay <- function(
        filename,
        useNames = TRUE,
        featureMakeUniq = TRUE
) {
    if (!requireNamespace("HDF5Array", quietly = TRUE)) {
        cli::cli_abort(c(
            x = "Package {.pkg HDF5Array} is required for reading 10X H5 data into DelayedArray.",
            i = "Please install with {.code BiocManager::install('HDF5Array')}."
        ))
    }
    infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
    genomes <- names(infile)
    output <- list()
    if (hdf5r::existsGroup(infile, 'matrix')) {
        # cellranger version 3
        featureSlot <- ifelse(useNames, 'features/name', 'features/id')
    } else {
        featureSlot <- ifelse(useNames, 'gene_names', 'genes')
    }
    for (genome in genomes) {
        mat <- HDF5Array::TENxMatrix(filepath = filename, group = genome)
        features <- infile[[paste0(genome, '/', featureSlot)]][]
        if (isTRUE(featureMakeUniq)) {
            features <- make.unique(features)
        }
        mat@seed@dimnames[[1]] <- features
        if (infile$exists(paste0(genome, '/features/feature_type'))) {
            types <- infile[[paste0(genome, '/features/feature_type')]][]
            typesUniq <- unique(types)
            if (length(typesUniq) > 1) {
                cli::cli_abort(
                    "LIGER currently does not support multimodal data in DelayedArray."
                )
            }
        }
        output[[genome]] <- mat
    }
    infile$close_all()
    if (length(output) == 1) {
        return(output[[genome]])
    } else{
        return(output)
    }
}

#' Read matrix from H5AD file
#' @rdname readH5AD
#' @export
#' @description
#' Read raw count matrix from H5AD file. By default, \code{readH5AD} load
#' specified layer into memory (\code{inMemory = TRUE}). To use LIGER in delayed
#' mode for handling large datasets, set \code{inMemory = FALSE} to load the
#' data as a \code{DelayedArray} object. Note that only CSR format is supported
#' for the matrix.
#' @details
#' Currently, the only supported H5AD AnnData encoding versions are as follows:
#'
#' \itemize{
#'  \item{\code{adata.X}, \code{adata.raw.X}, or \code{adata.layers['layer']} -
#'  csr_matrix 0.1.0}
#'  \item{\code{adata.obs} and \code{adata.var} - dataframe 0.2.0}
#'  \item{Categoricals in a data frame - categorical 0.2.0}
#' }
#'
#' If users possess H5AD files encoded with older specification, please either
#' open an issue on GitHub or use R package 'anndata' to manually extract
#' information.
#' @param filename Character string, path to the H5AD file.
#' @param layer Character string specifying the H5 path of raw count data to be
#' loaded. Use \code{'X'} for \code{adata.X}, \code{'raw/X'} for
#' \code{adata.raw.X}, or \code{'layers/layer_name'} for
#' \code{adata.layers['layer_name']}.
#' @param inMemory Logical, whether to load the data into memory. Default
#' \code{TRUE}. \code{FALSE} loads the data as a \code{DelayedArray} object.
#' @param obs Logical, whether to also load the cell metadata from
#' \code{adata.obs}. Default \code{FALSE}.
#' @return When loaded in memory, a sparse matrix of class \code{dgCMatrix} will
#' be returned. When loaded in delayed mode, a \code{TENxMatrix} object will be
#' returned. If \code{obs = TRUE}, a list containing the matrix and the cell
#' metadata will be returned.
#' @examples
#' tempH5AD <- tempfile(fileext = '.h5ad')
#' writeH5AD(pbmc, tempH5AD, overwrite = TRUE)
#' mat <- readH5AD(tempH5AD, layer = 'X')
#' delayMat <- readH5AD(tempH5AD, layer = 'X', inMemory = FALSE)
readH5AD <- function(
        filename,
        layer,
        inMemory = TRUE,
        obs = FALSE
) {
    if (isTRUE(inMemory)) readH5ADMem(filename, layer, obs)
    else readH5ADDelay(filename, layer, obs)
}

#' @rdname readH5AD
#' @export
readH5ADMem <- function(
        filename,
        layer,
        obs = FALSE
) {
    dfile <- hdf5r::H5File$new(filename, mode = 'r')
    layer <- .H5AD.checkLayer(dfile, layer)

    # Processing the matrix
    if (inherits(dfile[[layer]], 'H5Group')) {
        res <- switch(
            EXPR = hdf5r::h5attr(dfile[[layer]], 'encoding-type'),
            csr_matrix = .H5AD.readCSR(dfile, layer)
        )
    } else if (inherits(dfile[[layer]], 'H5D')) {
        res <- .H5AD.readDense(dfile, layer)
    }

    obs_index_col <- hdf5r::h5attr(dfile[['obs']], '_index')
    obs_names <- dfile[['obs']][[obs_index_col]][]
    var_index_col <- hdf5r::h5attr(dfile[['var']], '_index')
    var_names <- dfile[['var']][[var_index_col]][]
    dimnames(res) <- list(var_names, obs_names)

    if (isTRUE(obs)) {
        obs <- .H5AD.readObs(dfile)
        res <- list(
            matrix = res,
            obs = obs
        )
    }

    dfile$close_all()
    return(res)
}

.H5AD.checkLayer <- function(dfile, layer) {
    layerAvail <- character()
    if (dfile$exists('X')) layerAvail <- c('X')
    if (dfile$exists('raw')) {
        if (dfile$exists('raw/X')) {
            layerAvail <- c(layerAvail, 'raw/X')
        }
    }
    if (dfile$exists('layers')) {
        layerAvail <- c(layerAvail, paste0('layers/', names(dfile[['layers']])))
    }
    if (missing(layer)) {
        cli::cli_abort(
            c(x = '{.field layer} must be specified',
              i = 'Available options include: {.val {layerAvail}}')
        )
    }
    if (!layer %in% layerAvail) {
        cli::cli_abort(
            c(x = 'Cannot identify specified layer {.val {layer}}.',
              i = 'Valid options include: {.val {layerAvail}}')
        )
    }
    return(layer)
}

.H5AD.readCSR <- function(dfile, layer) {
    i <- dfile[[layer]][['indices']][] + 1
    p <- dfile[[layer]][['indptr']][]
    x <- dfile[[layer]][['data']][]
    dims <- rev(hdf5r::h5attr(dfile[[layer]], 'shape'))
    Matrix::sparseMatrix(i = i, p = p, x = x, dims = dims, repr = 'C')
}

.H5AD.readDense <- function(dfile, layer) {
    chunkSize <- getOption('ligerChunkSize', 1000L)
    nChunk <- ceiling(dfile[[layer]]$dims[2] / chunkSize)
    matList <- list()
    cli::cli_progress_bar("Loading dense data by chunk", total = nChunk)
    for (i in seq_len(nChunk)) {
        start <- (i - 1) * chunkSize + 1
        end <- min(i * chunkSize, dfile[[layer]]$dims[2])
        chunkMat <- dfile[[layer]][, start:end]
        chunkMat <- methods::as(chunkMat, 'CsparseMatrix')
        matList[[i]] <- chunkMat
        cli::cli_progress_update()
    }
    return(Reduce(cbind, matList))
}

.H5AD.readObs <- function(dfile) {
    obs <- dfile[['obs']]
    obs_index_col <- hdf5r::h5attr(obs, '_index')
    obs_names <- obs[[obs_index_col]][]
    df <- data.frame(row.names = obs_names)
    for (i in seq_along(names(obs))) {
        cn <- names(obs)[i]
        if (cn == '_index') next
        if (cn == '__categories') next
        if (inherits(obs[[cn]], 'H5D')) {
            if ('categories' %in% hdf5r::h5attr_names(obs[[cn]])) {
                # Old stupid encoding
                categories <- hdf5r::h5attr(obs[[cn]], 'categories')$dereference()[[1]][]
                ints <- obs[[cn]][] + 1
                df[[cn]] <- factor(categories[ints], levels = categories)
            } else {
                df[[cn]] <- obs[[cn]][]
            }
        }
        if (inherits(obs[[cn]], 'H5Group')) {
            # Latest simple encoding
            categories <- obs[[cn]][['categories']]
            codes <- obs[[cn]][['codes']]
            fct <- factor(categories[codes[] + 1], levels = categories[])
            df[[cn]] <- fct
        }
    }
    return(df)
}

#' @rdname readH5AD
#' @export
readH5ADDelay <- function(
        filename,
        layer,
        obs = FALSE
) {
    if (!requireNamespace("HDF5Array", quietly = TRUE)) {
        cli::cli_abort(c(
            x = "Package {.pkg HDF5Array} is required for reading 10X H5 data into DelayedArray.",
            i = "Please install with {.code BiocManager::install('HDF5Array')}."
        ))
    }
    dfile <- hdf5r::H5File$new(filename, mode = 'r')
    layer <- .H5AD.checkLayer(dfile, layer)
    obs_index_col <- hdf5r::h5attr(dfile[['obs']], '_index')
    obs_names <- dfile[['obs']][[obs_index_col]][]
    var_index_col <- hdf5r::h5attr(dfile[['var']], '_index')
    var_names <- dfile[['var']][[var_index_col]][]

    res <- HDF5Array::TENxMatrix(filepath = filename, group = layer)
    res@seed@dimnames <- list(var_names, obs_names)

    if (isTRUE(obs)) {
        obs <- .H5AD.readObs(dfile)
        res <- list(
            matrix = res,
            obs = obs
        )
    }
    dfile$close_all()
    return(res)
}

# nocov end

