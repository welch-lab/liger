setClassUnion("dgCMatrix_OR_NULL", c("dgCMatrix", "NULL"))
setClassUnion("matrix_OR_NULL", c("matrix", "NULL"))
setClassUnion("matrixLike", c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))
setClassUnion("matrixLike_OR_NULL", c("matrixLike", "NULL"))
setClassUnion("character_OR_NULL", c("character", "NULL"))
# It is quite hard to handle "H5D here, which is indeed defined as an R6 class.
# I'm not sure if this is a proper solution
setOldClass("H5D")
setOldClass("H5Group")
suppressWarnings(setClassUnion("dgCMatrix_OR_H5D_OR_NULL", c("dgCMatrix", "H5D", "NULL")))
setClassUnion("matrix_OR_H5D_OR_NULL", c("matrix", "H5D", "NULL"))
setClassUnion("matrixLike_OR_H5D_OR_H5Group_OR_NULL", c("matrixLike", "H5D", "H5Group", "NULL"))
setClassUnion("index",
              members = c("logical", "numeric", "character"))
setClassUnion("Number_or_NULL", c("integer", "numeric", "NULL"))
setClassUnion("dataframe", c("data.frame", "DataFrame", "NULL", "missing"))
setClassUnion("missing_OR_NULL", c("missing", "NULL"))

#' @importClassesFrom Matrix dgCMatrix dgTMatrix dgeMatrix
NULL


#' ligerDataset class
#'
#' Object for storing dastaset specific information. Will be embedded within a
#' higher level \linkS4class{liger} object
#' @docType class
#' @rdname ligerDataset-class
#' @slot rawData Raw data. Feature by cell matrix. Most of the time, sparse
#' matrix of integer numbers for RNA and ATAC data.
#' @slot normData Normalized data. Feature by cell matrix. Sparse if the
#' \code{rawData} it is normalized from is sparse.
#' @slot scaleData Scaled data, usually with subset shared variable features, by
#' cells. Most of the time sparse matrix of float numbers. This is the data used
#' for iNMF factorization.
#' @slot scaleUnsharedData Scaled data of variable features not shared with
#' other datasets. This is the data used for UINMF factorization.
#' @slot varUnsharedFeatures Variable features not shared with other datasets.
#' @slot V iNMF output matrix holding the dataset specific gene loading of each
#' factor. Feature by factor matrix.
#' @slot A Online iNMF intermediate product matrix.
#' @slot B Online iNMF intermediate product matrix.
#' @slot H iNMF output matrix holding the factor loading of each cell. Factor by
#' cell matrix.
#' @slot U UINMF output matrix holding the unshared variable gene loading of
#' each factor. Feature by factor matrix.
#' @slot h5fileInfo list of meta information of HDF5 file used for constructing
#' the object.
#' @slot featureMeta Feature metadata, DataFrame object.
#' @slot colnames Character vector of unique cell identifiers.
#' @slot rownames Character vector of unique feature names.
#' @importClassesFrom S4Vectors DataFrame
#' @exportClass ligerDataset
ligerDataset <- setClass(
    "ligerDataset",
    representation(
        rawData = "dgCMatrix_OR_H5D_OR_NULL",
        normData = "dgCMatrix_OR_H5D_OR_NULL",
        scaleData = "matrixLike_OR_H5D_OR_H5Group_OR_NULL",
        scaleUnsharedData = "matrixLike_OR_H5D_OR_H5Group_OR_NULL",
        varUnsharedFeatures = "character",
        H = "matrix_OR_NULL",
        V = "matrix_OR_NULL",
        A = "matrix_OR_NULL",
        B = "matrix_OR_NULL",
        U = "matrix_OR_NULL",
        h5fileInfo = "list",
        featureMeta = "DataFrame",
        colnames = "character",
        rownames = "character"
    )
)

.checkLigerDatasetBarcodes <- function(x) {
    # cell barcodes all consistant
    if (is.null(colnames(x))) {
        return(paste0("No valid cell barcode detected for ligerDataset.\n",
                      "Please create object with matrices with colnames."))
    }
    for (slot in c("rawData", "normData", "scaleData", "scaleUnsharedData",
                   "H")) {
        if (!slot %in% methods::slotNames(x)) next
        data <- methods::slot(x, slot)
        if (!is.null(data)) {
            barcodes.slot <- colnames(data)
            if (!identical(colnames(x), barcodes.slot)) {
                return(paste0("Inconsistant cell identifiers in `", slot,
                              "` slot."))
            }
        }
    }

    for (slot in c("scaleData", "V")) {
        featuresToCheck <- rownames(methods::slot(x, slot))
        check <- !featuresToCheck %in% rownames(x)
        if (any(check)) {
            msg <- paste0("Features in ", slot, " not found from dataset: ",
                          paste(featuresToCheck[check], collapse = ", "))
            return(msg)
        }
    }
    TRUE
}

.checkH5LigerDatasetLink <- function(x) {
    restoreGuide <- "Please try running `restoreH5Liger(object)`."
    if (!"H5File" %in% names(h5fileInfo(x))) {
        return(paste("`h5fileInfo` incomplete.", restoreGuide))
    }
    h5file <- getH5File(x)
    if (is.null(h5file)) {
        return(paste("`H5File` is NULL in `h5fileInfo` slot.", restoreGuide))
    }
    if (!h5file$is_valid) {
        return(paste("`H5File` is invalid in `h5fileInfo` slot.", restoreGuide))
    }
    if (!is.null(rawData(x))) {
        if (!rawData(x)$is_valid) {
            return(paste("`rawData` slot is invalid.", restoreGuide))
        }
    }
    if (!is.null(normData(x))) {
        if (!normData(x)$is_valid) {
            return(paste("`normData` slot is invalid.", restoreGuide))
        }
    }
    if (!is.null(scaleData(x))) {
        if (!scaleData(x)$is_valid) {
            return(paste("`scaleData` slot is invalid.", restoreGuide))
        }
    }
    TRUE
}

.valid.ligerDataset <- function(object) {
    if (isH5Liger(object)) {
        # message("Checking h5 ligerDataset validity")
        .checkH5LigerDatasetLink(object)
    } else {
        # message("Checking in memory ligerDataset validity")
        .checkLigerDatasetBarcodes(object)
    }
    # TODO more checks
    # TODO debating on whether to have check of the matching between scaleData
    # features and selected variable features.
}

setValidity("ligerDataset", .valid.ligerDataset)

#' @title liger class
#' @rdname liger-class
#' @docType class
#' @description \code{liger} object is the main data container for LIGER
#' analysis in R. The slot \code{datasets} is a list where each element should
#' be a \linkS4class{ligerDataset} object containing dataset specific
#' information, such as the expression matrices. The other parts of liger object
#' stores information that can be shared across the analysis, such as the cell
#' metadata.
#'
#' This manual provides explanation to the \code{liger} object structure as well
#' as usage of class-specific methods. Please see detail sections for more
#' information.
#'
#' For \code{liger} objects created with older versions of rliger package,
#' please try updating the objects individually with
#' \code{\link{convertOldLiger}}.
#' @slot datasets list of \linkS4class{ligerDataset} objects. Use generic
#' \code{dataset}, \code{dataset<-}, \code{datasets} or \code{datasets<-} to
#' interact with. See detailed section accordingly.
#' @slot cellMeta \link[S4Vectors]{DFrame} object for cell metadata. Pre-existing
#' metadata, QC metrics, cluster labeling and etc. are all stored here. Use
#' generic \code{cellMeta}, \code{cellMeta<-}, \code{$}, \code{[[]]} or
#' \code{[[]]<-} to interact with. See detailed section accordingly.
#' @slot varFeatures Character vector of names of variable features. Use generic
#' \code{varFeatures} or \code{varFeatures<-} to interact with. See detailed
#' section accordingly.
#' @slot W iNMF output matrix of shared gene loadings for each factor. See
#' \code{\link{runIntegration}}.
#' @slot H.norm Matrix of aligned factor loading for each cell. See
#' \code{\link{alignFactors}} and \code{\link{runIntegration}}.
#' @slot commands List of \linkS4class{ligerCommand} objects. Record of
#' analysis. Use \code{commands} to retrieve information. See detailed section
#' accordingly.
#' @slot uns List for unstructured meta-info of analyses or presets.
#' @slot version Record of version of rliger package
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom ggplot2 fortify
liger <- setClass(
    "liger",
    representation(
        datasets = "list",
        cellMeta = "DataFrame",
        varFeatures = "character_OR_NULL",
        W = "matrix_OR_NULL",
        H.norm = "matrix_OR_NULL",
        dimReds = "list",
        uns = "list",
        commands = "list",
        version = "ANY"
    ),
    methods::prototype(
        cellMeta = methods::new("DFrame"),
        version = utils::packageVersion("rliger")
    )
)


.checkAllDatasets <- function(x) {
    for (ld in datasets(x)) {
        methods::validObject(ld)
    }
    return(NULL)
}

.checkDatasetVar <- function(x) {
    if (!"dataset" %in% names(cellMeta(x))) {
        return("`datasets` variable missing in cellMeta(x)")
    }
    if (!is.factor(x@cellMeta$dataset)) {
        return("\"dataset\" variable in cellMeta is not a factor")
    }
    ds_cm <- as.character(levels(droplevels(x@cellMeta[["dataset"]])))
    ds_ds <- as.character(names(x@datasets))
    if (!identical(ds_cm, ds_ds)) {
        return("`levels(x$dataset)` does not match `names(x)`.")
    }

    datasetNamesFromDatasets <- as.character(rep(names(x), lapply(datasets(x), ncol)))
    names(datasetNamesFromDatasets) <- NULL
    if (!identical(datasetNamesFromDatasets, as.character(x$dataset))) {
        return("names of datasets do not match \"datasets\" variable in cellMeta")
    }
    return(NULL)
}

.checkDimReds <- function(x) {
    barcodes <- rownames(x@cellMeta)
    for (i in seq_along(x@dimReds)) {
        dr <- x@dimReds[[i]]
        drName <- names(x@dimReds[i])
        if (is.null(drName))
            return(paste("Unnamed dimReds at index", i))
        if (!inherits(dr, "matrix"))
            return(paste("DimReds", drName, "is not of matrix class"))
        if (!identical(rownames(dr), barcodes))
            return(paste("DimReds", drName, "does not match barcodes"))
    }
    return(NULL)
}

.checkLigerBarcodes <- function(x) {
    bcFromDatasets <- unlist(lapply(datasets(x), colnames), use.names = FALSE)
    if (!identical(colnames(x), bcFromDatasets)) {
        return("liger object barcodes do not match to barcodes in datasets")
    }
    if (!is.null(x@H.norm)) {
        if (!identical(rownames(x@H.norm), bcFromDatasets)) {
            return("H.norm barcodes do not match to barcodes in datasets.")
        }
    }

    return(NULL)
}

.checkLigerVarFeature <- function(x) {
    if (!is.null(varFeatures(x)) &&
        length(varFeatures(x)) > 0) {
        # if (!is.null(x@W))
        #     if (!identical(rownames(x@W), varFeatures(x)))
        #         return("Variable features do not match dimension of W matrix")
        for (d in names(x)) {
            ld <- dataset(x, d)
            # if (!is.null(ld@V)) {
            #     if (!identical(rownames(ld@V), varFeatures(x)))
            #         return(paste("Variable features do not match dimension",
            #                      "of V matrix in dataset", d))
            # }
            if (any(!varFeatures(x) %in% rownames(ld))) {
                nf <- setdiff(varFeatures(x), rownames(ld))
                return(cli::format_error(
                    "{length(nf)} variable feature{?s} do not exist in dataset {.val {d}}: {.val {nf}}"
                ))
            }
            # if (!is.null(scaleData(ld))) {
            #     if (!isH5Liger(ld)) {
            #         if (!identical(rownames(scaleData(ld)), varFeatures(x)))
            #             return(paste("Variable features do not match dimension",
            #                          "of scaleData in dataset", d))
            #     } else {
            #         if (inherits(scaleData(ld), "H5D")) {
            #             if (scaleData(ld)$dims[1] != length(varFeatures(x)))
            #                 return(paste("Variable features do not match ",
            #                              "dimension of scaleData in dataset ",
            #                              "(H5)", d))
            #         } else if (inherits(scaleData(ld), "H5Group")) {
            #             if (scaleData(ld)[["featureIdx"]]$dims != length(varFeatures(x))) {
            #                 return(paste("Variable features do not match ",
            #                              "dimension of scaleData in dataset ",
            #                              "(H5)", d))
            #             }
            #             scaleDataIdx <- scaleData(ld)[["featureIdx"]][]
            #             if (!identical(rownames(ld)[scaleDataIdx], varFeatures(x))) {
            #                 return("HDF5 scaled data feature index does not ",
            #                        "match variable features")
            #             }
            #         }
            #     }
            # }
        }
    }
    return(NULL)
}

.valid.liger <- function(object) {
    res <- .checkAllDatasets(object)
    if (!is.null(res)) return(res)
    res <- .checkDimReds(object)
    if (!is.null(res)) return(res)
    res <- .checkDatasetVar(object)
    if (!is.null(res)) return(res)
    res <- .checkLigerBarcodes(object)
    if (!is.null(res)) return(res)
    res <- .checkLigerVarFeature(object)
    if (!is.null(res)) return(res)
    # TODO more checks
}

setValidity("liger", .valid.liger)



setClassUnion("POSIXct_or_NULL", c("POSIXct", "NULL"))

#' ligerCommand object: Record the input and time of a LIGER function call
#' @slot funcName Name of the function
#' @slot time A time stamp object
#' @slot call A character string converted from system call
#' @slot parameters List of all arguments except the \linkS4class{liger} object.
#' Large object are summarized to short string.
#' @slot objSummary List of attributes of the \linkS4class{liger} object as a
#' snapshot when command is operated.
#' @slot ligerVersion Character string converted from
#' \code{packageVersion("rliger")}.
#' @slot dependencyVersion Named character vector of version number, if any
#' dependency library has a chance to be included by the function. A
#' dependency might only be invoked under certain conditions, such as using
#' an alternative algorithm, which a call does not actually reach to, but it
#' would still be included for this call.
#' @exportClass ligerCommand
#' @export
#' @rdname ligerCommand-class
ligerCommand <- setClass(
    Class = "ligerCommand",
    representation(
        funcName = "character",
        time = "POSIXct_or_NULL",
        call = "character",
        parameters = "list",
        objSummary = "list",
        ligerVersion = "character",
        dependencyVersion = "character"
    ),
    prototype(
        funcName = character(),
        time = NULL,
        parameters = list(),
        objSummary = list(
            datasets = character(),
            nCells = numeric(),
            nFeatures = numeric(),
            nVarFeatures = numeric(),
            cellMetaNames = character(),
            ligerVersion = character(),
            dependencyVersion = character()
        )
    )
)


################################################################################
# Developer guide for adding a new sub-class of `ligerDataset` for new modality
################################################################################
#
# Below is a check-list of the TODOs when new sub-classes need to be added.
# Please follow them carefully, and refer to existing code as examples.
#
# 1. Add `setClass` chunk for defining the new subclass, pay attention to:
#   a. Naming convention should be `liger{Modal}Dataset`, in camelCase
#   b. contains = "ligerDataset"
#   c. add new slots for modality specific information with `representation`
#   d. if the default new information could be empty, add `prototype`
# 2. In files `zzz.R`, `import.R`, `classConversion.R`, search for
#    text "modal". When seeing a multi-option vector argument, add a unique
#    abbreviation of this new data type to the vector. Don't forget updating
#    valid options in the manual documentaion as well.
# 3. If the new slot(s) added is thought to be retrieved by future developers
#    or users, getter and setter methods MUST be implemented.
# 4. Please go through the implementation of the following functions in file
#    `ligerDataset-class.R`, and make sure data in the new slot(s) is properly
#    handled.
#   a. .checkLigerDatasetBarcodes()
#   b. `dimnames<-()` (search: `setReplaceMethod("dimnames"`)
#   c. `[` (search: "[")
#
################################################################################

#' Subclass of ligerDataset for RNA modality
#' @description Inherits from \linkS4class{ligerDataset} class. Contained slots
#' can be referred with the link. This subclass does not have any different from
#' the default \code{ligerDataset} class except the class name.
#' @export
#' @exportClass ligerRNADataset
ligerRNADataset <- setClass(
    "ligerRNADataset", contains = "ligerDataset"
)


#' Subclass of ligerDataset for ATAC modality
#'
#' @description Inherits from \linkS4class{ligerDataset} class. Contained slots
#' can be referred with the link.
#' @slot rawPeak sparse matrix
#' @slot normPeak sparse matrix
#' @exportClass ligerATACDataset
#' @export
ligerATACDataset <- setClass(
    "ligerATACDataset",
    contains = "ligerDataset",
    representation = representation(rawPeak = "matrixLike_OR_NULL",
                                    normPeak = "matrixLike_OR_NULL"),
    prototype = prototype(rawPeak = NULL, normPeak = NULL)
)

.valid.ligerATACDataset <- function(object) {
    passSuperClassCheck <- .valid.ligerDataset(object)
    if (!isTRUE(passSuperClassCheck)) return(passSuperClassCheck)
    for (slot in c("rawPeak", "normPeak")) {
        data <- methods::slot(object, slot)
        if (!is.null(data)) {
            barcodes.slot <- colnames(data)
            if (!identical(object@colnames, barcodes.slot)) {
                return(paste0("Inconsistant cell identifiers in `", slot,
                              "` slot."))
            }
        }
    }
}

setValidity("ligerATACDataset", .valid.ligerATACDataset)

#' Subclass of ligerDataset for Spatial modality
#'
#' @description Inherits from \linkS4class{ligerDataset} class. Contained slots
#' can be referred with the link.
#' @slot coordinate dense matrix
#' @exportClass ligerSpatialDataset
#' @export
ligerSpatialDataset <- setClass(
    "ligerSpatialDataset",
    contains = "ligerDataset",
    representation = representation(coordinate = "matrix_OR_NULL"),
    prototype = prototype(coordinate = NULL)
)

.checkCoords <- function(ld, value) {
    if (is.null(rownames(value))) {
        cli::cli_alert_warning("No rownames with given spatial coordinate. Assuming they match with the cells.")
        rownames(value) <- colnames(ld)
    }
    if (is.null(colnames(value))) {
        if (ncol(value) <= 3) {
            colnames(value) <- c("x", "y", "z")[seq(ncol(value))]
        } else {
            cli::cli_abort("More than 3 dimensions for the coordinates but no colnames are given.")
        }
        cli::cli_alert_warning(
            "No colnames with given spatial coordinate. Setting to {.val {colnames(value)}}"
        )
    }
    full <- matrix(NA, nrow = ncol(ld), ncol = ncol(value),
                   dimnames = list(colnames(ld), colnames(value)))
    cellIsec <- intersect(rownames(value), colnames(ld))
    full[cellIsec, colnames(value)] <- value[cellIsec,]
    if (any(is.na(full))) {
        cli::cli_alert_warning("NA generated for missing cells.")
    }
    if (any(!rownames(value) %in% rownames(full))) {
        cli::cli_alert_warning("Cells in given coordinate not found in the dataset.")
    }
    return(full)
}

.valid.ligerSpatialDataset <- function(object) {
    passSuperClassCheck <- .valid.ligerDataset(object)
    if (!isTRUE(passSuperClassCheck)) return(passSuperClassCheck)
    coord <- object@coordinate
    if (!is.null(coord)) {
        barcodes.slot <- rownames(coord)
        if (!identical(object@colnames, barcodes.slot)) {
            return(paste0("Inconsistant cell identifiers in `coordinate` slot."))
        }
    }
}

setValidity("ligerSpatialDataset", .valid.ligerSpatialDataset)



#' Subclass of ligerDataset for Methylation modality
#'
#' @description Inherits from \linkS4class{ligerDataset} class. Contained slots
#' can be referred with the link. \code{\link{scaleNotCenter}} applied on
#' datasets of this class will automatically be taken by reversing the
#' normalized data instead of scaling the variable features.
#' @exportClass ligerMethDataset
#' @export
ligerMethDataset <- setClass(
    "ligerMethDataset",
    contains = "ligerDataset"
)

.valid.ligerMethDataset <- function(object) .valid.ligerDataset(object)

setValidity("ligerMethDataset", .valid.ligerMethDataset)
