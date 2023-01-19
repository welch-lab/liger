.collapseLongNames <- function(x) {
    if (length(x) < 5) {
        return(paste(x, collapse = ", "))
    } else {
        head <- paste(x[seq(3)], collapse = ", ")
        tail <- x[length(x)]
        return(paste0(head, ", ..., ", tail))
    }
}

# Returns character names of datasets
.checkUseDatasets <- function(object, useDatasets = NULL) {
    if (!inherits(object, "liger"))
        stop("A liger object is required.")
    if (is.null(useDatasets)) useDatasets <- names(object)
    else {
        if (is.numeric(useDatasets)) {
            if (max(useDatasets) > length(object)) {
                stop("Numeric dataset index out of bound. Only ",
                     length(object), " datasets exist.")
            }
            useDatasets <- unique(useDatasets)
            useDatasets <- names(object)[useDatasets]
        } else if (is.logical(useDatasets)) {
            if (length(useDatasets) != length(object)) {
                stop("Logical dataset subscription does not match the number ",
                     "of datasets (", length(object), ").")
            }
            useDatasets <- names(object)[useDatasets]
        } else if (is.character(useDatasets)) {
            if (any(!useDatasets %in% names(object))) {
                notFound <- useDatasets[!useDatasets %in% names(object)]
                stop("Specified dataset name(s) not found: ",
                     paste(notFound, collapse = ", "))
            }
        } else {
            stop("Please use a proper numeric/logical/character vector to ",
                 "select dataset to use.")
        }
    }
    useDatasets
}

# Used for checking the subscription of cells or features and returns numeric
# subscription
.idxCheck <- function(object, idx, orient) {
    if (orient == "feature") {
        getNames <- rownames
        getNumber <- nrow
    } else if (orient == "cell") {
        getNames <- colnames
        getNumber <- ncol
    } else {
        stop("`orient` should be either 'feature' or 'cell'.")
    }
    paramName <- paste0("`", orient, "Idx`")
    if (is.character(idx)) {
        if (any(!idx %in% getNames(object))) {
            notFound <- paste(idx[!idx %in% getNames(object)],
                              collapse = ", ")
            stop(paramName, " not found in object: ", notFound)
        }
        idx <- which(getNames(object) %in% idx)
    } else if (is.logical(idx)) {
        if (length(idx) != getNumber(object)) {
            stop("Length of logical ", paramName, " does not match to ",
                 "number of ", orient, "s in `object`.")
        }
        idx <- which(idx)
    } else if (is.numeric(idx)) {
        if (max(idx) > getNumber(object)) {
            stop("Numeric ", paramName, " subscription out of bound.")
        }
    } else if (is.null(idx)) {
        idx <- seq(getNumber(object))
    } else {
        stop("Please use character, logical or numeric subscription ",
             "for ", paramName, ".")
    }
    return(idx)
}

.checkLDSlot <- function(object, slot) {
    if (!inherits(object, "ligerDataset"))
        stop("Please use a ligerDataset object.")
    avail <- c("raw.data", "norm.data", "scale.data")
    if (is.null(slot)) slot <- avail
    else {
        if (any(!slot %in% avail)) {
            notFound <- slot[!slot %in% avail]
            stop("Specified slot not availalble: ",
                 paste(notFound, collapse = ", "), ". Use one or more from ",
                 '"raw.data", "norm.data" or "scale.data"')
        }
        if ("raw.data" %in% slot &
            is.null(raw.data(object))) {
            stop("`raw.data` is not available for use.")
        }
        if ("norm.data" %in% slot &
            is.null(norm.data(object))) {
            stop("`norm.data` is not available for use.")
        }
        if ("scale.data" %in% slot &
            is.null(scale.data(object))) {
            stop("`scale.data` is not available for use.")
        }
    }
    slot
}

