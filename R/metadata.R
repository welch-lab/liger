#' Liger cell/feature metadata
#' @rdname liger-metadata
#' @description
#' The rationale for defining these classes is to have a data frame extension
#' class holding the data, so that printed view looks tidy. \code{tibble} is
#' selected as the base class for metadata. However, \code{tibble} does not
#' allow \code{rownames}. Therefore, we lightly extend the class to preserve
#' a column for unique IDs of cells or features, and bring \code{rownames}
#' and character subscription of rows (`[]`) back to it.
#' @param x A data frame or tibble containing cell or feature metadata.
#' @param row.names If \code{x} is a tibble without the preserved unique ID
#' column (i.e., \code{.cellID} for cell metadata or \code{.featureID} for
#' feature metadata), which column should be used as the unique ID for cells or
#' features. This argument is ignored if \code{x} already has the preserved
#' unique ID column.
#' @param i,j Row and column subscription.
#' @param ... Additional arguments passed to `[.data.frame` method.
#' @param drop If \code{TRUE} the result is coerced to the lowest possible
#' dimension (see the examples). This only works for extracting elements, not
#' for the replacement.
#' @param value Replacement value when using the setter method.
#' @return
#' For \code{as.*} functions, a tibble with class \code{cellMeta} or
#' \code{featureMeta}. The first column is the preserved unique ID column (i.e.,
#' \code{.cellID} or \code{.featureID}).
#'
#' For \code{dimnames} method, a list of two character vectors: the first vector
#' is the unique IDs, and the second is the column names.
#'
#' For \code{[]} method, a subset of the original metadata tibble.
#' @export
#' @examples
#' as.cellMeta(pbmc@cellMeta)
as.cellMeta <- function(x, row.names = NULL) {
    if (is.null(x)) return(NULL)
    if (tibble::is_tibble(x)) {
        if ('.cellID' %in% colnames(x)) {
            x <- x %>%
                dplyr::relocate(.data[['.cellID']], .before = 1)
        } else {
            if (ncol(x) == 0) {
                nr <- nrow(x)
                x <- tibble::tibble(.cellID = character(nr))
            } else {
                if (is.null(row.names)) {
                    cli::cli_abort('{.cls tbl_df} input must have {.field .cellID} column, or a cell ID column must be specified for {.arg row.names}.')
                }
                x <- x %>%
                    dplyr::mutate(.cellID = .data[[row.names]]) %>%
                    dplyr::relocate(.data[['.cellID']], .before = 1)
            }
        }
    } else {
        x <- x %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = '.cellID') %>%
            tibble::as_tibble()
    }
    class(x) <- c('cellMeta', setdiff(class(x), 'cellMeta'))
    return(x)
}

#' @rdname liger-metadata
#' @export
as.featureMeta <- function(x, row.names = NULL) {
    if (is.null(x)) return(NULL)
    if (tibble::is_tibble(x)) {
        if ('.featureID' %in% colnames(x)) {
            x <- x %>%
                dplyr::relocate(.data[['.featureID']], .before = 1)
        } else {
            if (ncol(x) == 0) {
                nr <- nrow(x)
                x <- tibble::tibble(.featureID = character(nr))
            } else {
                if (is.null(row.names)) {
                    cli::cli_abort('{.cls tbl_df} input must have {.field .featureID} column, or a feature ID column must be specified for {.arg row.names}.')
                }
                x <- x %>%
                    dplyr::mutate(.featureID = .data[[row.names]]) %>%
                    dplyr::relocate(.data[['.featureID']], .before = 1)
            }
        }
    } else {
        x <- x %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = '.featureID') %>%
            tibble::as_tibble()
    }
    class(x) <- c('featureMeta', setdiff(class(x), 'featureMeta'))
    return(x)
}

#' @rdname liger-metadata
#' @export
#' @method dimnames cellMeta
dimnames.cellMeta <- function(x) {
    return(list(
        x$.cellID,
        colnames(x)
    ))
}

#' @rdname liger-metadata
#' @export
#' @method row.names<- cellMeta
`row.names<-.cellMeta` <- function(x, value) {
    x$.cellID <- value
    x
}

#' @rdname liger-metadata
#' @export
#' @method row.names<- featureMeta
`row.names<-.featureMeta` <- function(x, value) {
    x$.featureID <- value
    x
}

#' @rdname liger-metadata
#' @export
#' @method dimnames featureMeta
dimnames.featureMeta <- function(x) {
    return(list(
        x$.featureID,
        colnames(x)
    ))
}

#' @rdname liger-metadata
#' @export
`[.cellMeta` <- function(x, i, j, drop = FALSE) {
    if (!missing(i)) {
        if (is.character(i)) {
            i <- match(i, x$.cellID)
            if (any(is.na(i))) {
                cli::cli_abort('Cell ID subscription contains {sum(is.na(i))} ID{?s} not found in {.field .cellID} column.')
            }
        }
    }
    NextMethod()
}

#' @rdname liger-metadata
#' @export
`[.featureMeta` <- function(x, i, j, drop = FALSE) {
    if (!missing(i)) {
        if (is.character(i)) {
            i <- match(i, x$.featureID)
            if (any(is.na(i))) {
                cli::cli_abort('Feature ID subscription contains {sum(is.na(i))} ID{?s} not found in {.field .featureID} column.')
            }
        }
    }
    NextMethod()
}

#' @rdname liger-metadata
#' @export
`[<-.cellMeta` <- function(x, i, j, ..., value) {
    if (!missing(i)) {
        if (is.character(i)) {
            i <- match(i, x$.cellID)
            if (any(is.na(i))) {
                cli::cli_abort('Cell ID subscription contains {sum(is.na(i))} ID{?s} not found in {.field .cellID} column.')
            }
        }
    }
    NextMethod()
}

#' @rdname liger-metadata
#' @export
`[<-.featureMeta` <- function(x, i, j, ..., value) {
    if (!missing(i)) {
        if (is.character(i)) {
            i <- match(i, x$.featureID)
            if (any(is.na(i))) {
                cli::cli_abort('Feature ID subscription contains {sum(is.na(i))} ID{?s} not found in {.field .featureID} column.')
            }
        }
    }
    NextMethod()
}
