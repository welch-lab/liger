
#' @rdname peak
#' @export
setMethod("rawPeak", signature(x = "ligerATACDataset", dataset = "missing"),
          function(x, dataset = NULL) {
              x@rawPeak
          })

#' @rdname peak
#' @export
setReplaceMethod(
    "rawPeak",
    signature(x = "ligerATACDataset", dataset = "missing"),
    function(x, dataset = NULL, check = TRUE, value) {
        # It is very possible that the peak data inserted does not have the same
        # colnames as the ligerDataset object, where the latter could have a
        # prefix added when creating the liger object. We need to check for this
        # and adjust the colnames of the peak data accordingly.
        matching <- endsWith(colnames(x), colnames(value))
        if (!all(matching)) {
            cli::cli_abort(
                c("x" = "It seems that the cell identifiers from the inserted peak count do not (partially) match with those in the object.",
                  "i" = "The first three from the object: {.val {colnames(x)[1:3]}}",
                  "i" = "The first three from the peak: {.val {colnames(value)[1:3]}}")
            )
        }
        prefix <- gsub(colnames(value)[1], "", colnames(x)[1])
        if (nchar(prefix) > 0) {
            colnames(value) <- paste0(prefix, colnames(value))
        }
        x@rawPeak <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    })


#' @rdname peak
#' @export
setMethod("normPeak", signature(x = "ligerATACDataset", dataset = "missing"),
          function(x, dataset = NULL) {
              x@normPeak
          })

#' @rdname peak
#' @export
setReplaceMethod(
    "normPeak",
    signature(x = "ligerATACDataset", dataset = "missing"),
    function(x, dataset = NULL, check = TRUE, value) {
        matching <- endsWith(colnames(x), colnames(value))
        if (!all(matching)) {
            cli::cli_abort(
                c("x" = "It seems that the cell identifiers from the inserted peak count do not (partially) match with those in the object.",
                  "i" = "The first three from the object: {.val {colnames(x)[1:3]}}",
                  "i" = "The first three from the peak: {.val {colnames(value)[1:3]}}")
            )
        }
        prefix <- gsub(colnames(value)[1], "", colnames(x)[1])
        if (nchar(prefix) > 0) {
            colnames(value) <- paste0(prefix, colnames(value))
        }
        x@normPeak <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    })

#' @rdname coordinate
#' @export
setMethod("coordinate", signature(x = "ligerSpatialDataset", dataset = "missing"),
          function(x, dataset = NULL) {
              x@coordinate
          })

#' @rdname coordinate
#' @export
setReplaceMethod(
    "coordinate",
    signature(x = "ligerSpatialDataset", dataset = "missing"),
    function(x, dataset = NULL, check = TRUE, value) {
        value <- .checkCoords(ld = x, value = value)
        x@coordinate <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    })

