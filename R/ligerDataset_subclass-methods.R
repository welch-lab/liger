
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

