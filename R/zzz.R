#' @importFrom Rcpp evalCpp
#' @importFrom Matrix colSums rowSums t summary
#' @importFrom rlang .data %||%
#' @importFrom methods new show
#' @importFrom utils .DollarNames
#' @useDynLib rliger2, .registration = TRUE
NULL

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

.modalClassDict <- list(
    default = "ligerDataset",
    rna = "ligerRNADataset",
    atac = "ligerATACDataset",
    spatial = "ligerSpatialDataset"
)

.classModalDict <- list(
    ligerDataset = "default",
    ligerRNADataset = "rna",
    ligerATACDataset = "atac",
    ligerSpatialDataset = "spatial"
)

.ligerOptions <- list(
    ligerBaseSize = 10,
    ligerVerbose = TRUE,
    ligerDotSize = 1
)

.onLoad <- function(libname, pkgname) {
    toset <- setdiff(names(.ligerOptions), names(options()))
    if (length(toset) > 0) {
        options(.ligerOptions[toset])
    }
    return(invisible(NULL))
}

.onAttach <- function(libname, pkgname) {
    msg <- paste0(
        "Package `rliger` has been updated massively since version 1.99.0, ",
        "including the object structure which is not compatible with old ",
        "versions.\nWe recommand you backup your old analysis before ",
        "overwriting any existing result.\n`readLiger()` is provided for ",
        "reading RDS file storing an old object and it converts it to the ",
        "up-to-date structure."
    )
    ggrepelPath <- find.package("ggrepel", quiet = TRUE)
    if (length(ggrepelPath) == 0) {
        msg <- paste0(msg, "\n\nPackage \"ggrepel\" is highly recommended to ",
                      "be installed for better plotting quality.")
    }
    packageStartupMessage(msg)
    return(invisible(NULL))
}

.onUnload <- function(libpath) {
    library.dynam.unload("rliger2", libpath)
    return(invisible(NULL))
}
