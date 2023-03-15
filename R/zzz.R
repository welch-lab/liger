#' @importFrom Rcpp evalCpp
#' @importFrom Matrix colSums rowSums t summary
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom methods new show
NULL



.modalClassDict <- list(
    default = "ligerDataset",
    rna = "ligerRNADataset",
    atac = "ligerATACDataset"
)

.classModalDict <- list(
    ligerDataset = "default",
    ligerRNADataset = "rna",
    ligerATACDataset = "atac"
)

.ligerOptions <- list(
    ligerBaseSize = 10,
    ligerVerbose = TRUE
)

.onLoad <- function(libname, pkgname) {
    toset <- setdiff(names(.ligerOptions), names(options()))
    if (length(toset) > 0) {
        options(.ligerOptions[toset])
    }
    return(invisible(NULL))
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        "Package `rliger` has been updated massively since version 1.99.0, ",
        "including the object structure which is not compatible with old ",
        "versions.\nWe recommand you backup your old analysis before ",
        "overwriting any existing result.\n`readLiger()` is provided for ",
        "reading RDS file storing an old object and it converts it to the ",
        "up-to-date structure."
    )
    return(invisible(NULL))
}

.onUnload <- function(libpath) {
    library.dynam.unload("rliger", libpath)
    return(invisible(NULL))
}
