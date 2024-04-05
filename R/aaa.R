#' @importFrom Rcpp evalCpp
#' @importFrom Matrix colSums rowSums t summary
#' @importFrom rlang .data %||%
#' @importFrom methods new show
#' @importFrom utils .DollarNames
#' @useDynLib rliger, .registration = TRUE
NULL

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom magrittr %<>%
#' @export
magrittr::`%<>%`

scPalette <- c('#E41A1C', '#377EB8', '#4DAF4A', '#FFCF00', '#aa47b9', '#e67c14',
               '#e7a2b4', '#54B0E4', '#9a5831', '#BC9DCC', '#222F75', '#1B9E77',
               '#B2DF8A', '#E3BE00', '#FF6699', '#8f3c4d', '#01e1e6', '#591cc5',
               '#A6CEE3', '#CE1261', '#8CA77B', '#5E4FA2', '#08692e', '#DCF0B9',
               '#8DD3C7', '#AAAA77')

.modalClassDict <- list(
    default = "ligerDataset",
    rna = "ligerRNADataset",
    atac = "ligerATACDataset",
    spatial = "ligerSpatialDataset",
    meth = "ligerMethDataset"
)

.classModalDict <- list(
    ligerDataset = "default",
    ligerRNADataset = "rna",
    ligerATACDataset = "atac",
    ligerSpatialDataset = "spatial",
    ligerMethDataset = "meth"
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
        "versions.\n\n",
        "We recommand you backup your old analysis before overwriting any ",
        "existing result.\n\n",
        "`readLiger()` is provided for reading an RDS file storing an old ",
        "object and it converts the object to the up-to-date structure."
    )
    ggrepelPath <- find.package("ggrepel", quiet = TRUE)
    if (length(ggrepelPath) == 0) {
        msg <- paste0(msg, "\n\nPackage \"ggrepel\" is highly recommended to ",
                      "be installed for better plotting quality. Users can ",
                      "install it with `install.packages('ggrepel')`.")
    }
    rcppplancPath <- find.package("RcppPlanc", quiet = TRUE)
    if (length(rcppplancPath) == 0) {
        msg <- paste0(msg, "\n\nPackage \"RcppPlanc\" is required for ",
                      "performing integrative factorizations but is not ",
                      "detected. It is currently not on CRAN, so please ",
                      "install it with ",
                      "`devtools::install_github('welch-lab/RcppPlanc')`")
    }
    packageStartupMessage(msg)
    return(invisible(NULL))
}

.onUnload <- function(libpath) {
    library.dynam.unload("rliger", libpath)
    return(invisible(NULL))
}
