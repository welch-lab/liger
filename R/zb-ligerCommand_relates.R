# Call in a way like:
# object <- recordCommand(object, dependencies = ...)
# Conditionally, should be placed after input checks
# like `match.arg()` or `.checkUseDatasets()`
# `...` is for the ... arguments in real function call, so S3 arguments passed
# to downstream can be also captured
recordCommand <- function(
        object,
        ...,
        dependencies = NULL
) {
    # Generate time stamp
    time <- Sys.time()
    # Capture the call
    if (as.character(x = sys.calls()[[1]])[1] == "do.call") {
        # Using do.call to run
        callString <- deparse(expr = sys.calls()[[1]])
        funcName <- as.character(x = sys.calls()[[1]])[2]
    } else if (as.character(sys.call(-1))[1] == ".local") {
        # Probably using S4 methods, e.g. quantileNorm
        call <- call <- match.call(definition = sys.function(-3),
                                   call = sys.call(-3))
        callString <- deparse(call)
        funcName <- as.character(as.list(call)[[1]])
    } else {
        # Normal function call
        #call <- match.call(definition = sys.function(-1),
        #                   call = sys.call(-1), expand.dots = FALSE)
        callString <- deparse(sys.call(-1))
        funcName <- as.character(as.list(sys.call(-1))[[1]])
    }
    # Get all existing variable from the function environment
    # This would require that `recordCommand` has to be called at the very
    # beginning of working exported function but after `match.arg()` and etc.
    # is done.
    args <- as.list(parent.frame())
    # Capture more arguments
    moreArgs <- list(...)
    args <- c(args, moreArgs)
    args$object <- NULL
    objSummary <- list(
        datasets = names(object),
        nCells = sapply(datasets(object), ncol),
        nFeatures = sapply(datasets(object), nrow),
        nVarFeatures = length(varFeatures(object)),
        cellMetaNames = colnames(cellMeta(object)),
        dimW = dim(getMatrix(object, "W")),
        dimHNorm = dim(getMatrix(object, "H.norm"))
    )
    ## Replace other large object in `args` with text representation
    args <- lapply(args, function(x) {
        if (utils::object.size(x) <= 500) x
        else .makeTextRepr(x)
    })
    # Get dependencies' version recorded
    if (!is.null(dependencies))
        dependVer <- sapply(dependencies, function(x)
            as.character(utils::packageVersion(x))
        )
    else dependVer <- character()
    Command <- methods::new(
        "ligerCommand", funcName = funcName,
        objSummary = objSummary,
        parameters = args,
        ligerVersion = as.character(utils::packageVersion("rliger")),
        dependencyVersion = dependVer
    )
    # Do a hash tag to see if an identical operation has been done before
    logEntry <- paste0(funcName, "_", substr(rlang::hash(Command), 1, 10))
    # The following items are not included for a hash check
    Command@time <- time
    Command@call <- callString
    # If same command applied before, remove and append instead of
    # replacing in place, in order to reflect time order
    object@commands[[logEntry]] <- NULL
    object@commands[[logEntry]] <- Command
    return(object)
}

#' @param object A \code{ligerCommand} object
#' @export
#' @rdname ligerCommand-class
#' @examples
#' pbmc <- normalize(pbmc)
#' cmd <- commands(pbmc, "normalize")
#' cmd
setMethod(
    "show",
    signature(object = "ligerCommand"),
    function(object) {
        cat("A liger command record, performed at ")
        cat(format(object@time, "%m-%d-%Y %H:%M:%S", tz = "America/New_York",
                   usetz = TRUE), "\n")
        cat("Call:", object@call, "\n")
        cat("Parameters:\n")
        for (p in names(object@parameters)) {
            cat("   ", p, ":", .showParam(object@parameters[[p]]), "\n")
        }
        invisible(x = NULL)
    }
)

# Used only when need to have a text representation for large objects
.makeTextRepr <- function(x) {
    if (!is.null(dim(x))) {
        paste0("Large object of class `", class(x)[1], "`, size: ",
               paste(dim(x), collapse = " x "))
    } else if (is.numeric(x) || is.character(x) ||
               is.logical(x) || is.factor(x)) {
        paste0("Long ", class(x), " with ", length(x), " elements: ",
               .collapseLongNames(x))
    } else if (is.list(x)) {
        paste0("Large list object with ", length(x), " elements.")
    } else {
        paste0("Large object of class ", class(x))
    }
}

# Used when invoking "show" method of ligerCommand object
.showParam <- function(x) {
    if (is.null(x)) "NULL"
    else if (is.matrix(x))
        paste0("matrix(c(", paste0(x, collapse = ", "), "), nrow = ", nrow(x),
               ", ncol = ", ncol(x), ")")
    else if (is.numeric(x) | is.logical(x)) paste0(x, collapse = ", ")
    else if (is.character(x)) paste0("\"", x, "\"", collapse = ", ")
    else if (is.list(x)) paste0("list object with ", length(x), " elements.")
    else paste("Object of class", class(x))
}

#' Check difference of two liger command
#' @param object \linkS4class{liger} object
#' @param cmd1,cmd2 Exact string of command labels. Available options could be
#' viewed with running \code{commands(object)}.
#' @return If any difference found, a character vector summarizing all
#' differences
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- normalize(pbmc, log = TRUE, scaleFactor = 1e4)
#' cmds <- commands(pbmc)
#' commandDiff(pbmc, cmds[1], cmds[2])
commandDiff <- function(object, cmd1, cmd2) {
    cmd1 <- commands(object, cmd1)
    if (!inherits(cmd1, "ligerCommand"))
        cli::cli_abort("{.code cmd1} matching with multiple command records.
                       Availble options could be viewed with {.code commands(object)}.")
    cmd2 <- commands(object, cmd2)
    if (!inherits(cmd2, "ligerCommand"))
        cli::cli_abort("{.code cmd2} matching with multiple command records.
                       Availble options could be viewed with {.code commands(object)}.")
    .cmdDiff(cmd1, cmd2)
}

.cmdDiff <- function(x, y) {
    msg <- character()
    if (identical(x, y)) return(msg)
    if (x@funcName != y@funcName)
        return(paste0("Functions are different: ", x@funcName,
                      "() VS ", y@funcName, "()"))
    # When functions are the same, add detail
    # Object summary
    xObjSum <- x@objSummary
    yObjSum <- y@objSummary
    for (i in union(names(xObjSum), names(yObjSum))) {
        if (!i %in% names(xObjSum)) {
            msg <- c(msg, paste0("Entry not found in `x@objSummary`: ", i))
            next
        }
        if (!i %in% names(yObjSum)) {
            msg <- c(msg, paste0("Entry not found in `y@objSummary`: ", i))
            next
        }
        if (!identical(xObjSum[[i]], yObjSum[[i]])) {
            msg <- c(msg, paste0("Not identical: `x@objSummary$", i, "`: ",
                                 .collapseLongNames(xObjSum[[i]]),
                                 "; `y@objSummary$", i, "`: ",
                                 .collapseLongNames(yObjSum[[i]])))
        }
    }
    # Parameters
    xArg <- x@parameters
    yArg <- y@parameters
    for (i in union(names(xArg), names(yArg))) {
        if (i == "verbose") next
        if (!i %in% names(xArg)) {
            msg <- c(msg, paste0("Argument not found in `x`: ", i))
            next
        }
        if (!i %in% names(yArg)) {
            msg <- c(msg, paste0("Argument not found in `y`: ", i))
            next
        }
        if (!identical(xArg[[i]], yArg[[i]])) {
            msg <- c(msg, paste0("Argument not identical: ", i))
        }
    }
    # Versions
    if (x@ligerVersion != y@ligerVersion)
        msg <- c(msg, paste0("\"rliger\" versions differ: ", x@ligerVersion,
                             " VS ", y@ligerVersion))
    if (!identical(x@dependencyVersion, y@dependencyVersion))
        msg <- c(msg, paste0("Dependency versions differ. ",
                             "Please print `x@dependencyVersion` and ",
                             "`y@dependencyVersion` to see the differences."))
    return(msg)
}


