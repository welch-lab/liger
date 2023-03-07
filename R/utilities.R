# utility function for seuratToLiger function
# Compares colnames in reference matrix1 and adds back any missing
# column names to matrix.subset as rows
# Set transpose = TRUE if rownames of matrix1 should be referenced
addMissingCells <- function(matrix1, matrix.subset, transpose = F) {
  if (transpose) {
    matrix1 <- t(matrix1)
  }
  if (ncol(matrix1) != nrow(matrix.subset)) {
    extra <- matrix(
      NA,
      nrow = ncol(matrix1) - nrow(matrix.subset),
      ncol = ncol(matrix.subset)
    )
    colnames(extra) <- colnames(matrix.subset)
    rownames(extra) <-
      setdiff(colnames(matrix1), rownames(matrix.subset))
    matrix.subset <- rbind(matrix.subset, extra)
    # get original order
    matrix.subset <- matrix.subset[colnames(matrix1),]
  }
  return(matrix.subset)
}

# plyr::mapvalues
mapvalues <- function(x, from, to, warn_missing = TRUE) {
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x) && !is.null(x)) {
    stop("`x` must be an atomic vector or NULL.")
  }

  if (is.factor(x)) {
    # If x is a factor, call self but operate on the levels
    levels(x) <- mapvalues(levels(x), from, to, warn_missing)
    return(x)
  }

  mapidx <- match(x, from)
  mapidxNA  <- is.na(mapidx)

  # index of items in `from` that were found in `x`
  from_found <- sort(unique(mapidx))
  if (warn_missing && length(from_found) != length(from)) {
    message("The following `from` values were not present in `x`: ",
            paste(from[!(1:length(from) %in% from_found) ], collapse = ", "))
  }

  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}
