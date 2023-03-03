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
