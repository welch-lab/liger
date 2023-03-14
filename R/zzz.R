ligerOptions <- list(
  ligerBaseSize = 10,
  ligerVerbose = TRUE
)
.onLoad <- function(libname, pkgname) {
  toset <- setdiff(names(ligerOptions), names(options()))
  if (length(toset) > 0) {
    options(ligerOptions[toset])
  }
  return(invisible(NULL))
}
