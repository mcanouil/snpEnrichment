#' GC
#'
#' \code{GC} performs garbage collection until free memory indicators show no
#' change.
#'
#'
#' @param verbose A \code{logical}. If \code{TRUE}, the garbage collection
#' prints statistics about cons cells and the space allocated for vectors.
#' @param reset A \code{logical}. If \code{TRUE} the values for maximum space
#' used are reset to the current values.
#' @keywords internal
GC <- function(verbose = getOption("verbose"), reset = FALSE) {
  while (!identical(gc(verbose, reset)[, 4], gc(verbose, reset)[, 4])) {}
  gc(verbose, reset)
}
