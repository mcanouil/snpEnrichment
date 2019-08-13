#' .verbose
#'
#' @param expr
#' @keywords internal
.verbose <- function(expr) {
  invisible(utils::capture.output(expr))
}
