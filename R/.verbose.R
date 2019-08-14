.verbose <- function(expr) {
  invisible(utils::capture.output(expr))
}
