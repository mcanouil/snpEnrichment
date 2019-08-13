#' .enrichmentRatio
#'
#' @param table
#' @keywords internal
.enrichmentRatio <- function(table) {
  tmp <- as.numeric(table)
  (tmp[1] * tmp[4]) / (tmp[2] * tmp[3])
}
