#' .checkFilePath
#'
#' @param path
#' @keywords internal
.checkFilePath <- function(path) {
  gsub("/*$", "/", path)
}
