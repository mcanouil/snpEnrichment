#' .readSignal
#'
#' @param pattern
#' @keywords internal
.readSignal <- function(pattern) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  utils::read.delim(
    file = paste0(tmpDir, "/snpEnrichment/", pattern, ".signal"), header = TRUE, stringsAsFactors = FALSE,
    colClasses = c("NULL", "character", "numeric"), na.string = c("NA", ""),
    check.names = FALSE, strip.white = TRUE, col.names = c("", "SNP", "PVALUE")
  )
}
