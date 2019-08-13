#' .checkSnpListDir
#'
#' @param pattern
#' @keywords internal
.checkSnpListDir <- function(snpListDir, pattern) {
  snpListDir <- .checkFilePath(snpListDir)
  if (length(list.files(snpListDir, pattern = pattern)) == 0) {
    stop(paste0("[Enrichment:readEnrichment] No snp list file found when at least one is needed."), call. = FALSE)
  }
  invisible()
}
