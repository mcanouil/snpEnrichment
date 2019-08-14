.checkSignalFile <- function(signalFile) {
  if (!file.exists(signalFile)) {
    stop(paste0("[Enrichment:initFiles] ", signalFile, " doesn't exist."), call. = FALSE)
  }
  invisible()
}
