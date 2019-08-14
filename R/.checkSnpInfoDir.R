.checkSnpInfoDir <- function(snpInfoDir) {
  snpInfoDir <- .checkFilePath(snpInfoDir)
  if (length(list.files(snpInfoDir, pattern = "*.bim")) != 22) {
    stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.bim")), " 'bim' files found when 22 is needed."), call. = FALSE)
  }
  if (length(list.files(snpInfoDir, pattern = "*.bed")) != 22) {
    stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.bed")), " 'bed' files found when 22 is needed."), call. = FALSE)
  }
  if (length(list.files(snpInfoDir, pattern = "*.fam")) != 22) {
    stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.fam")), " 'fam' files found when 22 is needed."), call. = FALSE)
  }
  invisible()
}
