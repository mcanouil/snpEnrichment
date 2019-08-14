initFiles <- function(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1) {
  if (missing(snpInfoDir) | missing(signalFile)) {
    stop("[Enrichment:initFiles] argument(s) missing.", call. = FALSE)
  }
  snpInfoDir <- .checkFilePath(snpInfoDir)
  FILES <- list.files(snpInfoDir, pattern = ".bim")
  .checkSnpInfoDir(snpInfoDir)
  .checkSignalFile(signalFile)
  tmpDir <- gsub("\\\\", "/", tempdir())
  dir.create(paste0(tmpDir, "/snpEnrichment/"), showWarnings = FALSE)
  cat("All files are ready for chromosome:\n  ")
  resParallel <- mclapply2(X = seq_len(22), mc.cores = min(22, mc.cores), FUN = function(iChr) {
    newPattern <- gsub(".bim", "", grep(paste0(pattern, iChr, "[^0-9]"), FILES, value = TRUE))
    err1 <- try(.writeSignal(pattern = newPattern, snpInfoDir = snpInfoDir, signalFile = signalFile), silent = TRUE)
    err2 <- try(.writeFreq(pattern = newPattern, snpInfoDir = snpInfoDir), silent = TRUE)
    cat(iChr, " ", sep = "")
    if (class(err1) == "try-error" | class(err2) == "try-error") invisible("ERROR") else invisible()
  })
  if (any(unlist(resParallel, use.names = FALSE) == "ERROR")) {
    stop("[Enrichment:initFiles] initialize files failed.", call. = FALSE)
  }
  cat("\n\n")
  invisible()
}