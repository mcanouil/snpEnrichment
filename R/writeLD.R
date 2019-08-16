#' Linkage Disequilibrium (LD) computation with PLINK
#'
#' [writeLD] write a '.ld' file for each chromosomes which contains the LD (r^2).
#'
#' @inheritParams readEnrichment
#' @param ldThresh Threshold value for LD calculation.
#' @param depth This parameter is mandatory and controls the maximum lag between SNPs considered.
#'
#' @return One ".ld" file per chromosome is returned by [writeLD] in `snpInfoDir` directory.
#'
#' @note The LD computation can take a long time depending on number of SNPs in
#'   `signalFile`. It is recommended to save LD results in a directory
#'   (`ldDir`) which is not a temporary directory.
#'
#' @examples
#' if (interactive()) {
#'   signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment")
#'   snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment")
#'   writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1)
#' }
#'
#' @export
writeLD <- function(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, depth = 1000, mc.cores = 1) {
  if (missing(pattern) | missing(snpInfoDir) | missing(signalFile) | missing(ldThresh)) {
    stop("[Enrichment:writeLD] argument(s) missing.", call. = FALSE)
  }
  tmpDir <- gsub("\\\\", "/", tempdir())
  dir.create(paste0(tmpDir, "/snpEnrichment/"), showWarnings = FALSE)
  snpInfoDir <- .checkFilePath(snpInfoDir)
  .checkSnpInfoDir(snpInfoDir)
  .checkSignalFile(signalFile)
  if (missing(ldDir) | is.null(ldDir)) {
    ldDir <- paste0(tmpDir, "/snpEnrichment/")
  } else {
    ldDir <- .checkFilePath(ldDir)
  }
  FILES <- list.files(snpInfoDir, pattern = ".bim")
  cat("Compute LD for chromosome:\n  ")
  resParallel <- .mclapply(X = seq_len(22), mc.cores = min(22, mc.cores), FUN = function(iChr) {
    newPattern <- gsub(".bim", "", grep(paste0(pattern, iChr, "[^0-9]"), FILES, value = TRUE))
    isThereSignals <- list.files(paste0(tmpDir, "/snpEnrichment/"), full.names = TRUE, pattern = ".signal")
    if (length(isThereSignals) != 22) .writeSignal(pattern = newPattern, snpInfoDir, signalFile)
    IN <- paste0(snpInfoDir, newPattern)
    OUT <- paste0(tmpDir, "/snpEnrichment/", newPattern)
    plinkData <- snpStats::read.plink(bed = paste0(IN, ".bed"), bim = paste0(IN, ".bim"), fam = paste0(IN, ".fam"), select.snps = .readSignal(newPattern)[, "SNP"])
    ldData <- snpStats::ld(x = plinkData$genotypes, depth = min(ncol(plinkData$genotypes) - 1, depth), stats = "R.squared")
    if (any(isNA <- is.na(ldData))) ldData <- replace(ldData, grep(TRUE, isNA), 0)
    ldData <- apply(ldData, 1, function(li) which(li > ldThresh))
    resLD <- data.frame(matrix(unlist(strsplit(names(unlist(ldData)), "\\."), use.names = FALSE), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("SNP_A", "SNP_B"))), stringsAsFactors = FALSE)
    resLD <- resLD[which(resLD[, 1] != resLD[, 2]), ]
    utils::write.table(resLD[, c("SNP_A", "SNP_B")], file = paste0(ldDir, newPattern, ".ld"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat(iChr, " ", sep = "")
  })
  cat("\n\n")
  invisible()
}
