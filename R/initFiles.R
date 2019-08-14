#' Initialize files for enrichment analysis
#'
#' [initFiles] create several files needed to run [readEnrichment].
#' ".frq" and ".signal" are created with PLINK.
#' LD computation can be run with [writeLD] or with PLINK.
#'
#'
#' @param pattern A character string containing a expression to be matched with all
#'   chromosomes files (*e.g.*,"Chrom" for files which start by "Chrom" followed
#'   by the chromosome number).
#' @param snpInfoDir A character string naming a directory containing the reference data
#'   in a PLINK format (.bed, .bim and .fam).
#' @param signalFile The name of the signal file which the data are to be read from
#'   (2 columns: "SNP" and "PVALUE"). Each row of the table appears as one line of the file.
#'   If it does not contain an `_absolute_` path, the file name is `_relative_` to the current
#'   working directory, [getwd]. The fields separator character have to be a space `" "`
#'   or a tabulation `"\t"`.
#' @param mc.cores The number of cores to use (default is `1`), *i.e.*,
#'   at most how many child processes will be run simultaneously.
#'   Must be at least one, and parallelization requires at least two cores.
#'
#' @return This function writes several files, in the temporary directory (defined in `R_SESSION_TMPDIR`),
#'   nothing else is returned. These files are used to build an [Enrichment-class] object by [readEnrichment]
#'   in order to compute enrichment analysis ([reSample]).
#'
#' @examples
#' if (interactive()) {
#'   snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment")
#'   signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment")
#'   initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)
#' }
#'
#' @export
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
