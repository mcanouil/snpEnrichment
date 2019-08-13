#' .writeSignal
#'
#' @param signalFile
#' @keywords internal
.writeSignal <- function(pattern, snpInfoDir, signalFile) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  if (length(unlist(strsplit(readLines(signalFile, n = 1), split = "\t"), use.names = FALSE)) > 1) {
    signal <- utils::read.delim(
      file = signalFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
      colClasses = c("character", "numeric"), na.string = c("NA", ""),
      check.names = FALSE, strip.white = TRUE, col.names = c("SNP", "PVALUE")
    )
  } else {
    if (length(unlist(strsplit(readLines(signalFile, n = 1), split = " "), use.names = FALSE)) > 1) {
      signal <- utils::read.delim(
        file = signalFile, header = TRUE, sep = " ", stringsAsFactors = FALSE,
        colClasses = c("character", "numeric"), na.string = c("NA", ""),
        check.names = FALSE, strip.white = TRUE, col.names = c("SNP", "PVALUE")
      )
    } else {
      stop('[Enrichment:initFiles] only " " and "\t" are allowed as columns separator in Signal file.', call. = FALSE)
    }
  }

  chrom.bim <- utils::read.delim(
    file = paste0(snpInfoDir, pattern, ".bim"), header = FALSE, stringsAsFactors = FALSE,
    colClasses = c("numeric", "character", "NULL", "NULL", "NULL", "NULL"), na.string = c("NA", ""),
    check.names = FALSE, strip.white = TRUE, col.names = c("CHR", "SNP", "", "", "", "")
  )
  signal <- merge(signal, chrom.bim, by = "SNP")[, c(3, 1, 2)]
  eval(parse(text = paste0(
    'write.table(signal, file = "', tmpDir, "/snpEnrichment/", pattern,
    '.signal", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")'
  )))
  invisible()
}
