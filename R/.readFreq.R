.readFreq <- function(pattern, snpInfoDir) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  utils::read.delim(
    file = paste0(tmpDir, "/snpEnrichment/", pattern, ".all"), header = TRUE, stringsAsFactors = FALSE,
    colClasses = c("numeric", "character", "integer", "numeric"), na.string = c("NA", ""),
    check.names = FALSE, strip.white = TRUE, col.names = c("CHR", "SNP", "POS", "MAF")
  )
}
