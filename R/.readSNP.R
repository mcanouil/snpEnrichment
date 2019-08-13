#' .readSNP
#'
#' @param snpListDir
#' @keywords internal
.readSNP <- function(pattern, snpListDir) {
  snpListFile <- list.files(gsub("/$", "", snpListDir), pattern = paste0(pattern, "[^0-9]"), full.names = TRUE)
  if (is.na(snpListFile)) {
    snpList <- data.frame()
  } else {
    snpList <- utils::read.delim(
      file = snpListFile, header = FALSE, stringsAsFactors = FALSE,
      na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE
    )
    switch(EXPR = as.character(ncol(snpList)),
      "1" = colnames(snpList) <- "SNP",
      "2" = colnames(snpList) <- c("CHR", "SNP"),
      "3" = colnames(snpList) <- c("CHR", "SNP", "TRANSCRIPT"),
      colnames(snpList)[1:3] <- c("CHR", "SNP", "TRANSCRIPT")
    )
    if (ncol(snpList) >= 3) snpList <- snpList[, c("SNP", "TRANSCRIPT")]
    if (nrow(snpList) != 0) snpList[, "eSNP"] <- 1
  }
  snpList
}
