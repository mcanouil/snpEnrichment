.splitByChrom <- function(pattern, snpListFile, directory) {
  if (missing(snpListFile)) {
    stop("[Enrichment:readEnrichment] argument(s) missing.", call. = FALSE)
  }
  snpList <- utils::read.delim(
    file = snpListFile, header = FALSE, stringsAsFactors = FALSE,
    na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE
  )
  switch(EXPR = as.character(ncol(snpList)),
    "1" = stop('[Enrichment:readEnrichment] at least two columns are needed: "Chromosome" and "rs name".', call. = FALSE),
    "2" = colnames(snpList) <- c("CHR", "SNP"),
    "3" = colnames(snpList) <- c("CHR", "SNP", "TRANSCRIPT"),
    colnames(snpList)[1:3] <- c("CHR", "SNP", "TRANSCRIPT")
  )
  if (ncol(snpList) > 3) {
    snpList <- snpList[, c("CHR", "SNP", "TRANSCRIPT")]
  }
  filePathDetails <- unlist(strsplit(snpListFile, "/"), use.names = FALSE)
  if (is.null(directory)) {
    filePath <- paste0(c(filePathDetails[-length(filePathDetails)], ""), collapse = "/")
  } else {
    dir.create(paste0(directory, "snpList/"), showWarnings = FALSE)
    filePath <- paste0(directory, "snpList/")
  }
  by(snpList, snpList[, "CHR"], function(snpListChr) {
    utils::write.table(snpListChr, file = paste0(filePath, pattern, unique(snpListChr[, "CHR"]), "-", filePathDetails[length(filePathDetails)]), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  })
  filePath
}
