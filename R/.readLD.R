.readLD <- function(pattern, snpInfoDir, ldDir) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  if (missing(ldDir) | is.null(ldDir)) {
    IN <- paste0(tmpDir, "/snpEnrichment/", pattern, ".ld")
  } else {
    IN <- paste0(.checkFilePath(ldDir), pattern, ".ld")
  }
  nbCol <- as.character(ncol(utils::read.table(text = readLines(con = IN, n = 1), stringsAsFactors = FALSE)))
  switch(EXPR = nbCol,
    "2" = {
      utils::read.delim(
        file = IN, header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character"),
        na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, sep = "\t"
      )
    },
    "7" = {
      utils::read.delim(
        file = IN, header = TRUE, stringsAsFactors = FALSE, colClasses = c("NULL", "NULL", "character", "NULL", "NULL", "character", "NULL"),
        na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, sep = ""
      )
    },
    stop(paste0('[Enrichment:readEnrichment] "', pattern, '.ld" structure must be a matrix file with columns:\n       c("SNP_A", "SNP_B") or c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")".'), call. = FALSE)
  )
}
