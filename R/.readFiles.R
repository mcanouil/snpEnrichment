#' .readFiles
#'
#' @param distThresh
#' @keywords internal
.readFiles <- function(pattern, snpInfoDir, snpListDir, distThresh) {
  fullPattern <- gsub(".bim", "", grep(paste0(pattern, "[^0-9]"), list.files(snpInfoDir, pattern = ".bim"), value = TRUE))
  signalPattern <- gsub(".signal", "", grep(paste0(pattern, "[^0-9]"), list.files(paste0(gsub("\\\\", "/", tempdir()), "/snpEnrichment/"), pattern = ".signal"), value = TRUE))
  eSNP <- .readSNP(pattern = fullPattern, snpListDir = snpListDir)
  signal <- .readSignal(pattern = signalPattern)
  plinkData <- .readFreq(pattern = fullPattern, snpInfoDir = snpInfoDir)

  signalPlink <- merge(signal, plinkData, by = "SNP")
  if (nrow(eSNP) != 0) {
    eSNPunique <- eSNP[!duplicated(eSNP[, "SNP"]), ]
    snpSignal <- merge(signalPlink, eSNPunique[, c("SNP", "eSNP")], by = "SNP", all.x = TRUE)
    snpSignal[, "eSNP"][is.na(snpSignal[, "eSNP"])] <- 0
    snpLoss <- c(
      length(eSNP[, "SNP"]),
      length(unique(eSNP[, "SNP"])),
      length(unique(intersect(eSNP[, "SNP"], signal[!is.na(signal[, "PVALUE"]), "SNP"])))
    )
  } else {
    snpSignal <- signalPlink
    snpSignal[, "eSNP"] <- 0
    snpLoss <- c(0, 0, 0)
  }

  temp <- unique(snpSignal[, c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP")])
  data <- transform(temp, SNP = as.character(temp$SNP), PVALUE = as.numeric(temp$PVALUE), CHR = as.numeric(temp$CHR), POS = as.integer(temp$POS), MAF = as.numeric(temp$MAF), eSNP = as.numeric(temp$eSNP))
  data[, "xSNP"] <- 0
  if (any(duplicated(data[, "SNP"]))) {
    cat(data[, "SNP"][duplicated(data[, "SNP"])], "\n")
    stop("[Enrichment:readEnrichment] Duplicated SNPs in Signal.", call. = FALSE)
  } else {
    rownames(data) <- data[, "SNP"]
  }
  list(data = data, snpLoss = snpLoss)
}
