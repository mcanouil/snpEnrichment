.writeFreq <- function(pattern, snpInfoDir) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  IN <- paste0(snpInfoDir, pattern)
  OUT <- paste0(tmpDir, "/snpEnrichment/", pattern)
  plinkData <- snpStats::read.plink(bed = paste0(IN, ".bed"), bim = paste0(IN, ".bim"), fam = paste0(IN, ".fam"), select.snps = .readSignal(pattern)[, "SNP"])
  plinkFreq <- snpStats::col.summary(plinkData$genotypes)
  plinkFreq <- cbind(snp.name = rownames(plinkFreq), MAF = plinkFreq[, "MAF"])
  plinkRes <- merge(plinkData$map, plinkFreq, by = "snp.name")
  plinkRes <- plinkRes[, c("chromosome", "snp.name", "position", "MAF")]
  plinkRes[, "MAF"] <- as.numeric(as.character(plinkRes[, "MAF"]))
  colnames(plinkRes) <- c("CHR", "SNP", "POS", "MAF")
  utils::write.table(plinkRes, paste0(OUT, ".all"), row.names = FALSE, sep = "\t")
  invisible()
}
