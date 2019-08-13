methods::setGeneric(
  name = "computeER",
  def = function(object, sigThresh = 0.05, mc.cores = 1) standardGeneric("computeER")
)

methods::setMethod(f = "computeER", signature = "Chromosome", definition = function(object, sigThresh = 0.05, mc.cores = 1) {
  if (missing(object)) {
    stop('[Chromosome:computeER] "Chromosome" object is required.', call. = FALSE)
  }
  object@Chromosomes <- mclapply2(object@Chromosomes, mc.cores = mc.cores, function(chr) {
    data <- chr@Data
    chrLD <- length(chr@LD)
    for (type in c("eSNP", "xSNP")) {
      if (!(chrLD == 0 & type == "xSNP")) {
        snpEnrich <- table(factor(chr@Data[, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE)), factor(chr@Data[, type], levels = c(0, 1)))
        colnames(snpEnrich) <- c("otherSNP", type)
        rownames(snpEnrich) <- eval(parse(text = paste0('c("P>=', sigThresh, '", "P<', sigThresh, '")')))
        chr[type] <- enrichSNP(List = chr[type]@List, Table = unclass(snpEnrich), EnrichmentRatio = .enrichmentRatio(snpEnrich))
      }
    }
    chr
  })
  object
})

methods::setMethod(f = "computeER", signature = "Enrichment", definition = function(object, sigThresh = 0.05, mc.cores = 1) {
  if (missing(object)) {
    stop('[Enrichment:computeER] "Enrichment" object is required.', call. = FALSE)
  }
  rowNames <- c(paste0("P>=", sigThresh), paste0("P<", sigThresh))
  object@Chromosomes <- mclapply2(object@Chromosomes, mc.cores = mc.cores, function(chr) {
    data <- chr@Data
    chrLD <- length(chr@LD)
    pvalFactor <- factor(data[, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
    for (iType in c("eSNP", "xSNP")) {
      if (!(chrLD == 0 & iType == "xSNP")) {
        snpEnrich <- table(pvalFactor, factor(data[, iType], levels = c(0, 1)))
        colnames(snpEnrich) <- c("otherSNP", iType)
        rownames(snpEnrich) <- rowNames
        chr[iType] <- enrichSNP(List = chr[iType]@List, Table = unclass(snpEnrich), EnrichmentRatio = .enrichmentRatio(snpEnrich))
      }
    }
    chr
  })
  for (iType in c("eSNP", "xSNP")) {
    bigEnrichment <- matrix(rowSums(sapply(seq_len(22), function(jChr) {
      object@Chromosomes[[jChr]][iType]@Table
    })), nrow = 2, ncol = 2)
    object[iType]@Table <- bigEnrichment
    object[iType]@EnrichmentRatio <- .enrichmentRatio(bigEnrichment)
  }
  object
})
