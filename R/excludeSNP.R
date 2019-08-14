methods::setGeneric(
  name = "excludeSNP",
  def = function(object, excludeFile, mc.cores = 1) standardGeneric("excludeSNP")
)

methods::setMethod(f = "excludeSNP", signature = "ANY", definition = function(object, excludeFile, mc.cores = 1) {
  if (!is.enrichment(object)) {
    stop('[Method:excludeSNP] not available for "', class(object), '" object.', call. = FALSE)
  }
})

methods::setMethod(f = "excludeSNP", signature = "Enrichment", definition = function(object, excludeFile, mc.cores = 1) {
  if (missing(excludeFile)) {
    stop('[Enrichment:excludeSNP] argument "excludeFile" is missing.', call. = FALSE)
  }
  cat("########## Exclude SNP list Start ##########\n")
  if (all(class(try(close(file(excludeFile)), silent = TRUE)) != "try-error")) {
    eSNPexclude <- utils::read.delim(file = excludeFile, header = FALSE, na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, stringsAsFactors = FALSE, sep = "\t")
  } else {
    eSNPexclude <- excludeFile
  }
  if (class(eSNPexclude) %in% c("matrix", "data.frame")) {
    if (ncol(eSNPexclude) > 1) {
      eSNPexclude <- eSNPexclude[, 1]
    }
  }
  eSNPexclude <- unlist(eSNPexclude, use.names = FALSE)
  callLD <- object@Call$readEnrichment$LD
  resParallel <- mclapply2(X = seq_len(22), mc.cores = min(22, mc.cores), FUN = function(iChr) {
    chrObject <- eval(parse(text = paste0("object@Chromosomes$Chrom", iChr)))
    temp <- chrObject@Data
    if (callLD) {
      xSNPexclude <- intersect(temp[, "SNP"], unique(c(eSNPexclude, chrObject@LD[names(chrObject@LD) %in% eSNPexclude])))
    } else {
      xSNPexclude <- eSNPexclude
    }
    temp[temp[, "SNP"] %in% xSNPexclude, "eSNP"] <- 0
    temp[temp[, "SNP"] %in% xSNPexclude, "xSNP"] <- 0
    chromosome(Data = temp, LD = chrObject@LD)
  })
  names(resParallel) <- paste0("Chrom", seq_len(22))
  result <- enrichment(
    Loss = object@Loss,
    Chromosomes = resParallel,
    Call = list(
      readEnrichment = object@Call$readEnrichment,
      reSample = list(
        object = NULL, nSample = NULL, empiricPvalue = NULL,
        MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL
      )
    )
  )
  rm(resParallel)
  GC()

  result <- computeER(object = result, sigThresh = object@Call$readEnrichment$sigThresh, mc.cores = mc.cores)
  cat("########### Update SNP list END ############\n")
  for (iType in c("eSNP", "xSNP")) {
    cat("   ", length(setdiff(object[iType]@List, result[iType]@List)), " SNPs are removed from", iType, "list.\n")
  }
  result@Loss <- cbind(result@Loss,
    exclude = c(
      result@Loss["Signal", "CIS"],
      length(result["List", seq_len(22)][["eSNP"]]),
      sapply(seq_len(22), function(iChr) {
        length(result["List", iChr][["eSNP"]])
      })
    )
  )

  cat("########### Exclude SNP list END ###########\n")
  result
})
