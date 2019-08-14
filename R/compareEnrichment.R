methods::setGeneric(
  name = "compareEnrichment",
  def = function(object.x, object.y, pattern = "Chrom", nSample = 100, empiricPvalue = TRUE, mc.cores = 1, onlyGenome = TRUE) standardGeneric("compareEnrichment")
)

methods::setMethod(f = "compareEnrichment", signature = "ANY", definition = function(object.x, object.y, pattern = "Chrom", nSample = 100, empiricPvalue = TRUE, mc.cores = 1, onlyGenome = TRUE) {
  if (missing(object.x) | missing(object.y)) {
    stop('[Enrichment:compareEnrichment] "Enrichment" object is required.', call. = FALSE)
  }
  if (nSample < 10) {
    nSample <- 10
    warning('[Enrichment:compareEnrichment] "nSample" was increased to 10.', call. = FALSE)
  }
  if (is.enrichment(object.x) & is.enrichment(object.y)) {
    if (nrow(object.x["Data"]) == 0) {
      stop('[Enrichment:compareEnrichment] "Enrichment" data is empty for "object.x".', call. = FALSE)
    }
    if (nrow(object.y["Data"]) == 0) {
      stop('[Enrichment:compareEnrichment] "Enrichment" data is empty for "object.y".', call. = FALSE)
    }
  } else {
    stop('[Enrichment:compareEnrichment] "Enrichment" object is required.', call. = FALSE)
  }

  sigThresh.x <- object.x@Call$readEnrichment$sigThresh
  sigThresh.y <- object.y@Call$readEnrichment$sigThresh
  if (!identical(sigThresh.x, sigThresh.y)) {
    warning(paste0('[Enrichment:compareEnrichment] "sigThresh" differs from "object.x" to "object.y".\n         "object.x" parameter is: ', deparse(sigThresh.x)), call. = FALSE)
  }
  sigThresh <- sigThresh.x

  MAFpool.x <- object.x@Call$reSample$MAFpool
  MAFpool.y <- object.y@Call$reSample$MAFpool
  if (!identical(MAFpool.x, MAFpool.y)) {
    warning(paste0('[Enrichment:compareEnrichment] "MAFpool" differs from "object.x" to "object.y".\n         "object.x" parameter is: ', deparse(MAFpool.x)), call. = FALSE)
  }
  MAFpool <- eval(MAFpool.x)
  if (missing(MAFpool) | is.null(MAFpool)) MAFpool <- c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5)

  object.x@Call$reSample$empiricPvalue <- empiricPvalue
  object.y@Call$reSample$empiricPvalue <- empiricPvalue

  l1 <- object.x@eSNP@List
  l2 <- object.y@eSNP@List
  if (identical(l1, l2)) {
    stop("[Enrichment:compareEnrichment] Both lists are identical.", call. = FALSE)
  }

  if (is.null(object.x@Call$reSample$nSample) | is.null(object.y@Call$reSample$nSample)) {
    cat("########## Resample objects Start ##########\n")
    if (is.null(object.x@Call$reSample$nSample)) {
      cat("  Resampling object.x ...\n")
      .verbose(reSample(object = object.x, nSample = nSample, empiricPvalue = empiricPvalue, MAFpool = MAFpool, mc.cores = mc.cores, onlyGenome = onlyGenome))
    }
    if (is.null(object.y@Call$reSample$nSample)) {
      cat("  Resampling object.y ...\n")
      .verbose(reSample(object = object.y, nSample = nSample, empiricPvalue = empiricPvalue, MAFpool = MAFpool, mc.cores = mc.cores, onlyGenome = onlyGenome))
    }
    cat("########### Resample objects End ###########\n")
  }

  enrichObject1 <- object.x
  enrichObject2 <- object.y
  object.x <- reset(object.x, "Resampling")
  object.y <- reset(object.y, "Resampling")

  if (length(l1) < length(l2)) {
    object1 <- object.x
    object2 <- object.y
    namesRes <- c("Enrichment_1", "Enrichment_2", "PVALUE_1", "PVALUE_2", "PVALUE", "nSAMPLE", "SNP_1", "SNP_2")
  } else {
    object1 <- object.y
    object2 <- object.x
    namesRes <- c("Enrichment_2", "Enrichment_1", "PVALUE_2", "PVALUE_1", "PVALUE", "nSAMPLE", "SNP_2", "SNP_1")
  }
  rm(object.x, object.y)
  object2 <- reset(object2, "Data")
  object2 <- reset(object2, "LD")

  cat("############ Comparison Start ##############\n")
  if (onlyGenome == FALSE) {
    listRes <- eval(parse(text = paste0("list(", paste(paste0("Chrom", seq_len(22), " = NULL"), collapse = ", "), ")")))
    for (iChr in seq_len(22)) {
      cat("  Chromosome ", if (nchar(iChr) == 1) paste0("0", iChr) else iChr, ": ", sep = "")
      listRes[[iChr]] <- .compareEnrich(object1 = object1@Chromosomes[[iChr]], object2 = object2@Chromosomes[[iChr]], nSample = nSample, empiricPvalue = empiricPvalue, sigThresh = sigThresh, MAFpool = MAFpool, mc.cores = mc.cores)
      if (identical(sort(object1@Chromosomes[[iChr]]@eSNP@List), sort(object2@Chromosomes[[iChr]]@eSNP@List))) {
        listRes[[iChr]] <- reset(listRes[[iChr]], "Z")
        listRes[[iChr]] <- reset(listRes[[iChr]], "Resampling")
        listRes[[iChr]]@eSNP@PValue <- as.numeric(NA)
        listRes[[iChr]]@xSNP@PValue <- as.numeric(NA)
      }
      cat("END\n")
    }
  }

  cat("  Genome       : ")
  result <- .compareEnrich(object1 = object1, object2 = object2, nSample = nSample, empiricPvalue = empiricPvalue, sigThresh = sigThresh, MAFpool = MAFpool, mc.cores = mc.cores)
  if (onlyGenome == FALSE) result@Chromosomes <- listRes
  cat("END\n")

  res <- list(eSNP = NULL, xSNP = NULL)
  summaryObj1 <- print(object1, what = "All")
  summaryObj2 <- print(object2, what = "All")
  summaryRes <- print(result, what = "All")
  if (empiricPvalue) {
    res[["eSNP"]] <- cbind(summaryObj1[["eSNP"]][, "EnrichmentRatio"], summaryObj2[["eSNP"]][, "EnrichmentRatio"], summaryObj1[["eSNP"]][, "PValue"], summaryObj2[["eSNP"]][, "PValue"], summaryRes[["eSNP"]][, c("PValue", "nSample")], summaryObj1[["eSNP"]][, "eSNP"], summaryObj2[["eSNP"]][, "eSNP"])
    res[["xSNP"]] <- cbind(summaryObj1[["xSNP"]][, "EnrichmentRatio"], summaryObj2[["xSNP"]][, "EnrichmentRatio"], summaryObj1[["xSNP"]][, "PValue"], summaryObj2[["xSNP"]][, "PValue"], summaryRes[["xSNP"]][, c("PValue", "nSample")], summaryObj1[["xSNP"]][, "xSNP"], summaryObj2[["xSNP"]][, "xSNP"])
    colnames(res[["eSNP"]]) <- namesRes
    colnames(res[["xSNP"]]) <- namesRes
  } else {
    res[["eSNP"]] <- cbind(summaryObj1[["eSNP"]][, "EnrichmentRatio"], summaryObj2[["eSNP"]][, "EnrichmentRatio"], summaryObj1[["eSNP"]][, "PValue"], summaryObj2[["eSNP"]][, "PValue"], summaryRes[["eSNP"]][, c("PValue", "Z", "nSample")], summaryObj1[["eSNP"]][, "eSNP"], summaryObj2[["eSNP"]][, "eSNP"])
    res[["xSNP"]] <- cbind(summaryObj1[["xSNP"]][, "EnrichmentRatio"], summaryObj2[["xSNP"]][, "EnrichmentRatio"], summaryObj1[["xSNP"]][, "PValue"], summaryObj2[["xSNP"]][, "PValue"], summaryRes[["xSNP"]][, c("PValue", "Z", "nSample")], summaryObj1[["xSNP"]][, "xSNP"], summaryObj2[["xSNP"]][, "xSNP"])
    colnames(res[["eSNP"]]) <- c(namesRes[1:5], "Z", namesRes[6:8])
    colnames(res[["xSNP"]]) <- c(namesRes[1:5], "Z", namesRes[6:8])
  }

  cat("############# Comparison End ###############\n")
  warning("[Enrichment:compareEnrichment] This function is in development!", call. = FALSE)
  if (onlyGenome) {
    invisible(
      list(
        object.xy = print(result, what = "Genome"),
        object.x = print(enrichObject1, what = "Genome"),
        object.y = print(enrichObject2, what = "Genome")
      )
    )
  } else {
    invisible(
      list(
        object.xy = print(result),
        object.x = print(enrichObject1),
        object.y = print(enrichObject2)
      )
    )
  }
})
