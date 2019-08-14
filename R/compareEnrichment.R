#' Compare enrichment analysis between two SNPs list
#'
#' Compare the enrichment analysis between two set of SNPs.
#' [compareEnrichment] compare two [Enrichment] objects.
#'
#' @param object.x,object.y An [Enrichment] object fully filled (*e.g.*, [readEnrichment]).
#' @param pattern A character string containing a expression to be matched
#'   with all chromosomes files (*e.g.*, "Chrom" for files which start by "Chrom"
#'   followed by the chromosome number).
#' @param nSample The number of resampling done by [reSample] for p-values computation (minimum is 100).
#' @param empiricPvalue `empiricPvalue = TRUE` (default) compute PValue based on the null
#'   distribution (resampling). If `empiricPvalue = TRUE`, the empirical p-values are computed instead.
#' @param mc.cores The number of cores to use (default is `1`),
#'   *i.e.*, at most how many child processes will be run simultaneously.
#'   Must be at least one, and parallelization requires at least two cores.
#' @param onlyGenome `onlyGenome = TRUE` (default) compute resampling step for all chromosomes.
#'
#' @return Return a `list` of three elements:
#' + object.xy [Enrichment] object from the comparison between `object.x` and `object.y`.
#' + object.x [Enrichment] object passed in `object.x` with resampling data.
#' + object.y [Enrichment] object passed in `object.y` with resampling data.
#'
#' @note Still in development.
#'
#' @examples
#' if (interactive()) {
#'   data(toyEnrichment)
#'
#'   reSample(
#'     object = toyEnrichment,
#'     nSample = 10,
#'     empiricPvalue = TRUE,
#'     MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
#'     mc.cores = 1,
#'     onlyGenome = TRUE
#'   )
#'
#'   excludeFile <- c(
#'     "rs7897180", "rs4725479", "rs315404", "rs17390391", "rs1650670",
#'     "rs6783390", "rs1642009", "rs4756586", "rs11995037", "rs4977345",
#'     "rs13136448", "rs4233536", "rs11151079", "rs2299657", "rs4833930",
#'     "rs1384", "rs7168184", "rs6909895", "rs7972667", "rs2293229",
#'     "rs918216", "rs6040608", "rs2817715", "rs13233541", "rs4486743",
#'     "rs2127806", "rs10912854", "rs1869052", "rs9853549", "rs448658",
#'     "rs2451583", "rs17483288", "rs10962314", "rs9612059", "rs1384182",
#'     "rs8049208", "rs12215176", "rs2980996", "rs1736976", "rs8089268",
#'     "rs10832329", "rs12446540", "rs7676237", "rs869922", "rs16823426",
#'     "rs1374393", "rs13268781", "rs11134505", "rs7325241", "rs7520109"
#'   )
#'   # OR
#'   excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment")
#'
#'   toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)
#'
#'   compareResults <- compareEnrichment(
#'     object.x = toyEnrichment,
#'     object.y = toyEnrichment_exclude,
#'     pattern = "Chrom",
#'     nSample = 10,
#'     empiricPvalue = FALSE,
#'     mc.cores = 1,
#'     onlyGenome = TRUE
#'   )
#' }
#'
#' @export
compareEnrichment <- function(
  object.x,
  object.y,
  pattern = "Chrom",
  nSample = 100,
  empiricPvalue = TRUE,
  mc.cores = 1,
  onlyGenome = TRUE
) {
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
}
