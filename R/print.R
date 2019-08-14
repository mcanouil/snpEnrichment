methods::setMethod(f = "print", signature = "EnrichSNP", definition = function(x) {
  EnrichmentRatio <- x@EnrichmentRatio
  Z <- x@Z
  PValue <- x@PValue
  Resampling <- nrow(x@Resampling)
  Data <- sum(x@Table)
  List <- length(x@List)
  resTmp <- c(
    if (length(EnrichmentRatio) == 0) NA else EnrichmentRatio,
    if (length(Z) == 0) NA else Z,
    if (length(PValue) == 0) NA else PValue,
    Resampling,
    Data,
    List
  )
  names(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", "eSNP")
  res <- t(resTmp)
  rownames(res) <- "eSNP"
  res
})

methods::setMethod(f = "print", signature = "Chromosome", definition = function(x, type = c("eSNP", "xSNP")) {
  if (missing(x)) {
    stop('[Chromosome:print] "x" is missing.', call. = FALSE)
  }
  if (is.null(type) | any(!type %in% c("eSNP", "xSNP"))) {
    stop('[Chromosome:print] "type" must be: "eSNP" and/or "xSNP".', call. = FALSE)
  }
  res <- list()
  for (iType in type) {
    resTmp <- print(x[iType])
    colnames(resTmp) <- c("EnrichmentRatio", "Z", "PValue", "nSample", "TotalSNP", iType)
    rownames(resTmp) <- paste("Chrom", iType, sep = ":")
    res[[iType]] <- resTmp
  }
  if (length(type) == 1) {
    res <- res[[1]]
  } else {
    res <- do.call("rbind", res)
    rownames(res) <- paste("Chrom", type, sep = ":")
  }
  res
})

methods::setMethod(
  f = "print",
  signature = "Enrichment",
  definition = function(x, what = "Genome", type = c("eSNP", "xSNP")) {
    if (missing(x)) {
      stop('[Enrichment:print] "x" is missing.', call. = FALSE)
    }
    if (is.null(what) | any(!what %in% c("All", "Genome", seq_len(22)))) {
      stop('[Enrichment:print] "what" must be: "All", "Genome" or numeric value (atomic or vector).', call. = FALSE)
    }
    if (is.null(type) | any(!type %in% c("eSNP", "xSNP"))) {
      stop('[Enrichment:print] "type" must be: "eSNP" and/or "xSNP".', call. = FALSE)
    }
    empiricPvalue <- x["Call"][["reSample"]][["empiricPvalue"]]
    if (is.null(empiricPvalue)) empiricPvalue <- FALSE

    if (length(what) == 1) {
      switch(EXPR = as.character(what),
        "Genome" = {
          res <- list()
          for (iType in type) {
            resTmp <- print(x[iType])
            colnames(resTmp) <- c("EnrichmentRatio", "Z", "PValue", "nSample", "TotalSNP", iType)
            rownames(resTmp) <- paste("Chrom", iType, sep = ":")
            if (empiricPvalue) {
              resTmp <- resTmp[, -grep("Z", colnames(resTmp))]
            }
            res[[iType]] <- resTmp
          }
          if (length(type) == 1) {
            res <- res[[1]]
          } else {
            res <- do.call("rbind", res)
            rownames(res) <- paste(what, type, sep = ":")
          }
          res
        },
        "All" = {
          res <- list()
          for (iType in type) {
            resTmp <- print(x[iType])
            tmp <- do.call("rbind", lapply(seq_len(22), function(n) {
              print(x["Chromosomes", n], type = iType)
            }))
            resTmp <- rbind(resTmp, tmp)
            rownames(resTmp) <- c("Genome", paste0("Chrom", seq_len(22)))
            colnames(resTmp) <- c("EnrichmentRatio", "Z", "PValue", "nSample", "TotalSNP", iType)
            if (empiricPvalue) {
              resTmp <- resTmp[, -grep("Z", colnames(resTmp))]
            }
            res[[iType]] <- resTmp
          }
          if (length(type) == 1) res <- res[[1]]
          res
        },
        {
          res <- print(x["Chromosomes", what], type = type)
          rownames(res) <- paste0("Chrom", what, ":", type)
          res
        }
      )
    } else {
      if (!is.numeric(what)) {
        stop('[Enrichment:print] "what" must be: "Genome", "All" or a numeric vector.', call. = FALSE)
      }
      resTmp <- lapply(what, function(iWhat) print(x["Chromosomes", iWhat], type = type))
      res <- do.call("rbind", resTmp)
      whatNames <- sapply(what, function(iWhat) paste0("Chrom", iWhat))
      rownames(res) <- paste0(rep(whatNames, each = length(type)), ":", type)
      if (empiricPvalue) res <- res[, -grep("Z", colnames(res))]
      res
    }
  }
)
