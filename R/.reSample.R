#' .reSample
#'
#' @param mc.cores
#' @keywords internal
.reSample <- function(object, nSample, empiricPvalue, sigThresh, MAFpool, mc.cores) {
  nResampling <- nrow(object@eSNP@Resampling)
  data <- object["Data"]
  chrLD <- object["LD"]
  isLD <- length(chrLD) != 0
  eEnrichment <- object@eSNP@Table
  eEnrichRatio <- object@eSNP@EnrichmentRatio
  eSNPlist <- object@eSNP@List
  data[, "MAFpool"] <- NA
  data[, "MAFpool"] <- as.factor(cut(data[, "MAF"], breaks = MAFpool, labels = FALSE, include.lowest = TRUE))
  nPool <- nlevels(data[, "MAFpool"])
  eSNPlistPool <- table(data[data[, "eSNP"] == 1, "MAFpool"])

  if (nResampling == 0) {
    eSNPsample <- NULL
    nSampleMin <- min(nSample, max(1000, nSample * 10 / 100))
  } else {
    eSNPsample <- object@eSNP@Resampling
    nSampleMin <- nSample
  }

  if (isLD) {
    xSNPsample <- if (nResampling == 0) NULL else object@xSNP@Resampling
    xEnrichment <- object@xSNP@Table
    xEnrichRatio <- object@xSNP@EnrichmentRatio
    xSNPlist <- object@xSNP@List
    xSNPlistPool <- table(data[data[, "xSNP"] == 1, "MAFpool"], deparse.level = 0)
    popSNP4Sample <- split(data[data[, "xSNP"] != 1, "SNP"], data[data[, "xSNP"] != 1, "MAFpool"])
  } else {
    popSNP4Sample <- split(data[data[, "eSNP"] != 1, "SNP"], data[data[, "eSNP"] != 1, "MAFpool"])
  }
  assoc <- factor(data[, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
  SNPlist <- data[, "SNP"]
  rm(data)

  cat(0, ".. ", sep = "")
  resParallel <- mclapply2(X = seq_len(nSampleMin), mc.cores = mc.cores, FUN = function(i) {
    eSNPlistRandom <- unlist(sapply(seq_len(nPool), function(g) {
      sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))
    }), use.names = FALSE)
    eSNPenrichStats <- table(assoc, factor(SNPlist %in% eSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
    if (isLD) {
      xSNPlistRandom <- intersect(SNPlist, unique(chrLD[which(names(chrLD) %in% eSNPlistRandom)]))
      xSNPenrichStats <- table(assoc, factor(SNPlist %in% xSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
      xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
    } else {
      xTmp <- NULL
    }
    eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
    return(list(eSNP = eTmp, xSNP = xTmp))
  })
  eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function(l) {
    l$eSNP
  })))
  eMeanEnrichRatio <- sum(eSNPsample[, 5]) / length(eSNPsample[, 5])
  Ze <- (eEnrichRatio - eMeanEnrichRatio) / sqrt((sum((eSNPsample[, 5] - eMeanEnrichRatio)^2)) / (length(eSNPsample[, 5]) - 1))
  if (isLD) {
    xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l) {
      l$xSNP
    })))
    xMeanEnrichRatio <- sum(xSNPsample[, 5]) / length(xSNPsample[, 5])
    Zx <- (xEnrichRatio - xMeanEnrichRatio) / sqrt((sum((xSNPsample[, 5] - xMeanEnrichRatio)^2)) / (length(xSNPsample[, 5]) - 1))
  } else {
    Zx <- 1
  }
  iSample <- nSampleMin
  cat(iSample, ".. ", sep = "")
  catStep <- (nSample / 10)
  rm(resParallel)

  while (iSample < nSample) {
    resParallel <- mclapply2(X = seq_len(nSampleMin), mc.cores = mc.cores, FUN = function(i) {
      eSNPlistRandom <- unlist(sapply(seq_len(nPool), function(g) {
        sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))
      }), use.names = FALSE)
      eSNPenrichStats <- table(assoc, factor(SNPlist %in% eSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
      if (isLD) {
        xSNPlistRandom <- intersect(SNPlist, unique(chrLD[which(names(chrLD) %in% eSNPlistRandom)]))
        xSNPenrichStats <- table(assoc, factor(SNPlist %in% xSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
        xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
      } else {
        xTmp <- NULL
      }
      eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
      list(eSNP = eTmp, xSNP = xTmp)
    })
    eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function(l) l$eSNP)))
    eMeanEnrichRatio <- sum(eSNPsample[, 5]) / length(eSNPsample[, 5])
    Ze <- (eEnrichRatio - eMeanEnrichRatio) / sqrt((sum((eSNPsample[, 5] - eMeanEnrichRatio)^2)) / (length(eSNPsample[, 5]) - 1))
    if (isLD) {
      xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l) l$xSNP)))
      xMeanEnrichRatio <- sum(xSNPsample[, 5]) / length(xSNPsample[, 5])
      Zx <- (xEnrichRatio - xMeanEnrichRatio) / sqrt((sum((xSNPsample[, 5] - xMeanEnrichRatio)^2)) / (length(xSNPsample[, 5]) - 1))
    } else {
      Zx <- 1
    }
    iSample <- iSample + nSampleMin
    if ((iSample - (iSample %/% catStep * catStep)) == 0) {
      cat(iSample, ".. ", sep = "")
    }
    rm(resParallel)
  }

  whichPvalue <- if (empiricPvalue) 2 else whichPvalue <- 1

  empiricPvalue.eSNP <- sum(eEnrichRatio < eSNPsample[, 5]) / length(eSNPsample[, 5])
  statisticPvalue.eSNP <- pnorm(Ze, lower.tail = FALSE)
  object@eSNP <- enrichSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = c(Distribution = statisticPvalue.eSNP, Empirical = empiricPvalue.eSNP)[whichPvalue], Resampling = eSNPsample)
  if (isLD) {
    empiricPvalue.xSNP <- sum(xEnrichRatio < xSNPsample[, 5]) / length(xSNPsample[, 5])
    statisticPvalue.xSNP <- pnorm(Zx, lower.tail = FALSE)
    object@xSNP <- enrichSNP(List = xSNPlist, Table = xEnrichment, EnrichmentRatio = xEnrichRatio, Z = Zx, PValue = c(Distribution = statisticPvalue.xSNP, Empirical = empiricPvalue.xSNP)[whichPvalue], Resampling = xSNPsample)
  } else {
    object@xSNP <- enrichSNP()
  }
  object
}
