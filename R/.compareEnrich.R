.compareEnrich <- function(object1, object2, nSample, empiricPvalue, sigThresh, MAFpool, mc.cores) {
  DATA <- object1["Data"]
  chrLD <- object1["LD"]
  isLD <- length(chrLD) != 0
  eEnrichment <- object1@eSNP@Table
  eEnrichRatio <- object1@eSNP@EnrichmentRatio
  eSNPlist <- object1@eSNP@List
  eList <- union(object1@eSNP@List, object2@eSNP@List)
  if (isLD) {
    xList <- union(object1@xSNP@List, object2@xSNP@List)
    data <- DATA[DATA[, "SNP"] %in% union(eList, xList), ]
  } else {
    data <- DATA[DATA[, "SNP"] %in% eList, ]
  }

  data[, "MAFpool"] <- NA
  data[, "MAFpool"] <- as.factor(cut(data[, "MAF"], breaks = MAFpool, labels = FALSE, include.lowest = TRUE))
  nPool <- nlevels(data$MAFpool)
  eSNPlistPool <- table(data[eSNPlist, "MAFpool"])

  eSNPsample <- NULL
  nSampleMin <- min(nSample, max(1000, nSample * 10 / 100))
  if (isLD) {
    xEnrichment <- object1@xSNP@Table
    xEnrichRatio <- object1@xSNP@EnrichmentRatio
    xSNPlist <- object1@xSNP@List
    xSNPlistPool <- table(data[xSNPlist, "MAFpool"])
    xSNPsample <- NULL
  }

  popSNP4Sample <- split(eList, data[eList, "MAFpool"])
  assoc_eSNP <- factor(data[eList, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
  if (isLD) {
    assoc_xSNP <- factor(data[xList, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
  }
  rm(data)

  cat("0.. ")
  resParallel <- mclapply2(X = seq_len(nSampleMin), mc.cores = mc.cores, FUN = function(i) {
    eSNPlistRandom <- unlist(sapply(seq_len(nPool), function(g) {
      sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))
    }), use.names = FALSE)
    eSNPenrichStats <- table(assoc_eSNP, factor(eList %in% eSNPlistRandom, levels = c(FALSE, TRUE)))
    eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
    if (isLD) {
      xSNPlistRandom <- intersect(xList, unique(chrLD[which(names(chrLD) %in% eSNPlistRandom)]))
      xSNPenrichStats <- table(assoc_xSNP, factor(xList %in% xSNPlistRandom, levels = c(FALSE, TRUE)))
      xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
    } else {
      xTmp <- NULL
    }
    list(eSNP = eTmp, xSNP = xTmp)
  })
  eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function(l) l$eSNP)))

  eSNPanyDup <- eSNPsample[!duplicated(eSNPsample), ]
  if (is.matrix(eSNPanyDup)) {
    if (length(which(is.na(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.na(eSNPsample[, 5]))) > 0) {
      eSNPsample[is.na(eSNPsample[, 5]), 5] <- rep(mean(eSNPsample[(!is.na(eSNPsample[, 5]) & is.finite(eSNPsample[, 5])), 5]), length(eSNPsample[is.na(eSNPsample[, 5]), 5]))
    }
    if (length(which(is.infinite(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.infinite(eSNPsample[, 5]))) > 0) {
      eSNPsample[is.infinite(eSNPsample[, 5]), 5] <- rep(max(eSNPsample[!is.infinite(eSNPsample[, 5]), 5]) * 1.05, length(eSNPsample[is.infinite(eSNPsample[, 5]), 5]))
    }
    eMeanEnrichRatio <- sum(eSNPsample[, 5]) / length(eSNPsample[, 5])
    Ze <- (eEnrichRatio - eMeanEnrichRatio) / sqrt((sum((eSNPsample[, 5] - eMeanEnrichRatio)^2)) / (length(eSNPsample[, 5]) - 1))
  } else {
    Ze <- as.numeric(NA)
  }

  if (isLD) {
    xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l) l$xSNP)))
    xSNPanyDup <- xSNPsample[!duplicated(xSNPsample), ]
    if (is.matrix(xSNPanyDup)) {
      if (length(which(is.na(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.na(xSNPsample[, 5]))) > 0) {
        xSNPsample[is.na(xSNPsample[, 5]), 5] <- rep(mean(xSNPsample[(!is.na(xSNPsample[, 5]) & is.finite(xSNPsample[, 5])), 5]), length(xSNPsample[is.na(xSNPsample[, 5]), 5]))
      }
      if (length(which(is.infinite(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.infinite(xSNPsample[, 5]))) > 0) {
        xSNPsample[is.infinite(xSNPsample[, 5]), 5] <- rep(max(xSNPsample[!is.infinite(xSNPsample[, 5]), 5]) * 1.05, length(xSNPsample[is.infinite(xSNPsample[, 5]), 5]))
      }
      xMeanEnrichRatio <- sum(xSNPsample[, 5]) / length(xSNPsample[, 5])
      Zx <- (xEnrichRatio - xMeanEnrichRatio) / sqrt((sum((xSNPsample[, 5] - xMeanEnrichRatio)^2)) / (length(xSNPsample[, 5]) - 1))
    } else {
      Zx <- as.numeric(NA)
    }
  } else {
    Zx <- as.numeric(NA)
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
      eSNPenrichStats <- table(assoc_eSNP, factor(eList %in% eSNPlistRandom, levels = c(FALSE, TRUE)))
      eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
      if (isLD) {
        xSNPlistRandom <- intersect(xList, unique(chrLD[which(names(chrLD) %in% eSNPlistRandom)]))
        xSNPenrichStats <- table(assoc_xSNP, factor(xList %in% xSNPlistRandom, levels = c(FALSE, TRUE)))
        xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
      } else {
        xTmp <- NULL
      }
      return(list(eSNP = eTmp, xSNP = xTmp))
    })
    eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function(l) l$eSNP)))

    eSNPanyDup <- eSNPsample[!duplicated(eSNPsample), ]
    if (is.matrix(eSNPanyDup)) {
      if (length(which(is.na(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.na(eSNPsample[, 5]))) > 0) {
        eSNPsample[is.na(eSNPsample[, 5]), 5] <- rep(mean(eSNPsample[(!is.na(eSNPsample[, 5]) & is.finite(eSNPsample[, 5])), 5]), length(eSNPsample[is.na(eSNPsample[, 5]), 5]))
      }
      if (length(which(is.infinite(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.infinite(eSNPsample[, 5]))) > 0) {
        eSNPsample[is.infinite(eSNPsample[, 5]), 5] <- rep(max(eSNPsample[is.finite(eSNPsample[, 5]), 5]) * 1.05, length(eSNPsample[is.infinite(eSNPsample[, 5]), 5]))
      }
      eMeanEnrichRatio <- sum(eSNPsample[, 5]) / length(eSNPsample[, 5])
      Ze <- (eEnrichRatio - eMeanEnrichRatio) / sqrt((sum((eSNPsample[, 5] - eMeanEnrichRatio)^2)) / (length(eSNPsample[, 5]) - 1))
    } else {
      Ze <- as.numeric(NA)
    }

    if (isLD) {
      xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l) l$xSNP)))
      xSNPanyDup <- xSNPsample[!duplicated(xSNPsample), ]
      if (is.matrix(xSNPanyDup)) {
        if (length(which(is.na(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.na(xSNPsample[, 5]))) > 0) {
          xSNPsample[is.na(xSNPsample[, 5]), 5] <- rep(mean(xSNPsample[(!is.na(xSNPsample[, 5]) & is.finite(xSNPsample[, 5])), 5]), length(xSNPsample[is.na(xSNPsample[, 5]), 5]))
        }
        if (length(which(is.infinite(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.infinite(xSNPsample[, 5]))) > 0) {
          xSNPsample[is.infinite(xSNPsample[, 5]), 5] <- rep(max(xSNPsample[is.finite(xSNPsample[, 5]), 5]) * 1.05, length(xSNPsample[is.infinite(xSNPsample[, 5]), 5]))
        }
        xMeanEnrichRatio <- sum(xSNPsample[, 5]) / length(xSNPsample[, 5])
        Zx <- (xEnrichRatio - xMeanEnrichRatio) / sqrt((sum((xSNPsample[, 5] - xMeanEnrichRatio)^2)) / (length(xSNPsample[, 5]) - 1))
      } else {
        Zx <- as.numeric(NA)
      }
    } else {
      Zx <- as.numeric(NA)
    }
    iSample <- iSample + nSampleMin
    if ((iSample - (iSample %/% catStep * catStep)) == 0) {
      cat(iSample, ".. ", sep = "")
    }
    rm(resParallel)
  }

  whichPvalue <- if (empiricPvalue) 2 else 1

  empiricPvalue.eSNP <- sum(eEnrichRatio < eSNPsample[, 5]) / length(eSNPsample[, 5])
  statisticPvalue.eSNP <- pnorm(Ze, lower.tail = FALSE)
  object1@eSNP <- enrichSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = c(Distribution = statisticPvalue.eSNP, Empirical = empiricPvalue.eSNP)[whichPvalue], Resampling = eSNPsample)
  if (isLD) {
    empiricPvalue.xSNP <- sum(xEnrichRatio < xSNPsample[, 5]) / length(xSNPsample[, 5])
    statisticPvalue.xSNP <- pnorm(Zx, lower.tail = FALSE)
    object1@xSNP <- enrichSNP(List = xSNPlist, Table = xEnrichment, EnrichmentRatio = xEnrichRatio, Z = Zx, PValue = c(Distribution = statisticPvalue.xSNP, Empirical = empiricPvalue.xSNP)[whichPvalue], Resampling = xSNPsample)
  } else {
    object1@xSNP <- enrichSNP()
  }
  if (is.enrichment(object1)) {
    object1@Chromosomes <- lapply(object1@Chromosomes, reset, "Z")
    object1@Chromosomes <- lapply(object1@Chromosomes, reset, "PValue")
    object1@Chromosomes <- lapply(object1@Chromosomes, reset, "Resampling")
  }
  object1
}
