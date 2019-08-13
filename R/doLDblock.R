methods::setGeneric(
  name = "doLDblock",
  def = function(object, mc.cores = 1) standardGeneric("doLDblock")
)

methods::setMethod(f = "doLDblock", signature = "Chromosome", definition = function(object, mc.cores = 1) {
  if (missing(object)) {
    stop('[Chromosome:doLDblock] "Chromosome" Object is required.', call. = FALSE)
  }
  data <- object@Data

  nbCORES <- mc.cores
  nbCores <- max(1, round((nbCORES - 22) / 22))
  chrLD <- object@LD

  byBlock <- split(chrLD, names(chrLD))
  byBlock <- unique(mclapply2(byBlock, mc.cores = nbCores, function(i) {
    names(i) <- NULL
    i
  }))
  LDblockTmp <- mclapply2(seq_along(byBlock), mc.cores = nbCores, function(jBlock) {
    isIn <- which(data[, "SNP"] %in% byBlock[[jBlock]])
    if (length(isIn) > 0) {
      range(data[which(data[, "SNP"] %in% byBlock[[jBlock]]), "POS"])
    } else {
      NA
    }
  })
  LDblock <- do.call("rbind", unique(LDblockTmp))
  LDblock <- stats::na.exclude(LDblock)
  rm(LDblockTmp)

  names(byBlock) <- NULL
  LDblock <- LDblock[order(LDblock[, 1]), ]
  rm(chrLD, byBlock)

  blockLim <- NULL
  for (iBlock in seq_len(nrow(LDblock))) {
    if (iBlock == 1) {
      POS <- LDblock[iBlock, ]
      blockLim <- rbind(blockLim, POS)
      jBlock <- 1
    } else {
      POS <- LDblock[iBlock, ]
      if (max(blockLim[nrow(blockLim), ]) < min(POS)) {
        blockLim <- rbind(blockLim, POS)
      } else {
        blockLim[nrow(blockLim), ] <- range(blockLim[nrow(blockLim), ], POS)
      }
      if (iBlock == nrow(LDblock)) {
        iBlock <- iBlock + 1
        if (max(blockLim[nrow(blockLim), ]) < min(POS)) {
          blockLim <- rbind(blockLim, POS)
        } else {
          blockLim[nrow(blockLim), ] <- range(blockLim[nrow(blockLim), ], POS)
        }
      }
    }
  }
  rm(LDblock)

  blockLim <- cbind(blockLim, seq_len(nrow(blockLim)))
  colnames(blockLim) <- c("MIN", "MAX", "IDBLOCK")
  rownames(blockLim) <- seq_len(nrow(blockLim))
  blockLim <- cbind(blockLim, LENGTH = NA)
  resParallel <- mclapply2(seq_len(nrow(blockLim)), mc.cores = nbCores, function(li) {
    blockLim[li, "LENGTH"] <- as.integer(blockLim[li, "MAX"]) - as.integer(blockLim[li, "MIN"])
    blockLim[li, ]
  })
  blockLim <- do.call("rbind", resParallel)
  rm(resParallel)

  data <- data[order(data[, "POS"]), ]
  data[, c("MIN", "MAX", "IDBLOCK", "LENGTH", "MAFmedian")] <- as.numeric(NA)
  tmpChr <- mclapply2(seq_len(nrow(blockLim)), mc.cores = nbCores, function(i) {
    m <- blockLim[i, ]
    interv <- seq(from = which(data[, "POS"] == m["MIN"]), to = which(data[, "POS"] == m["MAX"]))
    interv
    data[interv, c("MIN", "MAX", "IDBLOCK", "LENGTH")] <- matrix(rep(m, length(interv)), nrow = length(interv), byrow = TRUE)
    data[which(data[, "IDBLOCK"] %in% m["IDBLOCK"]), "MAFmedian"] <- stats::median(data[data[, "IDBLOCK"] %in% m["IDBLOCK"], "MAF"])
    data[interv, ]
  })

  dataTmp <- do.call("rbind", tmpChr)
  rm(tmpChr, blockLim)

  missingData <- data[!data[, "SNP"] %in% dataTmp[, "SNP"], ]
  maxIDBLOCK <- max(dataTmp[, "IDBLOCK"])
  for (iRow in seq_len(nrow(missingData))) {
    missingData[iRow, "MIN"] <- missingData[iRow, "POS"]
    missingData[iRow, "MAX"] <- missingData[iRow, "POS"]
    missingData[iRow, "LENGTH"] <- missingData[iRow, "MAX"] - missingData[iRow, "MIN"]
    missingData[iRow, "MAFmedian"] <- missingData[iRow, "MAF"]
    missingData[iRow, "IDBLOCK"] <- maxIDBLOCK + iRow
  }

  data <- rbind(missingData, dataTmp)
  data <- data[!is.na(data[, "SNP"]), ]
  data <- data[order(data[, "POS"]), ]
  rownames(data) <- data[, "SNP"]
  object@Data <- data
  object
})
