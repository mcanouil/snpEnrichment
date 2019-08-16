#' ~ Internal: snpEnrichment objects and methods ~
#'
#' These are not intended to be called by the user.
#' If you have a specific need or for more details on snpEnrichment, feel free to contact the author.
#'
#' @name internal
#' @keywords internal
NULL

#' @rdname internal
#' @keywords internal
.EnrichSNP.show <- function(object) {
  cat(
    "   - List :",
    paste0("(", length(object@List), ")"),
    if (length(object@List) <= 5) {
      if (length(object@List) == 0) "NA" else object@List
    } else {
      paste(paste(object@List[seq(5)], collapse = " "), "...")
    }
  )
  cat("\n   - Table :", paste0("(", paste(dim(object@Table), collapse = "x"), ")"))
  if (all(object@Table == 0)) {
    cat(" NA")
  } else {
    cat("\n")
    resFormat <- cbind(c("", rownames(object@Table)), rbind(colnames(object@Table), apply(round(object@Table, digits = 4), 2, as.character)))
    cat(paste("     ", apply(apply(matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)), 2, format, justify = "centre"), 1, paste, collapse = ""), "\n", sep = "", collapse = ""))
  }
  cat("\n   - EnrichmentRatio :", ifelse(length(object@EnrichmentRatio) == 0, NA, object@EnrichmentRatio))
  cat("\n   - Z               :", ifelse(length(object@Z) == 0, NA, object@Z))
  cat("\n   - PValue          :", ifelse(length(object@PValue) == 0, NA, ifelse((object@PValue == 0 & nrow(object@Resampling) != 0), paste0("<", 1 / nrow(object@Resampling)), object@PValue)))
  cat("\n   - Resampling      :", paste0("(", paste(dim(object@Resampling), collapse = "x"), ")"), ifelse(nrow(object@Resampling) == 0, 0, nrow(object@Resampling)))
  cat("\n")
  invisible(object)
}

#' @rdname internal
#' @keywords internal
.Chromosome.show <- function(object) {
  cat("  ~ Data :", paste0("(", paste(dim(object@Data), collapse = "x"), ")"))
  nrowShow <- seq_len(min(5, nrow(object@Data)))
  ncolShow <- seq_len(min(10, ncol(object@Data)))
  if (nrow(object@Data) == 0) {
    cat(" NA")
  } else {
    cat("\n")
    resFormat <- cbind(c("", rownames(object@Data[nrowShow, ncolShow])), rbind(colnames(object@Data[nrowShow, ncolShow]), apply(object@Data[nrowShow, ncolShow], 2, as.character)))
    resFormat <- rbind(resFormat, ".....")
    cat(paste("     ", apply(apply(matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)), 2, format, justify = "centre"), 1, paste, collapse = ""), "\n", sep = "", collapse = ""))
  }
  cat("\n  ~ LD :")
  if (length(object@LD) == 0) {
    cat(" NA")
  } else {
    cat("\n")
    if (length(object@LD) > 5) {
      tmpLD <- object@LD[seq_len(5)]
    } else {
      tmpLD <- object@LD
    }
    resFormat <- cbind(c("SNP1", "SNP2"), rbind(names(tmpLD), tmpLD))
    resFormat <- cbind(resFormat, "...")
    cat(paste("     ", apply(apply(matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)), 2, format, justify = "centre"), 1, paste, collapse = ""), "\n", sep = "", collapse = ""))
  }
  cat("\n  ~ eSNP :\n")
  .EnrichSNP.show(object@eSNP)
  cat("\n  ~ xSNP :\n")
  .EnrichSNP.show(object@xSNP)
  cat("\n")
  invisible(object)
}

#' @rdname internal
#' @keywords internal
.Enrichment.show <- function(object) {
  .showArgs <- function(args) {
    for (iFuncArg in seq(args)) {
      if (is.null(unlist(args[iFuncArg], use.names = FALSE))) {
        cat(paste0("    ", names(args[iFuncArg]), "() : Not yet called."))
      } else {
        tmpArgs <- args[[names(args[iFuncArg])]]
        type <- lapply(tmpArgs, class)
        res <- NULL
        tmpArgs <- lapply(tmpArgs, function(li) if (length(li) > 1) deparse(li) else li)
        for (iArg in names(tmpArgs)) {
          res <- c(
            res,
            paste0(
              iArg,
              ifelse(type[[iArg]] == "character", '="', "="),
              ifelse(
                test = type[[iArg]] == "NULL" | type[[iArg]] == "name",
                yes = deparse(tmpArgs[[iArg]]),
                no = tmpArgs[[iArg]]
              ),
              ifelse(type[[iArg]] == "character", '"', "")
            )
          )
        }
        cat(paste0(
          "    ", names(args[iFuncArg]), "(",
          paste(
            res,
            collapse = paste0(", \n", paste(rep(" ", nchar(names(args[iFuncArg])) + 5), collapse = ""))
          ),
          ")"
        ))
      }
    }
  }
  cat(" ~ Loss :", paste0("(", paste(dim(object@Loss), collapse = "x"), ")"))
  if (nrow(object@Loss) == 0) {
    cat(" NA")
  } else {
    cat("\n")
    resFormat <- cbind(
      c("", rownames(object@Loss[c(1, 2), ])),
      rbind(
        colnames(object@Loss[c(1, 2), ]),
        apply(round(object@Loss[c(1, 2), ], digits = 4), 2, as.character)
      )
    )
    resFormat <- rbind(resFormat, ".....")
    cat(paste(
      "     ",
      apply(
        X = apply(
          X = matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)),
          MARGIN = 2,
          FUN = format, justify = "centre"
        ),
        MARGIN = 1,
        FUN = paste, collapse = ""
      ),
      "\n",
      sep = "",
      collapse = ""
    ))
  }
  cat("\n ~ Call :\n")
  .showArgs(object@Call)
  cat("\n ~ eSNP :\n")
  .EnrichSNP.show(object@eSNP)
  cat("\n ~ xSNP :\n")
  .EnrichSNP.show(object@xSNP)
  nbVoid <- 0
  for (slot in names(object@Chromosomes)) {
    eval(parse(text = paste0(
      "nbVoid <- nbVoid + as.numeric(identical(chromosome(), ", paste0("object@Chromosomes$", slot), "))"
    )))
  }
  cat("\n ~ Chromosomes :", paste(nbVoid, length(names(object@Chromosomes)), sep = "/"), "empty Chromosomes")
  cat("\n")
  invisible(object)
}

#' @rdname internal
#' @keywords internal
.checkFilePath <- function(path) {
  gsub("/*$", "/", path)
}

#' @rdname internal
#' @keywords internal
.checkSignalFile <- function(signalFile) {
  if (!file.exists(signalFile)) {
    stop(paste0("[Enrichment:initFiles] ", signalFile, " doesn't exist."), call. = FALSE)
  }
  invisible()
}

#' @rdname internal
#' @keywords internal
.checkSnpInfoDir <- function(snpInfoDir) {
  snpInfoDir <- .checkFilePath(snpInfoDir)
  if (length(list.files(snpInfoDir, pattern = "*.bim")) != 22) {
    stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.bim")), " 'bim' files found when 22 is needed."), call. = FALSE)
  }
  if (length(list.files(snpInfoDir, pattern = "*.bed")) != 22) {
    stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.bed")), " 'bed' files found when 22 is needed."), call. = FALSE)
  }
  if (length(list.files(snpInfoDir, pattern = "*.fam")) != 22) {
    stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.fam")), " 'fam' files found when 22 is needed."), call. = FALSE)
  }
  invisible()
}

#' @rdname internal
#' @keywords internal
.checkSnpListDir <- function(snpListDir, pattern) {
  snpListDir <- .checkFilePath(snpListDir)
  if (length(list.files(snpListDir, pattern = pattern)) == 0) {
    stop(paste0("[Enrichment:readEnrichment] No snp list file found when at least one is needed."), call. = FALSE)
  }
  invisible()
}

#' @rdname internal
#' @keywords internal
.checkTranscript <- function(data, transcriptFile, distThresh) {
  if (any(transcriptFile != FALSE)) {
    transcript <- .readTranscript(transcriptFile = transcriptFile)
    transcriptCHR <- stats::na.exclude(transcript[transcript[, "CHR"] == unique(data[, "CHR"]), c("START", "END")])

    cisFunc <- function(line, distThresh, dataTranscript) {
      position <- as.numeric(line[4])
      CIS <- any(position > (dataTranscript[, "START"] - distThresh) & (position < (dataTranscript[, "END"] + distThresh)))
      c(line, CIS = as.numeric(CIS))
    }

    temp <- unique(as.data.frame(t(apply(data, 1, cisFunc, distThresh = distThresh * 10^3, dataTranscript = transcriptCHR)), stringsAsFactors = FALSE))
    temp[temp[, "CIS"] == 1, -grep("CIS", colnames(temp))]
  } else {
    unique(data)
  }
}

#' @rdname internal
#' @keywords internal
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
  statisticPvalue.eSNP <- stats::pnorm(Ze, lower.tail = FALSE)
  object1@eSNP <- enrichSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = c(Distribution = statisticPvalue.eSNP, Empirical = empiricPvalue.eSNP)[whichPvalue], Resampling = eSNPsample)
  if (isLD) {
    empiricPvalue.xSNP <- sum(xEnrichRatio < xSNPsample[, 5]) / length(xSNPsample[, 5])
    statisticPvalue.xSNP <- stats::pnorm(Zx, lower.tail = FALSE)
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

#' @rdname internal
#' @keywords internal
.enrichmentRatio <- function(table) {
  tmp <- as.numeric(table)
  (tmp[1] * tmp[4]) / (tmp[2] * tmp[3])
}

#' @rdname internal
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

#' @rdname internal
#' @keywords internal
.readFreq <- function(pattern, snpInfoDir) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  utils::read.delim(
    file = paste0(tmpDir, "/snpEnrichment/", pattern, ".all"), header = TRUE, stringsAsFactors = FALSE,
    colClasses = c("numeric", "character", "integer", "numeric"), na.string = c("NA", ""),
    check.names = FALSE, strip.white = TRUE, col.names = c("CHR", "SNP", "POS", "MAF")
  )
}

#' @rdname internal
#' @keywords internal
.readLD <- function(pattern, snpInfoDir, ldDir) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  if (missing(ldDir) | is.null(ldDir)) {
    IN <- paste0(tmpDir, "/snpEnrichment/", pattern, ".ld")
  } else {
    IN <- paste0(.checkFilePath(ldDir), pattern, ".ld")
  }
  nbCol <- as.character(ncol(utils::read.table(text = readLines(con = IN, n = 1), stringsAsFactors = FALSE)))
  switch(EXPR = nbCol,
    "2" = {
      utils::read.delim(
        file = IN, header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character"),
        na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, sep = "\t"
      )
    },
    "7" = {
      utils::read.delim(
        file = IN, header = TRUE, stringsAsFactors = FALSE, colClasses = c("NULL", "NULL", "character", "NULL", "NULL", "character", "NULL"),
        na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, sep = ""
      )
    },
    stop(paste0('[Enrichment:readEnrichment] "', pattern, '.ld" structure must be a matrix file with columns:\n       c("SNP_A", "SNP_B") or c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")".'), call. = FALSE)
  )
}

#' @rdname internal
#' @keywords internal
.readSignal <- function(pattern) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  utils::read.delim(
    file = paste0(tmpDir, "/snpEnrichment/", pattern, ".signal"), header = TRUE, stringsAsFactors = FALSE,
    colClasses = c("NULL", "character", "numeric"), na.string = c("NA", ""),
    check.names = FALSE, strip.white = TRUE, col.names = c("", "SNP", "PVALUE")
  )
}

#' @rdname internal
#' @keywords internal
.readSNP <- function(pattern, snpListDir) {
  snpListFile <- list.files(gsub("/$", "", snpListDir), pattern = paste0(pattern, "[^0-9]"), full.names = TRUE)
  if (is.na(snpListFile)) {
    snpList <- data.frame()
  } else {
    snpList <- utils::read.delim(
      file = snpListFile, header = FALSE, stringsAsFactors = FALSE,
      na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE
    )
    switch(EXPR = as.character(ncol(snpList)),
      "1" = colnames(snpList) <- "SNP",
      "2" = colnames(snpList) <- c("CHR", "SNP"),
      "3" = colnames(snpList) <- c("CHR", "SNP", "TRANSCRIPT"),
      colnames(snpList)[1:3] <- c("CHR", "SNP", "TRANSCRIPT")
    )
    if (ncol(snpList) >= 3) snpList <- snpList[, c("SNP", "TRANSCRIPT")]
    if (nrow(snpList) != 0) snpList[, "eSNP"] <- 1
  }
  snpList
}

#' @rdname internal
#' @keywords internal
.readTranscript <- function(transcriptFile) {
  if (all(class(try(close(file(transcriptFile)), silent = TRUE)) != "try-error")) {
    transcript <- utils::read.delim(file = transcriptFile, header = TRUE, stringsAsFactors = FALSE, na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE)
    transcript <- transcript[, 1:4]
    colnames(transcript) <- c("TRANSCRIPT", "CHR", "START", "END")
  } else {
    if (class(transcriptFile) %in% c("matrix", "data.frame")) {
      transcript <- transcriptFile[, 1:4]
      colnames(transcript) <- c("TRANSCRIPT", "CHR", "START", "END")
    } else {
      stop('[Enrichment:readEnrichment] "transcriptFile" required "matrix", "data.frame" or "txt file".', call. = FALSE)
    }
  }
  transcript
}

#' @rdname internal
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
  statisticPvalue.eSNP <- stats::pnorm(Ze, lower.tail = FALSE)
  object@eSNP <- enrichSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = c(Distribution = statisticPvalue.eSNP, Empirical = empiricPvalue.eSNP)[whichPvalue], Resampling = eSNPsample)
  if (isLD) {
    empiricPvalue.xSNP <- sum(xEnrichRatio < xSNPsample[, 5]) / length(xSNPsample[, 5])
    statisticPvalue.xSNP <- stats::pnorm(Zx, lower.tail = FALSE)
    object@xSNP <- enrichSNP(List = xSNPlist, Table = xEnrichment, EnrichmentRatio = xEnrichRatio, Z = Zx, PValue = c(Distribution = statisticPvalue.xSNP, Empirical = empiricPvalue.xSNP)[whichPvalue], Resampling = xSNPsample)
  } else {
    object@xSNP <- enrichSNP()
  }
  object
}

#' @rdname internal
#' @keywords internal
.splitByChrom <- function(pattern, snpListFile, directory) {
  if (missing(snpListFile)) {
    stop("[Enrichment:readEnrichment] argument(s) missing.", call. = FALSE)
  }
  snpList <- utils::read.delim(
    file = snpListFile, header = FALSE, stringsAsFactors = FALSE,
    na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE
  )
  switch(EXPR = as.character(ncol(snpList)),
    "1" = stop('[Enrichment:readEnrichment] at least two columns are needed: "Chromosome" and "rs name".', call. = FALSE),
    "2" = colnames(snpList) <- c("CHR", "SNP"),
    "3" = colnames(snpList) <- c("CHR", "SNP", "TRANSCRIPT"),
    colnames(snpList)[1:3] <- c("CHR", "SNP", "TRANSCRIPT")
  )
  if (ncol(snpList) > 3) {
    snpList <- snpList[, c("CHR", "SNP", "TRANSCRIPT")]
  }
  filePathDetails <- unlist(strsplit(snpListFile, "/"), use.names = FALSE)
  if (is.null(directory)) {
    filePath <- paste0(c(filePathDetails[-length(filePathDetails)], ""), collapse = "/")
  } else {
    dir.create(paste0(directory, "snpList/"), showWarnings = FALSE)
    filePath <- paste0(directory, "snpList/")
  }
  by(snpList, snpList[, "CHR"], function(snpListChr) {
    utils::write.table(snpListChr, file = paste0(filePath, pattern, unique(snpListChr[, "CHR"]), "-", filePathDetails[length(filePathDetails)]), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  })
  filePath
}

#' @rdname internal
#' @keywords internal
.verbose <- function(expr) {
  invisible(utils::capture.output(expr))
}

#' @rdname internal
#' @keywords internal
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

#' @rdname internal
#' @keywords internal
.writeSignal <- function(pattern, snpInfoDir, signalFile) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  if (length(unlist(strsplit(readLines(signalFile, n = 1), split = "\t"), use.names = FALSE)) > 1) {
    signal <- utils::read.delim(
      file = signalFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
      colClasses = c("character", "numeric"), na.string = c("NA", ""),
      check.names = FALSE, strip.white = TRUE, col.names = c("SNP", "PVALUE")
    )
  } else {
    if (length(unlist(strsplit(readLines(signalFile, n = 1), split = " "), use.names = FALSE)) > 1) {
      signal <- utils::read.delim(
        file = signalFile, header = TRUE, sep = " ", stringsAsFactors = FALSE,
        colClasses = c("character", "numeric"), na.string = c("NA", ""),
        check.names = FALSE, strip.white = TRUE, col.names = c("SNP", "PVALUE")
      )
    } else {
      stop('[Enrichment:initFiles] only " " and "\t" are allowed as columns separator in Signal file.', call. = FALSE)
    }
  }

  chrom.bim <- utils::read.delim(
    file = paste0(snpInfoDir, pattern, ".bim"), header = FALSE, stringsAsFactors = FALSE,
    colClasses = c("numeric", "character", "NULL", "NULL", "NULL", "NULL"), na.string = c("NA", ""),
    check.names = FALSE, strip.white = TRUE, col.names = c("CHR", "SNP", "", "", "", "")
  )
  signal <- merge(signal, chrom.bim, by = "SNP")[, c(3, 1, 2)]
  eval(parse(text = paste0(
    'write.table(signal, file = "', tmpDir, "/snpEnrichment/", pattern,
    '.signal", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")'
  )))
  invisible()
}

#' @rdname internal
#' @keywords internal
mclapply2 <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE, mc.allow.recursive = TRUE) {
  if (Sys.info()[["sysname"]] != "Linux") {
    mc.cores <- 1
  } else {
    mc.cores <- min(parallel::detectCores(), mc.cores)
  }
  return(parallel::mclapply(
    X = X, FUN = FUN, ...,
    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.silent = mc.silent,
    mc.cores = maxCores(mc.cores), mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
  ))
}

#' @rdname internal
#' @keywords internal
maxCores <- function(mc.cores = 1) {
  if (Sys.info()[["sysname"]] == "Linux") {
    nbCores <- parallel::detectCores()
    mc.cores.old <- mc.cores
    if (file.exists("/proc/meminfo")) {
      memInfo <- readLines("/proc/meminfo")
      sysMemFree <- memInfo[grep("^MemFree:", memInfo)]
      sysMemCached <- memInfo[grep("^Cached:", memInfo)]
      sysMemAvailable <- 0.95 * (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", sysMemFree)) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", sysMemCached)))
      sysProc <- as.numeric(unlist(strsplit(system(paste("ps v", Sys.getpid()), intern = TRUE)[2], " +"), use.names = FALSE)[8])
      mc.cores <- max(min(as.numeric(mc.cores), floor(sysMemAvailable / sysProc)), 1)
      if (mc.cores > nbCores) mc.cores <- nbCores
      if (mc.cores != mc.cores.old) {
        warning(paste0('To avoid memory overload "mc.cores" was decreased to "', mc.cores, '".'), call. = FALSE)
      }
    } else {
      mc.cores <- ifelse(mc.cores.old > nbCores, nbCores, mc.cores.old)
    }
  } else {
    mc.cores <- 1
  }
  mc.cores
}

#' @rdname internal
#' @keywords internal
GC <- function(verbose = getOption("verbose"), reset = FALSE) {
  while (!identical(gc(verbose, reset)[, 4], gc(verbose, reset)[, 4])) {}
  gc(verbose, reset)
}

#' @rdname internal
#' @keywords internal
is.EnrichSNP <- function(object) {
  if (length(object) > 1) {
    sapply(object, is.EnrichSNP)
  } else {
    class(object) == "EnrichSNP"
  }
}

#' @rdname internal
#' @keywords internal
is.chromosome <- function(object) {
  if (length(object) > 1) {
    sapply(object, is.chromosome)
  } else {
    class(object) == "Chromosome"
  }
}

#' is.enrichment
#'
#' Function to test if an object is of class [Enrichment-class].
#'
#' @param object An object to test.
#'
#' @rdname Enrichment-class
is.enrichment <-  function(object) {
  if (length(object) > 1) {
    sapply(object, is.enrichment)
  } else {
    class(object) == "Enrichment"
  }
}
