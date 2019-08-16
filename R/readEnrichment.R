#' Read and create EnrichmentRatio object
#'
#' Read files created by [initFiles] and create an [Enrichment-class] object.
#'
#' @param pattern A character string containing a expression to be matched with all chromosomes files
#'   (*e.g.*,`"Chrom"` for files which start by `"Chrom"` followed by the chromosome number).
#' @param signalFile The name of the signal file which the data are to be read from
#'   (2 columns: "SNP" and "PVALUE"). Each row of the table appears as one line of the file.
#'   If it does not contain an `_absolute_` path, the file name is `_relative_` to the current
#'   working directory, [getwd]. The fields separator character have to be a space `" "`
#'   or a tabulation `"\t"`.
#' @param transcriptFile A character string naming a file or a [data.frame] with four columns:
#'   Chromomosome, trancript's name, Starting and Ending positions.
#'   `data(trancript)` can be use as parameters. Default is `FALSE`.
#' @param snpListDir [character]: character string naming a directory
#'   containing a list of SNPs for one or several chromosomes. `snpListDir`
#'   can be a single file with at least two columns: chromosome and rs name.
#' @param snpInfoDir [character]: character string naming a directory
#'   containing the reference data in a PLINK format (.bed, .bim and .fam).
#' @param distThresh [numeric]: maximal distance (kb) between SNP and gene.
#'   `distThresh` is used if `transcriptFile` is set.
#' @param sigThresh [numeric]: statistical threshold for signal (*e.g.*, `0.05` for a given GWAS signal)
#'   used to compute an Enrichment Ratio.
#' @param LD [logical]: `LD=TRUE` (default is `FALSE`) read LD compute with [writeLD] function
#'   or with PLINK. Note that, this setting can increase the computation's time,
#'   depending on number of SNPs in the signal file.
#' @param ldDir [character]: character string naming a directory where the linkage disequilibrium files
#'   should be read (default `NULL` is in temporary directory). LD files can be the LD output from plink.
#' @param mc.cores [numeric]: the number of cores to use (default is `1`),
#'   *i.e.*, at most how many child processes will be run simultaneously.
#'   Must be at least one, and parallelization requires at least two cores.
#'
#' @return Return an object of class [Enrichment-class] partly filled.
#'
#' @examples
#' if (interactive()) {
#'   snpListDir <- system.file("extdata/List", package = "snpEnrichment")
#'   signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment")
#'   snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment")
#'   data(transcript)
#'   transcriptFile <- transcript
#'
#'   initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)
#'   toyData <- readEnrichment(
#'     pattern = "Chrom",
#'     signalFile,
#'     transcriptFile,
#'     snpListDir,
#'     snpInfoDir,
#'     distThresh = 1000,
#'     sigThresh = 0.05,
#'     LD = FALSE,
#'     ldDir = NULL,
#'     mc.cores = 1
#'   )
#'   toyData
#' }
#'
#' @export
readEnrichment <- function(pattern = "Chrom", signalFile, transcriptFile = FALSE, snpListDir, snpInfoDir, distThresh = 1000, sigThresh = 0.05, LD = FALSE, ldDir = NULL, mc.cores = 1) {
  cat("############# Read Enrichment ##############\n")
  if (missing(signalFile) | missing(snpListDir) | missing(snpInfoDir)) {
    stop("[Enrichment:readEnrichment] argument(s) missing.", call. = FALSE)
  }
  tmpDir <- paste0(gsub("\\\\", "/", tempdir()), "/snpEnrichment/")
  if (length(list.files(tmpDir, pattern = "\\.signal")) != 22 & length(list.files(tmpDir, pattern = "\\.all")) != 22) {
    stop('[Enrichment:readEnrichment] "initFiles" must be run before.', call. = FALSE)
  }

  if (is.null(ldDir)) {
    if (LD & length(list.files(tmpDir, pattern = "\\.ld")) != 22) {
      stop('[Enrichment:readEnrichment] "writeLD" must be run before. (Or LD computation by PLINK.)', call. = FALSE)
    }
  } else {
    if (LD & length(list.files(ldDir, pattern = "\\.ld")) != 22) {
      stop('[Enrichment:readEnrichment] "writeLD" must be run before. (Or LD computation by PLINK.)', call. = FALSE)
    }
  }

  snpInfoDir <- .checkFilePath(snpInfoDir)
  .checkSnpInfoDir(snpInfoDir)
  .checkSignalFile(signalFile)

  if (length(list.files(snpListDir)) == 0) {
    snpListDir <- .splitByChrom(pattern = pattern, snpListFile = snpListDir, directory = tmpDir)
  } else {
    snpListDir <- .checkFilePath(snpListDir)
    .checkSnpListDir(snpListDir, pattern)
  }

  cat("  Read Chromosomes:\n    ")
  resParallel <- mclapply2(seq_len(22), mc.cores = min(22, mc.cores), FUN = function(iChr) {
    files <- .readFiles(pattern = paste0(pattern, iChr), snpInfoDir = snpInfoDir, snpListDir = snpListDir, distThresh = distThresh)
    if (LD) {
      newPattern <- gsub(".bim", "", grep(paste0(pattern, iChr, "[^0-9]"), list.files(snpInfoDir, pattern = ".bim"), value = TRUE))
      linkageData <- .readLD(pattern = newPattern, snpInfoDir = snpInfoDir, ldDir = ldDir)

      data <- files$data[order(files$data$POS), ]
      snpSignal <- data[, "SNP"]

      ldData <- linkageData[, 2]
      names(ldData) <- linkageData[, 1]
      eSNP <- data[data[, "eSNP"] == 1, "SNP"]
      ldSNP <- linkageData[linkageData[, "SNP_A"] %in% eSNP, "SNP_B"]

      xSNP <- unique(c(eSNP, ldSNP))
      data[data[, "SNP"] %in% xSNP, "xSNP"] <- 1
      data <- data[!is.na(data[, "PVALUE"]), ]
      resChr <- chromosome(Data = data, LD = ldData)
    } else {
      resChr <- chromosome(Data = files$data)
    }
    snpLoss <- files$snpLoss

    if (all(transcriptFile == FALSE)) {
      resCheckData <- resChr@Data
      resCheckData <- transform(unique(resCheckData),
        SNP = as.character(resCheckData$SNP),
        PVALUE = as.numeric(resCheckData$PVALUE),
        CHR = as.numeric(resCheckData$CHR),
        POS = as.numeric(resCheckData$POS),
        MAF = as.numeric(resCheckData$MAF),
        eSNP = as.numeric(resCheckData$eSNP),
        xSNP = as.numeric(resCheckData$xSNP)
      )
      signalLoss <- c(NA, length(unique(resChr@Data[, "SNP"])), length(unique(resChr@Data[!is.na(resChr@Data[, "PVALUE"]), "SNP"])))
    } else {
      resCheckData <- .checkTranscript(data = resChr@Data, transcriptFile = transcriptFile, distThresh = distThresh)
      snpLoss <- c(snpLoss, sum(as.numeric(resCheckData[, "eSNP"])))
      resCheckData <- transform(unique(resCheckData),
        SNP = as.character(resCheckData$SNP),
        PVALUE = as.numeric(resCheckData$PVALUE),
        CHR = as.numeric(resCheckData$CHR),
        POS = as.numeric(resCheckData$POS),
        MAF = as.numeric(resCheckData$MAF),
        eSNP = as.numeric(resCheckData$eSNP),
        xSNP = as.numeric(resCheckData$xSNP)
      )
      signalLoss <- c(NA, length(unique(resChr@Data[, "SNP"])), length(unique(resChr@Data[!is.na(resChr@Data[, "PVALUE"]), "SNP"])), length(unique(resCheckData[!is.na(resCheckData[, "PVALUE"]), "SNP"])))
    }

    resChr <- chromosome(Data = resCheckData, LD = resChr@LD)

    cat(iChr, " ", sep = "")
    list(resChr, snpLoss, signalLoss)
  })
  names(resParallel) <- paste0("Chrom", seq_len(22))

  result <- enrichment()
  result@Chromosomes <- lapply(resParallel, "[[", 1)
  snpLoss <- t(sapply(resParallel, "[[", 2))
  signalLoss <- rowSums(sapply(resParallel, "[[", 3))
  loss <- rbind(Signal = signalLoss, Genome = apply(snpLoss, 2, sum), snpLoss)
  colnames(loss) <- c("Rows", "Unique", "Intersect.Signal", "CIS")[seq_len(ncol(loss))]

  SNPs <- result["List", seq_len(22)]
  result@eSNP@List <- SNPs[["eSNP"]]
  result@xSNP@List <- SNPs[["xSNP"]]

  if (length(unlist(strsplit(readLines(signalFile, n = 1), split = "\t"), use.names = FALSE)) > 1) {
    signal <- utils::read.delim(file = signalFile, header = TRUE, sep = "\t", colClasses = c("character", "numeric"), na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, col.names = c("SNP", "PVALUE"), stringsAsFactors = FALSE)
  } else {
    if (length(unlist(strsplit(readLines(signalFile, n = 1), split = " "), use.names = FALSE)) > 1) {
      signal <- utils::read.delim(file = signalFile, header = TRUE, sep = " ", colClasses = c("character", "numeric"), na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, col.names = c("SNP", "PVALUE"), stringsAsFactors = FALSE)
    } else {
      stop('[Enrichment:readEnrichment] only " " and "\t" are allowed as columns separator in Signal file.', call. = FALSE)
    }
  }
  loss[1, c(1, 2)] <- c(length(signal[, "SNP"]), length(unique(signal[, "SNP"])))
  result@Loss <- as.data.frame(loss)

  sysCall <- sys.call(sys.parent())
  argsSNP <- as.list(sysCall[-1])
  formal <- as.list(names(formals(as.character(sysCall))))
  names(formal) <- formal
  if (is.null(names(argsSNP))) {
    names(argsSNP) <- names(formal)[seq_along(argsSNP)]
  } else {
    emptyNames <- which(names(argsSNP) == "")
    names(argsSNP)[emptyNames] <- names(formal)[emptyNames]
  }
  names(argsSNP)[grep("^$", names(argsSNP))] <- names(formal)[grep("^$", names(argsSNP))]
  argsSNP <- c(argsSNP, lapply(formal[!names(formal) %in% names(argsSNP)], as.name))[names(formal)]
  paramClass <- sapply(argsSNP, class)

  for (iArg in names(formal)) {
    if (iArg != "...") {
      paramPos <- grep(iArg, names(formal), fixed = TRUE)
      argTmp <- argsSNP[[paramPos]]
      classTmp <- paramClass[[paramPos]]
      switch(EXPR = classTmp,
        "character" = formal[[iArg]] <- argTmp,
        "logical" = formal[[iArg]] <- argTmp,
        "numeric" = formal[[iArg]] <- argTmp,
        "integer" = formal[[iArg]] <- argTmp,
        "NULL" = formal[[iArg]] <- "NULL",
        "call" = formal[[iArg]] <- eval(argTmp),
        "name" = {
          resEval <- eval(argTmp)
          switch(EXPR = class(resEval),
            "character" = formal[[iArg]] <- resEval,
            "logical" = formal[[iArg]] <- resEval,
            "numeric" = formal[[iArg]] <- resEval,
            "integer" = formal[[iArg]] <- resEval,
            "matrix" = formal[[iArg]] <- argTmp,
            "data.frame" = formal[[iArg]] <- argTmp,
            "Enrichment" = formal[[iArg]] <- argTmp
          )
        }
      )
    }
  }
  formal[["ldDir"]] <- if (is.null(ldDir)) tmpDir else .checkFilePath(formal[["ldDir"]])
  result@Call$readEnrichment <- formal

  cat("\n########### Read Enrichment Done ###########\n\n")
  computeER(object = result, sigThresh = sigThresh, mc.cores = mc.cores)
}
