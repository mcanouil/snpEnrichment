#' Class \code{\linkS4class{Enrichment}}
#'
#' This class is defined to summarize the enrichment analysis on each
#' chromosomes and the whole genome.
#'
#'
#' @name Enrichment-class
#' @aliases Enrichment-class Enrichment [,Enrichment-method
#' [,Enrichment,ANY,ANY,ANY-method [<-,Enrichment-method
#' [<-,Enrichment,ANY,ANY,ANY-method show,Enrichment-method enrichment
#' enrichment-methods enrichment,ANY-method
#' @docType class
#' @section Objects from the Class: \code{\link{enrichment}} is defined to
#' build an object of class \code{\linkS4class{Enrichment}} in order to compute
#' an enrichment analysis.  \code{Enrichment} is the object containing the
#' results for all \code{\linkS4class{Chromosome}} object and for the whole
#' genome.\cr
#'
#' When an \code{\linkS4class{Enrichment}} object is created, it contains a
#' list of SNPs (e.g. eSNPs).  All the others slots are "empty". After
#' \code{\link{reSample}} is ran on an \code{\linkS4class{Enrichment}} object,
#' the slots: Table, EnrichmentRatio, Z, PValue and Resampling are filled.
#'
#' Note that if \code{\link{reSample}} is executed on an
#' \code{\linkS4class{Enrichment}} every new resampling is added to the
#' original ones, pre-existing statistics are erased and computed again with
#' the new resampling set.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords classes class enrichment
#' @examples
#'
#' data(toyEnrichment)
#' toyEnrich <- enrichment()
#' show(toyEnrich)
#'
#' toyEnrich["Loss"] <- toyEnrichment["Loss"]
#' toyEnrich["Loss"]
#'
#' toyEnrich <- enrichment(Loss = toyEnrichment["Loss"],
#'                         eSNP = toyEnrichment["eSNP"])
#' toyEnrich <- enrichment(Loss = toyEnrichment["Loss"])
#'
#' \dontrun{reSample(object = toyEnrichment,
#'          nSample = 10,
#'          empiricPvalue = TRUE,
#'          MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
#'          mc.cores = 1,
#'          onlyGenome = TRUE)
#' print(toyEnrichment)
#'
#' excludeFile <- c(
#' "rs7897180", "rs4725479", "rs315404", "rs17390391", "rs1650670",
#'     "rs6783390", "rs1642009", "rs4756586", "rs11995037", "rs4977345",
#'     "rs13136448", "rs4233536", "rs11151079", "rs2299657", "rs4833930",
#'     "rs1384", "rs7168184", "rs6909895", "rs7972667", "rs2293229",
#'     "rs918216", "rs6040608", "rs2817715", "rs13233541", "rs4486743",
#'     "rs2127806", "rs10912854", "rs1869052", "rs9853549", "rs448658",
#'     "rs2451583", "rs17483288", "rs10962314", "rs9612059", "rs1384182",
#'     "rs8049208", "rs12215176", "rs2980996", "rs1736976", "rs8089268",
#'     "rs10832329", "rs12446540", "rs7676237", "rs869922", "rs16823426",
#'     "rs1374393", "rs13268781", "rs11134505", "rs7325241", "rs7520109"
#' )
#'
#' toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)
#' print(toyEnrichment_exclude)}
#'
methods::setClass(
  Class = "Enrichment",
  representation = methods::representation(
    Loss = "data.frame",
    Call = "list",
    eSNP = "EnrichSNP",
    xSNP = "EnrichSNP",
    Chromosomes = "list"
  ),
  prototype = methods::prototype(
    Loss = data.frame(),
    Call = list(
      readEnrichment = list(
        pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL,
        snpInfoDir = NULL, distThresh = NULL, sigThresh = NULL, LD = NULL,
        ldDir = NULL, mc.cores = NULL
      ),
      reSample = list(
        object = NULL, nSample = NULL, empiricPvalue = NULL,
        MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL
      )
    ),
    eSNP = enrichSNP(),
    xSNP = enrichSNP(),
    Chromosomes = eval(parse(
      text = paste0(
        "list(",
        paste(paste0("Chrom", seq_len(22), " = chromosome()"), collapse = ", "),
        ")"
      )
    ))
  )
)

methods::setGeneric(
  name = "enrichment",
  def = function(Loss, Call, eSNP, xSNP, Chromosomes) standardGeneric("enrichment")
)

methods::setMethod(
  f = "enrichment",
  signature = "ANY",
  definition = function(Loss, Call, eSNP, xSNP, Chromosomes) {
    if (missing(Loss)) {
      Loss <- data.frame()
    }
    if (missing(Call)) {
      Call <- list(
        readEnrichment = list(
          pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL,
          snpInfoDir = NULL, distThresh = NULL, sigThresh = NULL, LD = NULL,
          ldDir = NULL, mc.cores = NULL
        ),
        reSample = list(
          object = NULL, nSample = NULL, empiricPvalue = NULL,
          MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL
        )
      )
    }
    if (missing(Chromosomes)) {
      Chromosomes <- eval(parse(
        text = paste0(
          "list(",
          paste(paste0("Chrom", seq_len(22), " = chromosome()"), collapse = ", "),
          ")"
        )
      ))
    }
    if (missing(eSNP)) {
      List <- eval(parse(
        text = paste0(
          "c(",
          paste(paste0("Chromosomes$Chrom", seq_len(22), "@eSNP@List"), collapse = ", "),
          ")"
        )
      ))
      eSNP <- enrichSNP(List = List)
    }
    if (missing(xSNP)) {
      List <- eval(parse(
        text = paste0(
          "c(",
          paste(paste0("Chromosomes$Chrom", seq_len(22), "@xSNP@List"), collapse = ", "),
          ")"
        )
      ))
      xSNP <- enrichSNP(List = List)
    }
    methods::new(
      "Enrichment",
      Loss = Loss,
      Call = Call,
      eSNP = eSNP,
      xSNP = xSNP,
      Chromosomes = Chromosomes
    )
  }
)

methods::setMethod(f = "[", signature = "Enrichment", definition = function(x, i, j, drop) {
  nbChr <- length(x@Chromosomes)
  if (missing(j)) {
    switch(EXPR = i,
      "Loss" = x@Loss,
      "Data" = {
        resData <- mclapply2(
          X = seq_len(22),
          mc.cores = min(22, parallel::detectCores()),
          FUN = function(iChr) x@Chromosomes[[iChr]]@Data
        )
        do.call("rbind", resData)
      },
      "LD" = {
        resLD <- mclapply2(
          X = seq_len(22),
          mc.cores = min(22, parallel::detectCores()),
          FUN = function(iChr) x@Chromosomes[[iChr]]@LD
        )
        unlist(resLD, use.names = FALSE)
      },
      "Call" = x@Call,
      "eSNP" = x@eSNP,
      "xSNP" = x@xSNP,
      "Table" = {
        res <- list(eSNP = NULL, xSNP = NULL)
        for (iType in c("eSNP", "xSNP")) {
          res[[iType]] <- eval(parse(text = paste0("x@", iType, "@Table")))
        }
        res
      },
      "Chromosomes" = x@Chromosomes,
      "Stats" = {
        res <- list(eSNP = NULL, xSNP = NULL)
        for (iType in c("eSNP", "xSNP")) {
          EnrichmentRatio <- eval(parse(text = paste0("x@", iType, "@EnrichmentRatio")))
          Z <- eval(parse(text = paste0("x@", iType, "@Z")))
          PValue <- eval(parse(text = paste0("x@", iType, "@PValue")))
          Resampling <- eval(parse(text = paste0("nrow(x@", iType, "@Resampling)")))
          Data <- x@Loss["Signal", ncol(x@Loss)]
          List <- eval(parse(text = paste0("length(x@", iType, "@List)")))
          resTmp <- c(
            if (length(EnrichmentRatio) == 0) NA else EnrichmentRatio,
            if (length(Z) == 0) NA else Z,
            if (length(PValue) == 0) NA else PValue,
            if (length(Resampling) == 0) NA else Resampling,
            if (length(Data) == 0) NA else Data,
            if (length(List) == 0) NA else List
          )
          names(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", iType)
          res[[iType]] <- resTmp
        }
        res
      },
      stop("[Enrichment:get] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
    )
  } else {
    if (max(j) > nbChr) {
      if (j == "ALL") {
        switch(EXPR = i,
          "Loss" = x@Loss,
          "Data" = {
            resData <- mclapply2(
              X = seq_len(22),
              mc.cores = min(22, parallel::detectCores()),
              FUN = function(iChr) x@Chromosomes[[iChr]]@Data
            )
            do.call("rbind", resData)
          },
          "LD" = {
            resLD <- mclapply2(
              X = seq_len(22),
              mc.cores = min(22, parallel::detectCores()),
              FUN = function(iChr) x@Chromosomes[[iChr]]@LD
            )
            unlist(resLD, use.names = FALSE)
          },
          "Call" = x@Call,
          "List" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@", iType, "@List")))
            }
            res
          },
          "Table" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@", iType, "@Table")))
            }
            res
          },
          "EnrichmentRatio" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@", iType, "@EnrichmentRatio")))
            }
            res
          },
          "Z" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@", iType, "@Z")))
            }
            res
          },
          "PValue" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@", iType, "@PValue")))
            }
            res
          },
          "Resampling" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@", iType, "@Resampling")))
            }
            res
          },
          "Chromosomes" = x@Chromosomes,
          "Stats" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              EnrichmentRatio <- eval(parse(text = paste0("x@", iType, "@EnrichmentRatio")))
              Z <- eval(parse(text = paste0("x@", iType, "@Z")))
              PValue <- eval(parse(text = paste0("x@", iType, "@PValue")))
              Resampling <- eval(parse(text = paste0("nrow(x@", iType, "@Resampling)")))
              Data <- x@Loss["Signal", ncol(x@Loss)]
              List <- eval(parse(text = paste0("length(x@", iType, "@List)")))
              resTmp <- c(
                if (length(EnrichmentRatio) == 0) NA else EnrichmentRatio,
                if (length(Z) == 0) NA else Z,
                if (length(PValue) == 0) NA else PValue,
                if (length(Resampling) == 0) NA else Resampling,
                if (length(Data) == 0) NA else Data,
                if (length(List) == 0) NA else List
              )
              for (iChr in seq_len(22)) {
                EnrichmentRatio <- eval(parse(text = paste0("x@Chromosomes$Chrom", iChr, "@", iType, "@EnrichmentRatio")))
                Z <- eval(parse(text = paste0("x@Chromosomes$Chrom", iChr, "@", iType, "@Z")))
                PValue <- eval(parse(text = paste0("x@Chromosomes$Chrom", iChr, "@", iType, "@PValue", )))
                Resampling <- eval(parse(text = paste0("nrow(x@Chromosomes$Chrom", iChr, "@", iType, "@Resampling)")))
                Data <- eval(parse(text = paste0("nrow(x@Chromosomes$Chrom", iChr, "@Data)")))
                List <- eval(parse(text = paste0("length(x@Chromosomes$Chrom", iChr, "@", iType, "@List)")))
                tmp <- c(
                  if (length(EnrichmentRatio) == 0) NA else EnrichmentRatio,
                  if (length(Z) == 0) NA else Z,
                  if (length(PValue) == 0) NA else PValue,
                  if (length(Resampling) == 0) NA else Resampling,
                  if (length(Data) == 0) NA else Data,
                  if (length(List) == 0) NA else List
                )
                resTmp <- rbind(resTmp, tmp)
              }
              rownames(resTmp) <- c("Genome", paste0("Chrom", seq_len(22)))
              colnames(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", iType)
              res[[iType]] <- resTmp
            }
            res
          },
          stop("[Enrichment:get] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
        )
      } else {
        stop('[Enrichment:get] "j" is out of limits.', call. = FALSE)
      }
    } else {
      if (length(j) > 1) {
        switch(EXPR = i,
          "Signal" = {
            eval(parse(text = paste0("list(", paste(paste0("x@Chromosomes$Chrom", j, "@Data[, c('SNP', 'PVALUE')]"), collapse = ", "), ")")))
          },
          "Data" = {
            resData <- mclapply2(
              X = j,
              mc.cores = min(length(j), parallel::detectCores()),
              FUN = function(iChr) x@Chromosomes[[iChr]]@Data
            )
            do.call("rbind", resData)
          },
          "LD" = {
            resLD <- mclapply2(
              X = j,
              mc.cores = min(length(j), parallel::detectCores()),
              FUN = function(iChr) x@Chromosomes[[iChr]]@LD
            )
            unlist(resLD, use.names = FALSE)
          },
          "Call" = x@Call,
          "List" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("c(", paste(paste0("x@Chromosomes$Chrom", j, "@", iType, "@List"), collapse = ", "), ")")))
            }
            res
          },
          "Table" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("list(", paste(paste0("x@Chromosomes$Chrom", j, "@", iType, "@Table"), collapse = ", "), ")")))
            }
            res
          },
          "EnrichmentRatio" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("c(", paste(paste0("x@Chromosomes$Chrom", j, "@", iType, "@EnrichmentRatio"), collapse = ", "), ")")))
            }
            res
          },
          "Z" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("c(", paste(paste0("x@Chromosomes$Chrom", j, "@", iType, "@Z"), collapse = ", "), ")")))
            }
            res
          },
          "PValue" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("c(", paste(paste0("x@Chromosomes$Chrom", j, "@", iType, "@PValue"), collapse = ", "), ")")))
            }
            res
          },
          "Resampling" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("list(", paste(paste0("x@Chromosomes$Chrom", j, "@", iType, "@Resampling"), collapse = ", "), ")")))
            }
            res
          },
          "Chromosomes" = {
            x@Chromosomes[j]
          },
          "Stats" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              resTmp <- NULL
              for (iChr in j) {
                EnrichmentRatio <- eval(parse(text = paste0("x@Chromosomes$Chrom", iChr, "@", iType, "@EnrichmentRatio")))
                Z <- eval(parse(text = paste0("x@Chromosomes$Chrom", iChr, "@", iType, "@Z")))
                PValue <- eval(parse(text = paste0("x@Chromosomes$Chrom", iChr, "@", iType, "@PValue")))
                Resampling <- eval(parse(text = paste0("nrow(x@Chromosomes$Chrom", iChr, "@", iType, "@Resampling)")))
                Data <- eval(parse(text = paste0("nrow(x@Chromosomes$Chrom", iChr, "@Data)")))
                List <- eval(parse(text = paste0("length(x@Chromosomes$Chrom", iChr, "@", iType, "@List)")))
                tmp <- c(
                  if (length(EnrichmentRatio) == 0) {
                    NA
                  } else {
                    EnrichmentRatio
                  },
                  if (length(Z) == 0) NA else Z,
                  if (length(PValue) == 0) NA else PValue,
                  if (length(Resampling) == 0) NA else Resampling,
                  if (length(Data) == 0) NA else Data,
                  if (length(List) == 0) NA else List
                )
                resTmp <- rbind(resTmp, tmp)
              }
              rownames(resTmp) <- c(paste0("Chrom", j))
              colnames(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", iType)
              res[[iType]] <- resTmp
            }
            res
          },
          stop("[Enrichment:get] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
        )
      } else {
        switch(EXPR = i,
          "Signal" = eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@Data[, c('SNP', 'PVALUE')]"))),
          "Data" = eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@Data"))),
          "LD" = eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@LD"))),
          "Call" = x@Call,
          "List" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", iType, "@List")))
            }
            res
          },
          "Table" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", iType, "@Table")))
            }
            res
          },
          "EnrichmentRatio" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", iType, "@EnrichmentRatio")))
            }
            res
          },
          "Z" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", iType, "@Z")))
            }
            return(res)
          },
          "PValue" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", iType, "@PValue")))
            }
            res
          },
          "Resampling" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              res[[iType]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", iType, "@Resampling")))
            }
            res
          },
          "Chromosomes" = x@Chromosomes[[j]],
          "Stats" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (iType in c("eSNP", "xSNP")) {
              EnrichmentRatio <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", iType, "@EnrichmentRatio")))
              Z <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", iType, "@Z")))
              PValue <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", iType, "@PValue")))
              Resampling <- eval(parse(text = paste0("nrow(x@Chromosomes$Chrom", j, "@", iType, "@Resampling)")))
              Data <- eval(parse(text = paste0("nrow(x@Chromosomes$Chrom", j, "@Data)")))
              List <- eval(parse(text = paste0("length(x@Chromosomes$Chrom", j, "@", iType, "@List)")))
              tmp <- c(
                if (length(EnrichmentRatio) == 0) NA else EnrichmentRatio,
                if (length(Z) == 0) NA else Z,
                if (length(PValue) == 0) NA else PValue,
                if (length(Resampling) == 0) NA else Resampling,
                if (length(Data) == 0) NA else Data,
                if (length(List) == 0) NA else List
              )
              names(tmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", iType)
              res[[iType]] <- tmp
            }
            res
          },
          stop("[Enrichment:get] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
        )
      }
    }
  }
})

methods::setMethod(f = "[<-", signature = "Enrichment", definition = function(x, i, j, value) {
  nbChr <- length(x@Chromosomes)
  if (missing(j)) {
    switch(EXPR = i,
      "Loss" = x@Loss <- value,
      "Data" = stop('[Enrichment:set] "Data" is not available for Set function.', call. = FALSE),
      "LD" = stop('[Enrichment:set] "LD" is not available for Set function.', call. = FALSE),
      "Call" = x@Call <- value,
      "eSNP" = x@eSNP <- value,
      "xSNP" = x@xSNP <- value,
      "List" = stop('[Enrichment:set] "List" is not available for Set function.', call. = FALSE),
      "Table" = stop('[Enrichment:set] "Table" is not available for Set function.', call. = FALSE),
      "EnrichmentRatio" = stop('[Enrichment:set] "EnrichmentRatio" is not available for Set function.', call. = FALSE),
      "Z" = stop('[Enrichment:set] "Z" is not available for Set function.', call. = FALSE),
      "PValue" = stop('[Enrichment:set] "PValue" is not available for Set function.', call. = FALSE),
      "Resampling" = stop('[Enrichment:set] "Resampling" is not available for Set function.', call. = FALSE),
      "Chromosomes" = x@Chromosomes <- value,
      "Stats" = stop('[Enrichment:set] "Stats" is not available for Set function.', call. = FALSE),
      stop("[Enrichment:set] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
    )
  } else {
    if (max(j) > nbChr) {
      stop('[Enrichment:set] "j" is out of limits.', call. = FALSE)
    } else {
      if (length(j) > 1) {
        stop('[Enrichment:set] "j" must be atomic.', call. = FALSE)
      } else {
        switch(EXPR = i,
          "Loss" = x@Loss <- value,
          "Data" = stop('[Enrichment:set] "Data" is not available for Set function.', call. = FALSE),
          "LD" = stop('[Enrichment:set] "LD" is not available for Set function.', call. = FALSE),
          "Call" = x@Call <- value,
          "eSNP" = x@Chromosomes[[j]]@eSNP <- value,
          "xSNP" = x@Chromosomes[[j]]@xSNP <- value,
          "List" = stop('[Enrichment:set] "List" is not available for Set function.', call. = FALSE),
          "Table" = stop('[Enrichment:set] "Table" is not available for Set function.', call. = FALSE),
          "EnrichmentRatio" = stop('[Enrichment:set] "EnrichmentRatio" is not available for Set function.', call. = FALSE),
          "Z" = stop('[Enrichment:set] "Z" is not available for Set function.', call. = FALSE),
          "PValue" = stop('[Enrichment:set] "PValue" is not available for Set function.', call. = FALSE),
          "Resampling" = stop('[Enrichment:set] "Resampling" is not available for Set function.', call. = FALSE),
          "Chromosomes" = x@Chromosomes[[j]] <- value,
          "Stats" = stop('[Enrichment:set] "Stats" is not available for Set functions.', call. = FALSE),
          stop("[Enrichment:set] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
        )
      }
    }
  }
  methods::validObject(x)
  x
})
