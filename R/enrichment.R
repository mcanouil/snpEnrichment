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


methods::setMethod(f = "is.enrichment", signature = "ANY", definition = function(object) {
  if (length(object) > 1) {
    sapply(object, is.enrichment)
  } else {
    class(object) == "Enrichment"
  }
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
methods::setMethod(f = "show", signature = "Enrichment", definition = function(object) {
  cat("    ~~~ Class:", class(object), "~~~ \n")
  .Enrichment.show(object)
  invisible(object)
})


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


methods::setMethod(f = "computeER", signature = "Enrichment", definition = function(object, sigThresh = 0.05, mc.cores = 1) {
  if (missing(object)) {
    stop('[Enrichment:computeER] "Enrichment" object is required.', call. = FALSE)
  }
  rowNames <- c(paste0("P>=", sigThresh), paste0("P<", sigThresh))
  object@Chromosomes <- mclapply2(object@Chromosomes, mc.cores = mc.cores, function(chr) {
    data <- chr@Data
    chrLD <- length(chr@LD)
    pvalFactor <- factor(data[, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
    for (iType in c("eSNP", "xSNP")) {
      if (!(chrLD == 0 & iType == "xSNP")) {
        snpEnrich <- table(pvalFactor, factor(data[, iType], levels = c(0, 1)))
        colnames(snpEnrich) <- c("otherSNP", iType)
        rownames(snpEnrich) <- rowNames
        chr[iType] <- enrichSNP(List = chr[iType]@List, Table = unclass(snpEnrich), EnrichmentRatio = .enrichmentRatio(snpEnrich))
      }
    }
    chr
  })
  for (iType in c("eSNP", "xSNP")) {
    bigEnrichment <- matrix(rowSums(sapply(seq_len(22), function(jChr) {
      object@Chromosomes[[jChr]][iType]@Table
    })), nrow = 2, ncol = 2)
    object[iType]@Table <- bigEnrichment
    object[iType]@EnrichmentRatio <- .enrichmentRatio(bigEnrichment)
  }
  object
})


methods::setMethod(f = "reSample", signature = "Enrichment", definition = function(object, nSample = 100, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = TRUE) {
  if (missing(object)) {
    stop('[Enrichment:reSample] "Enrichment" object is required.', call. = FALSE)
  }
  if (nSample < 10) {
    nSample <- 10
    warning("[Enrichment:reSample] nSample was increased to 10.", call. = FALSE)
  }
  sigThresh <- object@Call$readEnrichment$sigThresh

  cat("########### Resample Enrichment ############\n")
  warnings.env <- new.env()
  assign("minCores", mc.cores, envir = warnings.env)
  assign("maxCores", 0, envir = warnings.env)
  nSampleOld <- object@Call$reSample$nSample
  if (onlyGenome == FALSE) {
    listRes <- eval(parse(text = paste0("list(", paste(paste0("Chrom", seq_len(22), " = NULL"), collapse = ", "), ")")))
    for (iChr in seq_len(22)) {
      cat("  Chromosome ", if (nchar(iChr) == 1) paste0("0", iChr) else iChr, ": ", sep = "")
      nbCores <- suppressWarnings(maxCores(mc.cores))
      assign("minCores", min(get("minCores", envir = warnings.env), nbCores), envir = warnings.env)
      assign("maxCores", max(get("maxCores", envir = warnings.env), nbCores), envir = warnings.env)
      suppressWarnings(listRes[[iChr]] <- reSample(object = object@Chromosomes[[iChr]], nSample = nSample, empiricPvalue = empiricPvalue, sigThresh = sigThresh, MAFpool = MAFpool, mc.cores = mc.cores))
      cat("END\n")
    }
    object@Chromosomes <- listRes
    rm(listRes)
  }

  cat("  Genome       : ")
  nbCores <- suppressWarnings(maxCores(mc.cores))
  assign("minCores", min(get("minCores", envir = warnings.env), nbCores), envir = warnings.env)
  assign("maxCores", max(get("maxCores", envir = warnings.env), nbCores), envir = warnings.env)
  suppressWarnings(result <- .reSample(object = object, nSample = nSample, empiricPvalue = empiricPvalue, sigThresh = sigThresh, MAFpool = MAFpool, mc.cores = mc.cores))
  cat("END\n")
  rm(object)

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
          if (class(try(resEval <- eval(argTmp), silent = TRUE)) == "try-error") {
            formal[[iArg]] <- argTmp
          } else {
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
        }
      )
    }
  }

  if (is.numeric(nSampleOld)) {
    formal$nSample <- nSampleOld + formal$nSample
  }
  result@Call$reSample <- formal[c("object", "nSample", "empiricPvalue", "MAFpool", "mc.cores", "onlyGenome")]
  nameObject <- deparse(result@Call$reSample[["object"]])
  assign(nameObject, result, inherits = TRUE, envir = parent.frame(2))

  assign("maxCores", min(get("maxCores", envir = warnings.env), mc.cores), envir = warnings.env)
  if (get("minCores", envir = warnings.env) == get("maxCores", envir = warnings.env)) {
    if (get("minCores", envir = warnings.env) != mc.cores) {
      warning(paste0(
        '[Enrichment:reSample] To avoid memory overload "mc.cores" was decreased to ',
        get("minCores", envir = warnings.env), "."
      ), call. = FALSE)
    }
  } else {
    warning(paste0(
      '[Enrichment:reSample] To avoid memory overload "mc.cores" was decreased to min=',
      get("minCores", envir = warnings.env), " and max=",
      get("maxCores", envir = warnings.env), "."
    ), call. = FALSE)
  }
  cat("######## Resample Enrichment Done ##########\n")
  cat(paste0('*** Object "', nameObject, '" has been updated. ***\n\n'))
  invisible(result)
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


methods::setMethod(f = "reset", signature = "Enrichment", definition = function(object, i) {
  switch(EXPR = i,
    "Loss" = object@Loss <- data.frame(),
    "Data" = object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Data"),
    "LD" = object@Chromosomes <- lapply(object@Chromosomes, reset, i = "LD"),
    "Call" = {
      object@Call <- list(
        readEnrichment = list(
          pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL, snpInfoDir = NULL,
          distThresh = NULL, sigThresh = NULL, LD = NULL, ldDir = NULL, mc.cores = NULL
        ),
        reSample = list(
          object = NULL, nSample = NULL, empiricPvalue = NULL,
          MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL
        )
      )
    },
    "eSNP" = {
      object@eSNP <- enrichSNP()
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "eSNP")
    },
    "xSNP" = {
      object@xSNP <- enrichSNP()
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "xSNP")
    },
    "List" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@List <- character()
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "List")
    },
    "Table" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@Table <- matrix(0, ncol = 2, nrow = 2)
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Table")
    },
    "EnrichmentRatio" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@EnrichmentRatio <- numeric()
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "EnrichmentRatio")
    },
    "Z" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@Z <- numeric()
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Z")
    },
    "PValue" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@PValue <- numeric()
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "PValue")
    },
    "Resampling" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@Resampling <- matrix(0, ncol = 5, nrow = 0)
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Resampling")
    },
    "Chromosomes" = {
      object@Chromosomes <- eval(parse(text = paste0("list(", paste(paste0("Chrom", seq_len(22), " = chromosome()"), collapse = ", "), ")")))
    },
    stop("[Enrichment:reset] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
  )
  object
})


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


methods::setMethod(f = "plot", signature = "Enrichment", definition = function(x, what = "Genome", type = c("eSNP", "xSNP"), ggplot = FALSE, pvalue = TRUE, ...) {
  if (is.null(unlist(x@Call$reSample, use.names = FALSE))) {
    stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
  }

  if (is.null(what) | any(!what %in% c("All", "Genome", seq_len(22)))) {
    stop('[Enrichment:plot] "what" must be: "All", "Genome" or numeric value (atomic or vector).', call. = FALSE)
  }
  if (is.null(type) | any(!type %in% c("eSNP", "xSNP"))) {
    stop('[Enrichment:plot] "type" must be: "eSNP" and/or "xSNP".', call. = FALSE)
  }
  if (any(type %in% "xSNP") & length(x["xSNP"]["List"]) == 0) {
    type <- "eSNP"
  }

  .computeER4plot <- function(EnrichSNPObject) {
    ER <- EnrichSNPObject@EnrichmentRatio
    if (nrow(EnrichSNPObject@Resampling) == 0) {
      stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
    }
    resampling <- EnrichSNPObject@Resampling[, 5]
    ERsample <- NULL
    size <- length(resampling)
    if (size >= 1000) {
      interv <- unique(c(seq(from = min(1000, floor(0.1 * size)), to = size, by = floor(size / 1000)), size))
    } else {
      interv <- unique(c(seq(from = max(floor(0.1 * size), 3), to = size, by = 1), size))
    }
    ERsample <- sapply(interv, function(k) {
      resamplingInterv <- resampling[1:k]
      resamplingClean <- resamplingInterv[!(is.na(resamplingInterv) | is.infinite(resamplingInterv))]
      mu <- mean(resamplingClean)
      sigma <- sqrt(stats::var(resamplingClean))

      if (sigma == 0 | is.na(sigma)) {
        if (mu == 0) 0 else ER - mu
      } else {
        (ER - mu) / sigma
      }
    })
    names(ERsample) <- interv
    as.matrix(ERsample)
  }
  .computeEmpP4plot <- function(EnrichSNPObject) {
    ER <- EnrichSNPObject@EnrichmentRatio
    if (nrow(EnrichSNPObject@Resampling) == 0) {
      stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
    }
    resampling <- EnrichSNPObject@Resampling[, 5]
    ERsample <- NULL
    size <- length(resampling)
    if (size >= 1000) {
      interv <- unique(c(seq(from = min(1000, floor(0.1 * size)), to = size, by = floor(size / 1000)), size))
    } else {
      interv <- unique(c(seq(from = max(floor(0.1 * size), 3), to = size, by = 1), size))
    }
    ERsample <- sapply(interv, function(k) {
      resamplingInterv <- resampling[1:k]
      resamplingClean <- resamplingInterv[!(is.na(resamplingInterv) | is.infinite(resamplingInterv))]
      sum(EnrichSNPObject@EnrichmentRatio < resamplingClean) / length(resamplingClean)
    })
    names(ERsample) <- interv
    as.matrix(ERsample)
  }

  if (x@Call$reSample$empiricPvalue) {
    matrixER <- list(eSNP = NULL, xSNP = NULL)
    for (iType in type) {
      if (length(what) == 1) {
        switch(EXPR = as.character(what),
          "Genome" = {
            matrixER[[iType]] <- .computeEmpP4plot(x[iType])
            colnames(matrixER[[iType]]) <- "Genome"
          },
          "All" = {
            matrixER[[iType]] <- .computeEmpP4plot(x[iType])
            for (j in seq_len(22)) {
              matrixER[[iType]] <- cbind(matrixER[[iType]], .computeEmpP4plot(x["Chromosomes", j][iType]))
            }
            colnames(matrixER[[iType]]) <- c("Genome", paste0("Chrom", seq_len(22)))
          }, {
            for (j in what) {
              matrixER[[iType]] <- cbind(matrixER[[iType]], .computeEmpP4plot(x["Chromosomes", j][iType]))
            }
            colnames(matrixER[[iType]]) <- paste0("Chrom", what)
          }
        )
      } else {
        for (j in what) {
          matrixER[[iType]] <- cbind(matrixER[[iType]], .computeEmpP4plot(x["Chromosomes", j][iType]))
        }
        colnames(matrixER[[iType]]) <- paste0("Chrom", what)
      }
    }
  } else {
    matrixER <- list(eSNP = NULL, xSNP = NULL)
    for (iType in type) {
      if (length(what) == 1) {
        switch(EXPR = as.character(what),
          "Genome" = {
            matrixER[[iType]] <- .computeER4plot(x[iType])
            colnames(matrixER[[iType]]) <- "Genome"
          },
          "All" = {
            matrixER[[iType]] <- .computeER4plot(x[iType])
            for (j in seq_len(22)) {
              matrixER[[iType]] <- cbind(matrixER[[iType]], .computeER4plot(x["Chromosomes", j][iType]))
            }
            colnames(matrixER[[iType]]) <- c("Genome", paste0("Chrom", seq_len(22)))
          }, {
            for (j in what) {
              matrixER[[iType]] <- cbind(matrixER[[iType]], .computeER4plot(x["Chromosomes", j][iType]))
            }
            colnames(matrixER[[iType]]) <- paste0("Chrom", what)
          }
        )
      } else {
        for (j in what) {
          matrixER[[iType]] <- cbind(matrixER[[iType]], .computeER4plot(x["Chromosomes", j][iType]))
        }
        colnames(matrixER[[iType]]) <- paste0("Chrom", what)
      }
    }
  }
  if (ggplot) {
    is.installed <- function(mypkg) is.element(mypkg, utils::installed.packages()[, 1])
    if (!require("ggplot2") | !require("grid")) {
      stop('[Enrichment:plot] "ggPlot2" and "grid" packages must be installed with "ggplot=TRUE".', call. = FALSE)
    }
    multiplot <- function(..., plotlist = NULL, cols = 1, rows = 1, layout = NULL) {
      plots <- c(list(...), plotlist)
      numPlots <- length(plots)
      if (is.null(layout)) {
        layout <- matrix(seq_len(cols * ceiling(numPlots / cols)), ncol = cols, nrow = rows, byrow = TRUE)
      }
      if (numPlots == 1) {
        print(plots[[1]])
      } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:numPlots) {
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
        }
      }
    }
    .ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
      if ((diff(h) %% 360) < 1) {
        h[2] <- h[2] - 360 / n
      }
      grDevices::hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
    }
    listPlots <- list()
    for (iType in type) {
      if (pvalue) {
        if (x@Call$reSample$empiricPvalue) {
          ylab <- "P-Value (Empirical)"
        } else {
          matrixER[[iType]] <- apply(matrixER[[iType]], 2, stats::pnorm, lower.tail = FALSE)
          ylab <- "P-Value (From Z-statistic)"
        }
      } else {
        ylab <- "Z statistic"
      }
      if (ncol(matrixER[[iType]]) > 1) {
        matrixER[[iType]] <- apply(matrixER[[iType]], 2, scale)
      }
      tmpDF <- as.data.frame(t(matrixER[[iType]]))
      cnames <- colnames(tmpDF)
      colnames(tmpDF) <- paste0("R", colnames(tmpDF))
      tmpDF$IID <- factor(colnames(matrixER[[iType]]), levels = c("Genome", paste0("Chrom", seq_len(22))), labels = c("Genome", paste0("Chrom", seq_len(22))))
      tmp <- stats::reshape(tmpDF, idvar = "IID", direction = "long", varying = list(grep("R", colnames(tmpDF))), times = cnames)
      colnames(tmp) <- c("IID", "Resampling", "Z")
      tmp[, "Resampling"] <- as.numeric(tmp[, "Resampling"])

      p <- ggplot2::ggplot(data = tmp, ggplot2::aes_string(x = "Resampling", y = "Z", colour = "IID")) + ggplot2::geom_line()
      noGridColour <- "transparent"
      base_size <- 12
      base_family <- ""
      p <- p + ggplot2::theme(
        line = ggplot2::element_line(colour = "black", size = 0.5, linetype = 1, lineend = "butt"),
        rect = ggplot2::element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
        text = ggplot2::element_text(family = base_family, face = "plain", colour = "black", size = base_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9),
        axis.text = ggplot2::element_text(size = ggplot2::rel(0.8), colour = "black"),
        strip.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(vjust = 1),
        axis.text.y = ggplot2::element_text(hjust = 1),
        axis.ticks = ggplot2::element_line(colour = "black"),
        axis.title.x = ggplot2::element_text(),
        axis.title.y = ggplot2::element_text(angle = 90),
        axis.ticks.length = ggplot2::unit(0.15, "cm"),
        axis.ticks.margin = ggplot2::unit(0.1, "cm"),
        legend.background = ggplot2::element_rect(fill = "white", colour = "black"),
        legend.margin = ggplot2::unit(0.2, "cm"),
        legend.key = ggplot2::element_rect(fill = "white", colour = "black"),
        legend.key.size = ggplot2::unit(1.2, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
        legend.text.align = NULL,
        legend.title = ggplot2::element_text(size = ggplot2::rel(0.8), face = "bold", hjust = 0),
        legend.title.align = NULL,
        legend.position = "right",
        legend.direction = NULL,
        legend.justification = "center",
        legend.box = NULL,
        panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(colour = noGridColour[1]),
        panel.grid.minor = ggplot2::element_line(colour = noGridColour[length(noGridColour)], size = 0.25),
        panel.margin = ggplot2::unit(0.25, "lines"),
        strip.background = ggplot2::element_rect(fill = "black", colour = "black"),
        strip.text.x = ggplot2::element_text(colour = "white"),
        strip.text.y = ggplot2::element_text(angle = -90, colour = "white"),
        plot.background = ggplot2::element_rect(colour = "white"),
        plot.title = ggplot2::element_text(size = ggplot2::rel(1.2)),
        plot.margin = ggplot2::unit(c(1, 1, 0.5, 0.5), "lines"),
        complete = TRUE
      )
      p <- p + ggplot2::xlab(paste(iType, "Resampling"))
      if (ncol(matrixER[[iType]]) > 1) {
        p <- p + ggplot2::ylab(paste(ylab, "(scale and center)"))
      } else {
        p <- p + ggplot2::ylab(ylab)
      }
      p <- p + ggplot2::theme(legend.title = ggplot2::element_blank())
      if (length(what) >= 5 | what == "All") {
        p <- p + ggplot2::theme(legend.background = ggplot2::element_rect(fill = "gray90", linetype = "dotted"))
      } else {
        quarter <- floor(3 / 4 * nrow(tmp)):nrow(tmp)
        rangeZtmp <- range(tmp[, "Z"])
        if (rangeZtmp[1] == rangeZtmp[2]) {
          if (all(rangeZtmp == 0)) {
            rangeZtmp[1] <- -1
            rangeZtmp[2] <- 1
          } else {
            rangeZtmp[1] <- rangeZtmp[1] * 0.90
            rangeZtmp[2] <- rangeZtmp[2] * 1.10
          }
        }
        rangeZ <- seq(rangeZtmp[1], rangeZtmp[2], by = diff(rangeZtmp) * 1 / 3)
        names(rangeZ) <- c("0%", "33%", "66%", "100%")
        rangeQuarter <- range(tmp[quarter, "Z"])
        inf <- apply(sapply(rangeQuarter, function(lim) lim < rangeZ), 1, all)
        sup <- apply(sapply(rangeQuarter, function(lim) lim > rangeZ), 1, all)
        if (sum(inf) <= sum(sup)) {
          p <- p + ggplot2::theme(legend.justification = c(1, 0), legend.position = c(1, 0))
        } else {
          p <- p + ggplot2::theme(legend.justification = c(1, 1), legend.position = c(1, 1))
        }
      }

      if ("Genome" %in% unique(tmp$IID)) {
        p <- p + ggplot2::scale_colour_manual(values = c("black", .ggplotColours(ifelse(length(unique(tmp$IID)) - 1 > 0, length(unique(tmp$IID)) - 1, 1))))
      }
      listPlots[[iType]] <- p
    }
    multiplot(plotlist = listPlots, cols = length(listPlots))
    invisible(listPlots)
  } else {
    graphics::par(mfrow = c(1, length(type)))
    for (iType in type) {
      if (pvalue) {
        if (x@Call$reSample$empiricPvalue) {
          matrixER[[iType]] <- .computeEmpP4plot(x[iType])
          ylab <- "P-Value (Empirical)"
        } else {
          matrixER[[iType]] <- apply(matrixER[[iType]], 2, stats::pnorm, lower.tail = FALSE)
          ylab <- "P-Value (From Z-statistic)"
        }
      } else {
        ylab <- "Z statistic"
      }
      nbCol <- ncol(matrixER[[iType]])
      ylim <- range(stats::na.exclude(matrixER[[iType]]))
      xNames <- rownames(matrixER[[iType]])
      colors <- grDevices::rainbow(nbCol)
      if (nbCol > 1) {
        matrixER[[iType]] <- apply(matrixER[[iType]], 2, scale)
        ylab <- paste(ylab, "(scale and center)")
        graphics::plot(x = xNames, y = matrixER[[iType]][, 1], ylab = ylab, xlab = iType, type = "l", ylim = ylim)
        res <- sapply(seq_len(ncol(matrixER[[iType]][, -1])), function(iER) {
          graphics::lines(x = xNames, y = matrixER[[iType]][, iER + 1], iType = "l", ylim = ylim, col = colors[iER + 1])
        })
      } else {
        graphics::plot(x = xNames, y = matrixER[[iType]][, 1], ylab = ylab, xlab = iType, type = "l", ylim = ylim)
      }
    }
    invisible()
  }
})


methods::setMethod(f = "getEnrichSNP", signature = "Enrichment", definition = function(object, type = "eSNP") {
  alpha <- object["Call"][["readEnrichment"]][["sigThresh"]]
  switch(type,
    "eSNP" = object["Data"][object["Data"][, "PVALUE"] < alpha & object["Data"][, type] == 1, ],
    "xSNP" = {
      if (object["Call"][["readEnrichment"]][["LD"]]) {
        message("Loading ...")
        dataSNP <- object["Data"]
        dataLD <- object["LD"]
        xSNP <- dataSNP[dataSNP[, "PVALUE"] < alpha & dataSNP[, type] == 1, ]
        dataLDtmp <- dataLD[dataLD %in% xSNP[, "SNP"]]
        dataLDtmp <- cbind(SNP_A = names(dataLDtmp), SNP_B = dataLDtmp)
        dataLDtmp <- dataLDtmp[dataLDtmp[, "SNP_A"] %in% dataSNP[dataSNP[, "eSNP"] %in% 1, "SNP"], ]
        xSNPld <- do.call("rbind", by(dataLDtmp, dataLDtmp[, "SNP_B"], function(iDta) {
          cbind(xSNP = unique(as.character(iDta[, "SNP_B"])), LD_with_eSNP = paste(iDta[, "SNP_A"], collapse = ";"))
        }))
        merge(xSNP, xSNPld, by.x = "SNP", by.y = "xSNP")
      } else {
        warning('[Enrichment:getEnrichSNP] significant "eSNP" are returned instead of "xSNP",\n    "readEnrichment" should be run with "LD=TRUE".', call. = FALSE)
        object["Data"][object["Data"][, "PVALUE"] < alpha & object["Data"][, type] == 1, ]
      }
    },
    stop('[Enrichment:getEnrichSNP] "type" should be equal to "eSNP" or "xSNP".', call. = FALSE)
  )
})
