#' Class [EnrichSNP]
#'
#' @slot List A list of SNPs used to compute enrichment (*e.g.*, eSNP or xSNP).
#' @slot Table A contingency table with SNPs (columns) and P-Values from signal (rows).
#' @slot EnrichmentRatio Enrichment Ratio is computed on the contingency table (`Table` slot).
#' @slot Z A statistic computed from `EnrichmentRatio` and resampling results.
#' @slot PValue P-Value associated with the statistic `Z`.
#' @slot Resampling A matrix with by row, the contingency table and the odds ratio for each resampling.
#'
#' @name class-EnrichSNP
#' @rdname internal
#' @keywords internal
methods::setClass(
  Class = "EnrichSNP",
  representation = methods::representation(
    List = "character",
    Table = "matrix",
    EnrichmentRatio = "numeric",
    Z = "numeric",
    PValue = "numeric",
    Resampling = "matrix"
  ),
  prototype = methods::prototype(
    List = character(),
    Table = matrix(0, ncol = 2, nrow = 2),
    EnrichmentRatio = numeric(),
    Z = numeric(),
    PValue = numeric(),
    Resampling = matrix(0, ncol = 5, nrow = 0)
  )
)

#' Constructor generic for [EnrichSNP]
#'
#' @name EnrichSNP
#' @rdname internal
#' @keywords internal
methods::setGeneric(name = "enrichSNP", def = function(List, Table, EnrichmentRatio, Z, PValue, Resampling) standardGeneric("enrichSNP"))

#' Constructor method for [EnrichSNP]
#'
#' @name EnrichSNP
#' @rdname internal
#' @keywords internal
methods::setMethod(f = "enrichSNP", signature = "ANY", definition = function(List, Table, EnrichmentRatio, Z, PValue, Resampling) {
  if (missing(List)) List <- character()
  if (missing(Table)) Table <- matrix(0, ncol = 2, nrow = 2)
  if (missing(EnrichmentRatio)) EnrichmentRatio <- numeric()
  if (missing(Z)) Z <- numeric()
  if (missing(PValue)) PValue <- numeric()
  if (missing(Resampling)) Resampling <- matrix(0, ncol = 5, nrow = 0)
  methods::new(
    "EnrichSNP",
    List = List,
    Table = Table,
    EnrichmentRatio = EnrichmentRatio,
    Z = Z,
    PValue = PValue,
    Resampling = Resampling
  )
})

#' Getter method for [EnrichSNP]
#'
#' @name EnrichSNP
#' @rdname internal
#' @keywords internal
methods::setMethod(f = "[", signature = "EnrichSNP", definition = function(x, i, j, drop) {
  switch(EXPR = i,
    "List" = x@List,
    "Table" = x@Table,
    "EnrichmentRatio" = x@EnrichmentRatio,
    "Z" = x@Z,
    "PValue" = x@PValue,
    "Resampling" = x@Resampling,
    stop("[EnrichSNP:get] ", i, ' is not a "EnrichSNP" slot.', call. = FALSE)
  )
})

#' Setter method for [EnrichSNP]
#'
#' @name EnrichSNP
#' @rdname internal
#' @keywords internal
methods::setMethod(f = "[<-", signature = "EnrichSNP", definition = function(x, i, j, value) {
  switch(EXPR = i,
    "List" = x@List <- value,
    "Table" = x@Table <- value,
    "EnrichmentRatio" = x@EnrichmentRatio <- value,
    "Z" = x@Z <- value,
    "PValue" = x@PValue <- value,
    "Resampling" = x@Resampling <- value,
    stop("[EnrichSNP:set] ", i, ' is not a "EnrichSNP" slot.', call. = FALSE)
  )
  methods::validObject(x)
  x
})

#' Print method for [EnrichSNP]
#'
#' @name EnrichSNP
#' @rdname internal
#' @keywords internal
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

#' Show method for [EnrichSNP]
#'
#' @name EnrichSNP
#' @rdname internal
#' @keywords internal
methods::setMethod(f = "show", signature = "EnrichSNP", definition = function(object) {
  cat("     ~~~ Class:", class(object), "~~~\n")
  .EnrichSNP.show(object)
  invisible(object)
})

#' Class [Chromosome]
#'
#' @slot Data A [data.frame] with 6 columns (`"SNP"`, `"PVALUE"`, `"CHR"`, `"MAF"`, `"eSNP"`, `"xSNP"`).
#'   Where `"eSNP"` and `"xSNP"` are logical columns defining the lists of SNPs (extended or not).
#' @slot LD A data.frame which contains LD informations between SNPs (computed with [writeLD] or PLINK).
#' @slot eSNP Contain a [EnrichSNP] object for a list of SNPs (eSNP).
#' @slot xSNP Contain a [EnrichSNP] object for a extended list of SNPs (xSNP).
#'
#' @name class-Chromosome
#' @rdname internal
#' @keywords internal
methods::setClass(
  Class = "Chromosome",
  representation = methods::representation(
    Data = "data.frame",
    LD = "character",
    eSNP = "EnrichSNP",
    xSNP = "EnrichSNP"
  ),
  prototype = methods::prototype(
    Data = data.frame(),
    LD = character(),
    eSNP = enrichSNP(),
    xSNP = enrichSNP()
  )
)

#' Constructor generic for [Chromosome]
#'
#' @name Chromosome
#' @rdname internal
#' @keywords internal
methods::setGeneric(name = "chromosome", def = function(Data, LD, eSNP, xSNP) standardGeneric("chromosome"))

#' Constructor method for [Chromosome]
#'
#' @name Chromosome
#' @rdname internal
#' @keywords internal
methods::setMethod(f = "chromosome", signature = "ANY", definition = function(Data, LD, eSNP, xSNP) {
  if (missing(Data)) {
    Data <- data.frame()
    if (missing(eSNP)) {
      eSNP <- enrichSNP()
      xSNP <- enrichSNP()
    }
  } else {
    if (missing(eSNP)) {
      eSNP <- enrichSNP(List = Data[Data[, "eSNP"] == 1, "SNP"])
      xSNP <- enrichSNP(List = Data[Data[, "xSNP"] == 1, "SNP"])
    }
  }
  if (missing(LD)) LD <- character()
  methods::new("Chromosome", Data = Data, LD = LD, eSNP = eSNP, xSNP = xSNP)
})

#' Getter method for [Chromosome]
#'
#' @name Chromosome
#' @rdname internal
#' @keywords internal
methods::setMethod(f = "[", signature = "Chromosome", definition = function(x, i, j, drop) {
  switch(EXPR = i,
    "Data" = x@Data,
    "LD" = x@LD,
    "eSNP" = x@eSNP,
    "xSNP" = x@xSNP,
    stop("[Chromosome:get] ", i, ' is not a "Chromosome" slot.', call. = FALSE)
  )
})

#' Setter method for [Chromosome]
#'
#' @name Chromosome
#' @rdname internal
#' @keywords internal
methods::setMethod(f = "[<-", signature = "Chromosome", definition = function(x, i, j, value) {
  switch(EXPR = i,
    "Data" = x@Data <- value,
    "LD" = x@LD <- value,
    "eSNP" = x@eSNP <- value,
    "xSNP" = x@xSNP <- value,
    stop("[Chromosome:get] ", i, ' is not a "Chromosome" slot.', call. = FALSE)
  )
  methods::validObject(x)
  invisible(x)
})

#' Print method for [Chromosome]
#'
#' @name Chromosome
#' @rdname internal
#' @keywords internal
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

#' Show method for [Chromosome]
#'
#' @name Chromosome
#' @rdname internal
#' @keywords internal
methods::setMethod(f = "show", signature = "Chromosome", definition = function(object) {
  cat("     ~~~ Class:", class(object), "~~~\n")
  .Chromosome.show(object)
  invisible(object)
})

#' An S4 class [Enrichment-class]
#'
#' This class is defined to summarize the enrichment analysis on each chromosomes and the whole genome.
#'
#' @section Objects from the Class: [Enrichment-class] is defined to build an object
#'   of class [Enrichment-class] in order to compute an enrichment analysis.
#'   [Enrichment-class] is the object containing the results for all [Chromosome] object
#'   and for the whole genome.
#'
#'   When an [Enrichment-class] object is created, it contains a list of SNPs (*e.g.*, eSNPs).
#'   All the others slots are "empty".
#'   After [reSample] is ran on an [Enrichment-class] object, the slots:
#'   Table, EnrichmentRatio, Z, PValue and Resampling are filled.
#'
#'   Note that if [reSample] is executed on an [Enrichment-class] every new resampling is added to the
#'   original ones, pre-existing statistics are erased and computed again with the new resampling set.
#'
#' @slot Signal A three columns `data.frame`: "SNP", "PVALUE" and "IN" (*e.g.*, GWAS).
#'   "IN" is computed during the reading step and gives informations about which SNPs are kept
#'   for the enrichment analysis.
#' @slot Loss A four columns data.frame: "Rows", "Unique", "Intersect.Ref.Signal" and "CIS".
#'   This slot gives information on data losses.
#' @slot Call Each parameters used for the reading or resampling step are stored in this slot.
#' @slot eSNP Contain a [EnrichSNP] object for a list of SNPs (eSNP).
#' @slot xSNP Contain a [EnrichSNP] object for a extended list of SNPs (xSNP).
#' @slot Chromosomes A list of 22 [Chromosome] objects.
#'
#' @examples
#'
#' data(toyEnrichment)
#' toyEnrich <- enrichment()
#' show(toyEnrich)
#'
#' toyEnrich["Loss"] <- toyEnrichment["Loss"]
#' toyEnrich["Loss"]
#'
#' toyEnrich <- enrichment(Loss = toyEnrichment["Loss"], eSNP = toyEnrichment["eSNP"])
#' toyEnrich <- enrichment(Loss = toyEnrichment["Loss"])
#'
#' if (interactive()) {
#'   reSample(
#'     object = toyEnrichment,
#'     nSample = 10,
#'     empiricPvalue = TRUE,
#'     MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
#'     mc.cores = 1,
#'     onlyGenome = TRUE
#'   )
#'   print(toyEnrichment)
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
#'
#'   toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)
#'   print(toyEnrichment_exclude)
#' }
#'
#' @name Enrichment-class
#' @exportClass Enrichment
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

#' Constructor generic for [Enrichment-class]
#'
#' @name enrichment
#' @rdname Enrichment-class
#' @exportMethod enrichment
methods::setGeneric(
  name = "enrichment",
  def = function(Loss, Call, eSNP, xSNP, Chromosomes) standardGeneric("enrichment")
)

#' Constructor method for [Enrichment-class]
#'
#' @name enrichment
#' @rdname Enrichment-class
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

#' Getter method for [Enrichment-class]
#'
#' @name enrichment
#' @rdname Enrichment-class
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

#' Setter method for [Enrichment-class]
#'
#' @name enrichment
#' @rdname Enrichment-class
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

#' Print method for [Enrichment-class]
#'
#' @name enrichment
#' @rdname Enrichment-class
#' @exportMethod print
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

#' Show method for [Enrichment-class]
#'
#' @name enrichment
#' @rdname Enrichment-class
methods::setMethod(f = "show", signature = "Enrichment", definition = function(object) {
  cat("    ~~~ Class:", class(object), "~~~ \n")
  .Enrichment.show(object)
  invisible(object)
})

#' @name computeER
#' @rdname internal
#' @keywords internal
methods::setGeneric(name = "computeER", def = function(object, sigThresh = 0.05, mc.cores = 1) standardGeneric("computeER"))
#' @name computeER-Chromosome
#' @rdname internal
#' @keywords internal
methods::setMethod(f = "computeER", signature = "Chromosome", definition = function(object, sigThresh = 0.05, mc.cores = 1) {
  if (missing(object)) {
    stop('[Chromosome:computeER] "Chromosome" object is required.', call. = FALSE)
  }
  object@Chromosomes <- mclapply2(object@Chromosomes, mc.cores = mc.cores, function(chr) {
    data <- chr@Data
    chrLD <- length(chr@LD)
    for (type in c("eSNP", "xSNP")) {
      if (!(chrLD == 0 & type == "xSNP")) {
        snpEnrich <- table(factor(chr@Data[, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE)), factor(chr@Data[, type], levels = c(0, 1)))
        colnames(snpEnrich) <- c("otherSNP", type)
        rownames(snpEnrich) <- eval(parse(text = paste0('c("P>=', sigThresh, '", "P<', sigThresh, '")')))
        chr[type] <- enrichSNP(List = chr[type]@List, Table = unclass(snpEnrich), EnrichmentRatio = .enrichmentRatio(snpEnrich))
      }
    }
    chr
  })
  object
})
#' @name computeER-Enrichment
#' @rdname internal
#' @keywords internal
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

#' @name doLDblock
#' @rdname internal
#' @keywords internal
methods::setGeneric(name = "doLDblock", def = function(object, mc.cores = 1) standardGeneric("doLDblock"))
#' @name doLDblock-Chromosome
#' @rdname internal
#' @keywords internal
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

#' @name reset
#' @rdname internal
#' @exportMethod reset
methods::setGeneric(name = "reset", def = function(object, i) standardGeneric("reset"))
#' @name reset-ANY
#' @rdname internal
methods::setMethod(f = "reset", signature = "ANY", definition = function(object, i) {
  if (!(is.enrichment(object) & is.chromosome(object))) {
    stop('[Method:reset] not available for "', class(object), '" object.', call. = FALSE)
  }
})
#' @name reset-EnrichSNP
#' @rdname internal
methods::setMethod(f = "reset", signature = "EnrichSNP", definition = function(object, i) {
  switch(EXPR = i,
    "List" = object@List <- character(),
    "Table" = object@Table <- matrix(0, ncol = 2, nrow = 2),
    "EnrichmentRatio" = object@EnrichmentRatio <- numeric(),
    "Z" = object@Z <- numeric(),
    "PValue" = object@PValue <- numeric(),
    "Resampling" = object@Resampling <- matrix(0, ncol = 5, nrow = 0),
    stop("[EnrichSNP:reset] ", i, ' is not a "EnrichSNP" slot.', call. = FALSE)
  )
  object
})
#' @name reset-Chromosome
#' @rdname internal
methods::setMethod(f = "reset", signature = "Chromosome", definition = function(object, i) {
  switch(EXPR = i,
    "Data" = object@Data <- data.frame(),
    "LD" = object@LD <- character(),
    "eSNP" = object@eSNP <- enrichSNP(),
    "xSNP" = object@xSNP <- enrichSNP(),
    "List" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@List <- reset(object[type], "List")
      }
    },
    "Table" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@Table <- matrix(0, ncol = 2, nrow = 2)
      }
    },
    "EnrichmentRatio" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@EnrichmentRatio <- numeric()
      }
    },
    "Z" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@Z <- numeric()
      }
    },
    "PValue" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@PValue <- numeric()
      }
    },
    "Resampling" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@Resampling <- matrix(0, ncol = 5, nrow = 0)
      }
    },
    stop("[Enrichment:reset] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
  )
  object
})
#' @name reset-Enrichment
#' @rdname internal
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

#' Compute enrichment analysis on an [Enrichment-class] object
#'
#' After [initFiles] and [readEnrichment] has been run.
#' [reSample] computes a statistic value and a p-value for each chromosomes and for the whole genome.
#'
#' @param object An object to be updated. It is intended, an object returned by the [readEnrichment] function.
#' @param nSample The number of resampling done by [reSample] for p-values computation (minimum is 100).
#' @param empiricPvalue `TRUE` (default) compute PValue based on the null distribution (resampling).
#'   If `TRUE`, the empirical p-values are computed instead.
#' @param sigThresh Statistical threshold for signal (*e.g.*,`0.05` for a given GWAS signal)
#'   used to compute an Enrichment Ratio.
#' @param MAFpool Either a numeric vector giving the breaks points of intervals into which SNP's MAF
#'   (Minor Allele Frequency) is to be split.
#' @param mc.cores The number of cores to use (default is `1`, *i.e.*, at most how many child processes
#'   will be run simultaneously. Must be at least one, and parallelization requires at least two cores.
#' @param onlyGenome Compute resampling step for all chromosomes (default `TRUE`).
#'
#' @return Return the object given in argument, updated by the resampling results.
#'
#' @examples
#' if (interactive()) {
#'   data(toyEnrichment)
#'   reSample(
#'     object = toyEnrichment,
#'     nSample = 10,
#'     empiricPvalue = TRUE,
#'     MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
#'     onlyGenome = TRUE
#'   )
#'   toyEnrichment
#' }
#'
#' @name reSample
#' @exportMethod reSample
methods::setGeneric(
  name = "reSample",
  def = function(object, nSample = 100, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = TRUE, ...) standardGeneric("reSample")
)
#' @name reSample-Chromosome
#' @rdname reSample
methods::setMethod(
  f = "reSample",
  signature = "Chromosome",
  definition = function(
    object, nSample = 100, empiricPvalue = TRUE, sigThresh = 0.05,
    MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1
  ) {
    if (missing(object)) {
      stop('[Enrichment:reSample] "Enrichment" object is required.', call. = FALSE)
    }
    if (nSample < 10) {
      nSample <- 10
      warning("[Enrichment:reSample] nSample was increased to 10.", call. = FALSE)
    }
    .reSample(
      object = object,
      nSample = nSample,
      empiricPvalue = empiricPvalue,
      sigThresh = sigThresh,
      MAFpool = MAFpool,
      mc.cores = mc.cores
    )
  }
)
#' @name reSample-ANY
#' @rdname reSample
methods::setMethod(f = "reSample", signature = "ANY", definition = function(object, nSample, sigThresh, MAFpool, mc.cores) {
  if (!(is.enrichment(object) & is.chromosome(object))) {
    stop('[Method:reSample] not available for "', class(object), '" object.', call. = FALSE)
  }
})
#' @name reSample-Enrichment
#' @rdname reSample
#' @rdname Enrichment-class
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

#' Get all eSNP/xSNP which are enriched
#'
#' [getEnrichSNP] get all eSNP/xSNP in a [Enrichment-class] object which are significant in the signal
#' according to `sigThresh` defined in [readEnrichment].
#'
#' @param object An object of class [Enrichment-class].
#' @param type Extract `"eSNP"` or `"xSNP"`" data.
#'
#' @return Return a `data.frame` with eSNP/xSNP which are enriched in signal given to `signalFile`
#'   in function [readEnrichment].
#'
#' @examples
#' if (interactive()) {
#'   data(toyEnrichment)
#'   eSNPenriched <- getEnrichSNP(object = toyEnrichment, type = "eSNP")
#'   head(eSNPenriched)
#' }
#'
#' @name getEnrichSNP
#' @exportMethod getEnrichSNP
methods::setGeneric(
  name = "getEnrichSNP",
  def = function(object, type = "eSNP") standardGeneric("getEnrichSNP")
)
#' @name getEnrichSNP-ANY
#' @rdname getEnrichSNP
methods::setMethod(f = "getEnrichSNP", signature = "ANY", definition = function(object, type = "eSNP") {
  if (!(is.enrichment(object))) {
    stop('[Method:getEnrichSNP] not available for "', class(object), '" object.', call. = FALSE)
  }
})
#' @name getEnrichSNP-Enrichment
#' @rdname Enrichment-class
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

#' Plot method (S4) for [Enrichment-class] object
#'
#' [plot] is a generic function for plotting of R objects.
#' The function invokes particular [methods] which depend on the [class] of the first argument.
#'
#' @param x An object of class [Enrichment-class] which the Z statistics or p-values have to be drawn.
#' @param what Default `"Genome"`. Plot Z statistics or p-values for genome only
#'   (what must be: `"All"`, `"Genome"` or numeric vector).
#' @param type Plot the selected analysis for `"eSNP"` and/or `"xSNP"`.
#' @param ggplot Use ggplot (default `FALSE`) instead of classic plot method.
#' @param pvalue If `TRUE`, p-value convergense is plotted. Otherwise, Z statistic is plotted.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see [par]).
#'
#' @examples
#' if (interactive()) {
#'   data(toyEnrichment)
#'   reSample(toyEnrichment, 10)
#'   plot(toyEnrichment)
#' }
#'
#' @name plot
#' @exportMethod plot
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

