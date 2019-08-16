#' Class [EnrichSNP-class]
#'
#' @slot List A list of SNPs used to compute enrichment (*e.g.*, eSNP or xSNP).
#' @slot Table A contingency table with SNPs (columns) and P-Values from signal (rows).
#' @slot EnrichmentRatio Enrichment Ratio is computed on the contingency table (`Table` slot).
#' @slot Z A statistic computed from `EnrichmentRatio` and resampling results.
#' @slot PValue P-Value associated with the statistic `Z`.
#' @slot Resampling A matrix with by row, the contingency table and the odds ratio for each resampling.
#'
#' @name EnrichSNP-class
#' @exportClass EnrichSNP
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

#' Constructor generic for [EnrichSNP-class]
#'
#' @name EnrichSNP
#' @rdname EnrichSNP-class
#' @keywords internal
methods::setGeneric(name = "enrichSNP", def = function(List, Table, EnrichmentRatio, Z, PValue, Resampling) standardGeneric("enrichSNP"))

#' Constructor method for [EnrichSNP-class]
#'
#' @name enrichSNP
#' @rdname EnrichSNP-class
#' @aliases enrichSNP,ANY-method
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

#' Getter method for [EnrichSNP-class]
#'
#' @name [
#' @rdname internal
#' @aliases [,EnrichSNP-method
#' @exportMethod [
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

#' Setter method for [EnrichSNP-class]
#'
#' @name [<-
#' @rdname internal
#' @aliases [<-,EnrichSNP-method
#' @exportMethod [<-
methods::setMethod(f = '[<-', signature = "EnrichSNP", definition = function(x, i, j, value) {
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

#' Print method for [EnrichSNP-class]
#'
#' @name print
#' @rdname internal
#' @aliases print,EnrichSNP-method
#' @exportMethod print
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

#' Show method for [EnrichSNP-class]
#'
#' @name show
#' @rdname internal
#' @aliases show,EnrichSNP-method
#' @exportMethod show
methods::setMethod(f = "show", signature = "EnrichSNP", definition = function(object) {
  cat("     ~~~ Class:", class(object), "~~~\n")
  .EnrichSNP.show(object)
  invisible(object)
})

#' Class [Chromosome-class]
#'
#' @slot Data A [data.frame] with 6 columns (`"SNP"`, `"PVALUE"`, `"CHR"`, `"MAF"`, `"eSNP"`, `"xSNP"`).
#'   Where `"eSNP"` and `"xSNP"` are logical columns defining the lists of SNPs (extended or not).
#' @slot LD A data.frame which contains LD informations between SNPs (computed with [writeLD] or PLINK).
#' @slot eSNP Contain a [EnrichSNP-class] object for a list of SNPs (eSNP).
#' @slot xSNP Contain a [EnrichSNP-class] object for a extended list of SNPs (xSNP).
#'
#' @name Chromosome-class
#' @exportClass Chromosome
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

#' Constructor generic for [Chromosome-class]
#'
#' @name chromosome
#' @rdname Chromosome-class
#' @aliases chromosome
#' @exportMethod chromosome
methods::setGeneric(name = "chromosome", def = function(Data, LD, eSNP, xSNP) standardGeneric("chromosome"))

#' Constructor method for [Chromosome-class]
#'
#' @name chromosome
#' @rdname Chromosome-class
#' @aliases chromosome,ANY-method
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

#' Getter method for [Chromosome-class]
#'
#' @name [
#' @rdname internal
#' @aliases [,Chromosome-method
#' @exportMethod [
methods::setMethod(f = "[", signature = "Chromosome", definition = function(x, i, j, drop) {
  switch(EXPR = i,
    "Data" = x@Data,
    "LD" = x@LD,
    "eSNP" = x@eSNP,
    "xSNP" = x@xSNP,
    stop("[Chromosome:get] ", i, ' is not a "Chromosome" slot.', call. = FALSE)
  )
})

#' Setter method for [Chromosome-class]
#'
#' @name [<-
#' @rdname internal
#' @aliases [<-,Chromosome-method
#' @exportMethod [<-
methods::setMethod(f = '[<-', signature = "Chromosome", definition = function(x, i, j, value) {
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

#' Print method for [Chromosome-class]
#'
#' @name print
#' @rdname internal
#' @aliases print,Chromosome-method
#' @exportMethod print
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

#' Show method for [Chromosome-class]
#'
#' @name show
#' @rdname internal
#' @aliases show,Chromosome-method
#' @exportMethod show
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
#'   [Enrichment-class] is the object containing the results for all [Chromosome-class] object
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
#' @slot eSNP Contain a [EnrichSNP-class] object for a list of SNPs (eSNP).
#' @slot xSNP Contain a [EnrichSNP-class] object for a extended list of SNPs (xSNP).
#' @slot Chromosomes A list of 22 [Chromosome-class] objects.
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
#' @aliases enrichment,ANY-method
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
#' @name [
#' @rdname internal
#' @aliases [,Enrichment-method
#' @exportMethod [
methods::setMethod(f = "[", signature = "Enrichment", definition = function(x, i, j, drop) {
  nbChr <- length(x@Chromosomes)
  if (missing(j)) {
    switch(EXPR = i,
      "Loss" = x@Loss,
      "Data" = {
        resData <- .mclapply(
          X = seq_len(22),
          mc.cores = min(22, parallel::detectCores()),
          FUN = function(iChr) x@Chromosomes[[iChr]]@Data
        )
        do.call("rbind", resData)
      },
      "LD" = {
        resLD <- .mclapply(
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
            resData <- .mclapply(
              X = seq_len(22),
              mc.cores = min(22, parallel::detectCores()),
              FUN = function(iChr) x@Chromosomes[[iChr]]@Data
            )
            do.call("rbind", resData)
          },
          "LD" = {
            resLD <- .mclapply(
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
            resData <- .mclapply(
              X = j,
              mc.cores = min(length(j), parallel::detectCores()),
              FUN = function(iChr) x@Chromosomes[[iChr]]@Data
            )
            do.call("rbind", resData)
          },
          "LD" = {
            resLD <- .mclapply(
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
#' @name [<-
#' @rdname internal
#' @aliases [<-,Enrichment-method
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
#' @name print
#' @rdname internal
#' @aliases print,Enrichment-method
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
#' @name show
#' @rdname internal
#' @aliases show,Enrichment-method
#' @exportMethod show
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
#' @aliases computeER,Chromosome-method
#' @keywords internal
methods::setMethod(f = "computeER", signature = "Chromosome", definition = function(object, sigThresh = 0.05, mc.cores = 1) {
  if (missing(object)) {
    stop('[Chromosome:computeER] "Chromosome" object is required.', call. = FALSE)
  }
  object@Chromosomes <- .mclapply(object@Chromosomes, mc.cores = mc.cores, function(chr) {
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
#' @aliases computeER,Enrichment-method
#' @keywords internal
methods::setMethod(f = "computeER", signature = "Enrichment", definition = function(object, sigThresh = 0.05, mc.cores = 1) {
  if (missing(object)) {
    stop('[Enrichment:computeER] "Enrichment" object is required.', call. = FALSE)
  }
  rowNames <- c(paste0("P>=", sigThresh), paste0("P<", sigThresh))
  object@Chromosomes <- .mclapply(object@Chromosomes, mc.cores = mc.cores, function(chr) {
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
#' @aliases doLDblock,Chromosome-method
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
  byBlock <- unique(.mclapply(byBlock, mc.cores = nbCores, function(i) {
    names(i) <- NULL
    i
  }))
  LDblockTmp <- .mclapply(seq_along(byBlock), mc.cores = nbCores, function(jBlock) {
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
  resParallel <- .mclapply(seq_len(nrow(blockLim)), mc.cores = nbCores, function(li) {
    blockLim[li, "LENGTH"] <- as.integer(blockLim[li, "MAX"]) - as.integer(blockLim[li, "MIN"])
    blockLim[li, ]
  })
  blockLim <- do.call("rbind", resParallel)
  rm(resParallel)

  data <- data[order(data[, "POS"]), ]
  data[, c("MIN", "MAX", "IDBLOCK", "LENGTH", "MAFmedian")] <- as.numeric(NA)
  tmpChr <- .mclapply(seq_len(nrow(blockLim)), mc.cores = nbCores, function(i) {
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
