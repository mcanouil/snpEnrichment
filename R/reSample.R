#' Compute enrichment analysis on an \code{\linkS4class{Enrichment}} object
#'
#' After \code{\link{initFiles}} and \code{\link{readEnrichment}} has been run.
#' \code{\link{reSample}} computes a statistic value and a p-value for each
#' chromosomes and for the whole genome.
#'
#'
#' @name reSample
#' @aliases reSample reSample-methods reSample,Chromosome-method
#' reSample,Enrichment-method reSample,ANY-method
#' @docType methods
#' @param object [Enrichment or Chromosome]: an object to be updated. It is
#' intended, an object returned by the \code{\link{readEnrichment}} function.
#' @param nSample [numeric]: the number of resampling done by
#' \code{\link{reSample}} for p-values computation (minimum is 100).
#' @param empiricPvalue [logical]: \code{empiricPvalue=TRUE} (default) compute
#' PValue based on the null distribution (resampling).  If
#' \code{empiricPvalue=TRUE}, the empirical p-values are computed instead.
#' @param sigThresh [numeric]: statistical threshold for signal (e.g.
#' \code{sigThresh = 0.05} for a given GWAS signal) used to compute an
#' Enrichment Ratio.
#' @param MAFpool [vector(numeric)]: either a numeric vector giving the breaks
#' points of intervals into which SNP's MAF (Minor Allele Frequency) is to be
#' split.
#' @param mc.cores [numeric]: the number of cores to use (default is
#' \code{mc.cores=1}), i.e. at most how many child processes will be run
#' simultaneously.  Must be at least one, and parallelization requires at least
#' two cores.
#' @param onlyGenome [logical]: \code{onlyGenome=TRUE} (default) compute
#' resampling step for all chromosomes.
#' @return Return the object given in argument, updated by the resampling
#' results.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords reSample methods
#' @examples
#'
#' \dontrun{data(toyEnrichment)
#' reSample(object = toyEnrichment,
#'          nSample = 10,
#'          empiricPvalue = TRUE,
#'          MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
#'          onlyGenome = TRUE)
#' toyEnrichment}
#'
methods::setGeneric(
  name = "reSample",
  def = function(object, nSample = 100, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = TRUE, ...) standardGeneric("reSample")
)

methods::setMethod(f = "reSample", signature = "ANY", definition = function(object, nSample, sigThresh, MAFpool, mc.cores) {
  if (!(is.enrichment(object) & is.chromosome(object))) {
    stop('[Method:reSample] not available for "', class(object), '" object.', call. = FALSE)
  }
})

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
