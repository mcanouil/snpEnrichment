#' Class \code{\linkS4class{Chromosome}}
#'
#' This class is defined to summarize the enrichment analysis about a
#' chromosome.
#'
#'
#' @name Chromosome-class
#' @aliases Chromosome-class Chromosome [,Chromosome-method
#' [,Chromosome,ANY,ANY,ANY-method [<-,Chromosome-method
#' [<-,Chromosome,ANY,ANY,ANY-method show,Chromosome-method chromosome
#' chromosome-methods chromosome,ANY-method
#' @docType class
#' @note \code{\linkS4class{Chromosome}} object is not intended to be used
#' alone on this version.\ It is a part of the \code{\linkS4class{Enrichment}}
#' object.
#' @section Objects from the Class: \code{\link{chromosome}} is defined to
#' build an object of class \code{\linkS4class{Chromosome}} in order to compute
#' an enrichment analysis.  A \code{\linkS4class{Chromosome}} object contains
#' the original data, a list of SNPs, some results and resampling data.\cr
#'
#' When created, an \code{\linkS4class{Chromosome}} object is "empty".
#' \code{\link{readEnrichment}} initializes a \code{\linkS4class{Chromosome}}
#' object with value from PLINK computation and user's files.  In this step,
#' only the fields "Data", "LD", "SNP" are filled.  \code{reSample} fills the
#' fields: Table, EnrichmentRatio, Z, PValue and Resampling of a
#' \code{\linkS4class{Chromosome}}.
#'
#' Note that if \code{\link{reSample}} is executed on an
#' \code{\linkS4class{Chromosome}} every new resampling is added to the
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
#' @keywords classes class chromosome
#' @examples
#'
#' Data <- structure(
#'     list(
#'         SNP = c("rs4970420", "rs3766178",
#'                 "rs10910030", "rs10910036",
#'                 "rs2234167", "rs6661861"),
#'         PVALUE = c(0.9244, 0.167, 0.01177, 0.4267, 0.9728, 0.4063),
#'         CHR = c(1, 1, 1, 1, 1, 1),
#'         POS = c(1106473, 1478180, 2035684, 2183754, 2494330, 3043121),
#'         MAF = c(0.2149, 0.3117, 0.374, 0.3753, 0.1565, 0.06101),
#'         eSNP = c(0, 1, 1, 0, 0, 0),
#'         xSNP = c(0, 1, 1, 0, 0, 0)
#'     ),
#'     .Names = c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP", "xSNP"),
#'     row.names = c("rs4970420", "rs3766178",
#'                   "rs10910030", "rs10910036",
#'                   "rs2234167", "rs6661861"),
#' class = "data.frame")
#'
#' toyChr <- chromosome(Data = Data)
#' show(toyChr)
#' toyChr
#'
#' toyChr <- chromosome()
#' toyChr["Data"] <- Data
#' toyChr
#'
#' is.chromosome(toyChr) # TRUE
#'
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

methods::setGeneric(
  name = "chromosome",
  def = function(Data, LD, eSNP, xSNP) standardGeneric("chromosome")
)

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


methods::setMethod(f = "[", signature = "Chromosome", definition = function(x, i, j, drop) {
  switch(EXPR = i,
    "Data" = x@Data,
    "LD" = x@LD,
    "eSNP" = x@eSNP,
    "xSNP" = x@xSNP,
    stop("[Chromosome:get] ", i, ' is not a "Chromosome" slot.', call. = FALSE)
  )
})


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
