#' Class \code{"\linkS4class{EnrichSNP}"}
#'
#' This class is defined to summarize the enrichment analysis. It's a part of
#' \code{\linkS4class{Chromosome}} and \code{\linkS4class{Enrichment}} classes.
#'
#'
#' @name EnrichSNP-class
#' @aliases EnrichSNP-class EnrichSNP [,EnrichSNP-method
#' [,EnrichSNP,ANY,ANY,ANY-method [<-,EnrichSNP-method
#' [<-,EnrichSNP,ANY,ANY,ANY-method show,EnrichSNP-method
#' print,EnrichSNP-method
#' @docType class
#' @note \code{\linkS4class{EnrichSNP}} object is not intended to be use
#' directly by user. It is a part of the \code{\linkS4class{Enrichment}} and
#' \code{\linkS4class{Chromosome}} object.
#' @section Slots: \describe{ \item{List}{[vector(character)]: a list of SNPs
#' used to compute enrichment (e.g. eSNP or xSNP).} \item{Table}{[matrix]:
#' Contingency table with SNPs (columns) and P-Values from signal (rows).}
#' \item{EnrichmentRatio}{[numeric]: Enrichment Ratio is computed on the
#' contingency table (\code{Table} slot).} \item{Z}{[numeric]: A statistic
#' computed from \code{EnrichmentRatio} and resampling results.}
#' \item{PValue}{[numeric]: P-Value associated with the statistic \code{Z}.}
#' \item{Resampling}{[matrix]: A matrix with by row, the contingency table and
#' the odds ratio for each resampling.} }
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords classes class enrichSNP
#'
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

methods::setGeneric(
  name = "enrichSNP",
  def = function(List, Table, EnrichmentRatio, Z, PValue, Resampling) standardGeneric("enrichSNP")
)

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
