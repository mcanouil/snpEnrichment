#' Is an Chromosome object
#'
#' 'is.chromosome' returns 'TRUE' if 'x' is an \code{\linkS4class{Chromosome}}
#' object and 'FALSE' otherwise.
#'
#'
#' @aliases is.chromosome is.chromosome-methods is.chromosome,ANY-method
#' @param object [ANY]: object to be tested.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords chromosome is is.chromosome
#' @examples
#'
#' a <- chromosome()
#' c <- chromosome()
#' is.chromosome(list())                # FALSE
#' is.chromosome(1)                     # FALSE
#' is.chromosome(a)                     # TRUE
#' is.chromosome(c(a, c))               # TRUE TRUE
#' is.chromosome(list(a, b = "char"))   # TRUE FALSE
#' is.chromosome(c(a, b = list(12, c))) # TRUE FALSE TRUE
#'
methods::setGeneric(
  name = "is.chromosome",
  def = function(object) standardGeneric("is.chromosome")
)
methods::setMethod(f = "is.chromosome", signature = "ANY", definition = function(object) {
  if (length(object) > 1) {
    sapply(object, is.chromosome)
  } else {
    class(object) == "Chromosome"
  }
})
