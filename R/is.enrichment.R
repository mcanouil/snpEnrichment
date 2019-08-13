#' Is an Enrichment object
#'
#' 'is.enrichment' returns 'TRUE' if 'x' is an \code{\linkS4class{Enrichment}}
#' object and 'FALSE' otherwise.
#'
#'
#' @aliases is.enrichment is.enrichment-methods is.enrichment,ANY-method
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
#' @keywords enrichment is is.enrichment
#' @examples
#'
#' a <- enrichment()
#' c <- enrichment()
#' is.enrichment(list())                # FALSE
#' is.enrichment(1)                     # FALSE
#' is.enrichment(a)                     # TRUE
#' is.enrichment(c(a, c))               # TRUE TRUE
#' is.enrichment(list(a, b = "char"))   # TRUE FALSE
#' is.enrichment(c(a, b = list(12, c))) # TRUE FALSE TRUE
#'
methods::setGeneric(
  name = "is.enrichment",
  def = function(object) standardGeneric("is.enrichment")
)
methods::setMethod(f = "is.enrichment", signature = "ANY", definition = function(object) {
  if (length(object) > 1) {
    sapply(object, is.enrichment)
  } else {
    class(object) == "Enrichment"
  }
})
