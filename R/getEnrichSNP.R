#' Get all eSNP/xSNP which are enriched
#'
#' \code{\link{getEnrichSNP}} get all eSNP/xSNP in a
#' \code{\linkS4class{Enrichment}} object which are significant in the signal
#' according to \code{sigThresh} defined in \code{\link{readEnrichment}}.
#'
#'
#' @name getEnrichSNP
#' @aliases getEnrichSNP getEnrichSNP-methods getEnrichSNP,Enrichment-method
#' getEnrichSNP,ANY-method
#' @docType methods
#' @param object [Enrichment]: an object of class
#' \code{\linkS4class{Enrichment}}.
#' @param type [character]: extract \code{eSNP} or \code{xSNP} data.
#' @return Return a \code{data.frame} with eSNP/xSNP which are enriched in
#' signal given to \code{signalFile} in function \code{\link{readEnrichment}}.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords getEnrichSNP methods
#' @examples
#'
#' \dontrun{data(toyEnrichment)
#' eSNPenriched <- getEnrichSNP(object = toyEnrichment, type = "eSNP")
#' head(eSNPenriched)}
#'
methods::setGeneric(
  name = "getEnrichSNP",
  def = function(object, type = "eSNP") standardGeneric("getEnrichSNP")
)

methods::setMethod(f = "getEnrichSNP", signature = "ANY", definition = function(object, type = "eSNP") {
  if (!(is.enrichment(object))) {
    stop('[Method:getEnrichSNP] not available for "', class(object), '" object.', call. = FALSE)
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
