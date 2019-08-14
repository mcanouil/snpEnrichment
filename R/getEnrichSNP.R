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
