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
