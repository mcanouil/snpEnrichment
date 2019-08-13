methods::setGeneric(
  name = "is.EnrichSNP",
  def = function(object) standardGeneric("is.EnrichSNP")
)
methods::setMethod(f = "is.EnrichSNP", signature = "ANY", definition = function(object) {
  if (length(object) > 1) {
    sapply(object, is.EnrichSNP)
  } else {
    class(object) == "EnrichSNP"
  }
})
