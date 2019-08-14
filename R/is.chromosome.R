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
