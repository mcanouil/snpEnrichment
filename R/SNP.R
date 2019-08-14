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
