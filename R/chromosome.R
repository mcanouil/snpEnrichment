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
