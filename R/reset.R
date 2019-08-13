methods::setGeneric(
  name = "reset",
  def = function(object, i) standardGeneric("reset")
)

methods::setMethod(f = "reset", signature = "ANY", definition = function(object, i) {
  if (!(is.enrichment(object) & is.chromosome(object))) {
    stop('[Method:reset] not available for "', class(object), '" object.', call. = FALSE)
  }
})

methods::setMethod(f = "reset", signature = "EnrichSNP", definition = function(object, i) {
  switch(EXPR = i,
    "List" = object@List <- character(),
    "Table" = object@Table <- matrix(0, ncol = 2, nrow = 2),
    "EnrichmentRatio" = object@EnrichmentRatio <- numeric(),
    "Z" = object@Z <- numeric(),
    "PValue" = object@PValue <- numeric(),
    "Resampling" = object@Resampling <- matrix(0, ncol = 5, nrow = 0),
    stop("[EnrichSNP:reset] ", i, ' is not a "EnrichSNP" slot.', call. = FALSE)
  )
  object
})

methods::setMethod(f = "reset", signature = "Chromosome", definition = function(object, i) {
  switch(EXPR = i,
    "Data" = object@Data <- data.frame(),
    "LD" = object@LD <- character(),
    "eSNP" = object@eSNP <- enrichSNP(),
    "xSNP" = object@xSNP <- enrichSNP(),
    "List" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@List <- reset(object[type], "List")
      }
    },
    "Table" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@Table <- matrix(0, ncol = 2, nrow = 2)
      }
    },
    "EnrichmentRatio" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@EnrichmentRatio <- numeric()
      }
    },
    "Z" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@Z <- numeric()
      }
    },
    "PValue" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@PValue <- numeric()
      }
    },
    "Resampling" = {
      for (type in c("eSNP", "xSNP")) {
        object[type]@Resampling <- matrix(0, ncol = 5, nrow = 0)
      }
    },
    stop("[Enrichment:reset] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
  )
  object
})

methods::setMethod(f = "reset", signature = "Enrichment", definition = function(object, i) {
  switch(EXPR = i,
    "Loss" = object@Loss <- data.frame(),
    "Data" = object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Data"),
    "LD" = object@Chromosomes <- lapply(object@Chromosomes, reset, i = "LD"),
    "Call" = {
      object@Call <- list(
        readEnrichment = list(
          pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL, snpInfoDir = NULL,
          distThresh = NULL, sigThresh = NULL, LD = NULL, ldDir = NULL, mc.cores = NULL
        ),
        reSample = list(
          object = NULL, nSample = NULL, empiricPvalue = NULL,
          MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL
        )
      )
    },
    "eSNP" = {
      object@eSNP <- enrichSNP()
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "eSNP")
    },
    "xSNP" = {
      object@xSNP <- enrichSNP()
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "xSNP")
    },
    "List" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@List <- character()
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "List")
    },
    "Table" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@Table <- matrix(0, ncol = 2, nrow = 2)
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Table")
    },
    "EnrichmentRatio" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@EnrichmentRatio <- numeric()
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "EnrichmentRatio")
    },
    "Z" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@Z <- numeric()
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Z")
    },
    "PValue" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@PValue <- numeric()
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "PValue")
    },
    "Resampling" = {
      for (iType in c("eSNP", "xSNP")) {
        object[iType]@Resampling <- matrix(0, ncol = 5, nrow = 0)
      }
      object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Resampling")
    },
    "Chromosomes" = {
      object@Chromosomes <- eval(parse(text = paste0("list(", paste(paste0("Chrom", seq_len(22), " = chromosome()"), collapse = ", "), ")")))
    },
    stop("[Enrichment:reset] ", i, ' is not a "Enrichment" slot.', call. = FALSE)
  )
  object
})
