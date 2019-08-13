setClass(
  Class = "EnrichSNP",
  representation = representation(
    List = "character",
    Table = "matrix",
    EnrichmentRatio = "numeric",
    Z = "numeric",
    PValue = "numeric",
    Resampling = "matrix"
  ),
  prototype = prototype(
    List = character(),
    Table = matrix(0, ncol = 2, nrow = 2),
    EnrichmentRatio = numeric(),
    Z = numeric(),
    PValue = numeric(),
    Resampling = matrix(0, ncol = 5, nrow = 0)
  )
)


setMethod(f = "enrichSNP", signature = "ANY", definition = function(List, Table, EnrichmentRatio, Z, PValue, Resampling) {
  if (missing(List)) List <- character()
  if (missing(Table)) Table <- matrix(0, ncol = 2, nrow = 2)
  if (missing(EnrichmentRatio)) EnrichmentRatio <- numeric()
  if (missing(Z)) Z <- numeric()
  if (missing(PValue)) PValue <- numeric()
  if (missing(Resampling)) Resampling <- matrix(0, ncol = 5, nrow = 0)
  new(
    "EnrichSNP",
    List = List,
    Table = Table,
    EnrichmentRatio = EnrichmentRatio,
    Z = Z,
    PValue = PValue,
    Resampling = Resampling
  )
})


setMethod(f = "is.EnrichSNP", signature = "ANY", definition = function(object) {
  if (length(object) > 1) {
    sapply(object, is.EnrichSNP)
  } else {
    class(object) == "EnrichSNP"
  }
})


setMethod(f = "print", signature = "EnrichSNP", definition = function(x) {
  EnrichmentRatio <- x@EnrichmentRatio
  Z <- x@Z
  PValue <- x@PValue
  Resampling <- nrow(x@Resampling)
  Data <- sum(x@Table)
  List <- length(x@List)
  resTmp <- c(
    if (length(EnrichmentRatio) == 0) NA else EnrichmentRatio,
    if (length(Z) == 0) NA else Z,
    if (length(PValue) == 0) NA else PValue,
    Resampling,
    Data,
    List
  )
  names(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", "eSNP")
  res <- t(resTmp)
  rownames(res) <- "eSNP"
  res
})


.EnrichSNP.show <- function(object) {
  cat(
    "   - List :",
    paste0("(", length(object@List), ")"),
    if (length(object@List) <= 5) {
      if (length(object@List) == 0) "NA" else object@List
    } else {
      paste(paste(object@List[seq(5)], collapse = " "), "...")
    }
  )
  cat("\n   - Table :", paste0("(", paste(dim(object@Table), collapse = "x"), ")"))
  if (all(object@Table == 0)) {
    cat(" NA")
  } else {
    cat("\n")
    resFormat <- cbind(c("", rownames(object@Table)), rbind(colnames(object@Table), apply(round(object@Table, digits = 4), 2, as.character)))
    cat(paste("     ", apply(apply(matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)), 2, format, justify = "centre"), 1, paste, collapse = ""), "\n", sep = "", collapse = ""))
  }
  cat("\n   - EnrichmentRatio :", ifelse(length(object@EnrichmentRatio) == 0, NA, object@EnrichmentRatio))
  cat("\n   - Z               :", ifelse(length(object@Z) == 0, NA, object@Z))
  cat("\n   - PValue          :", ifelse(length(object@PValue) == 0, NA, ifelse((object@PValue == 0 & nrow(object@Resampling) != 0), paste0("<", 1 / nrow(object@Resampling)), object@PValue)))
  cat("\n   - Resampling      :", paste0("(", paste(dim(object@Resampling), collapse = "x"), ")"), ifelse(nrow(object@Resampling) == 0, 0, nrow(object@Resampling)))
  cat("\n")
  invisible(object)
}
setMethod(f = "show", signature = "EnrichSNP", definition = function(object) {
  cat("     ~~~ Class:", class(object), "~~~\n")
  .EnrichSNP.show(object)
  invisible(object)
})


setMethod(f = "[", signature = "EnrichSNP", definition = function(x, i, j, drop) {
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


setMethod(f = "[<-", signature = "EnrichSNP", definition = function(x, i, j, value) {
  switch(EXPR = i,
    "List" = x@List <- value,
    "Table" = x@Table <- value,
    "EnrichmentRatio" = x@EnrichmentRatio <- value,
    "Z" = x@Z <- value,
    "PValue" = x@PValue <- value,
    "Resampling" = x@Resampling <- value,
    stop("[EnrichSNP:set] ", i, ' is not a "EnrichSNP" slot.', call. = FALSE)
  )
  validObject(x)
  x
})


setMethod(f = "reset", signature = "EnrichSNP", definition = function(object, i) {
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
