methods::setMethod(f = "show", signature = "EnrichSNP", definition = function(object) {
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
  cat("     ~~~ Class:", class(object), "~~~\n")
  .EnrichSNP.show(object)
  invisible(object)
})

methods::setMethod(f = "show", signature = "Chromosome", definition = function(object) {
  .Chromosome.show <- function(object) {
    cat("  ~ Data :", paste0("(", paste(dim(object@Data), collapse = "x"), ")"))
    nrowShow <- seq_len(min(5, nrow(object@Data)))
    ncolShow <- seq_len(min(10, ncol(object@Data)))
    if (nrow(object@Data) == 0) {
      cat(" NA")
    } else {
      cat("\n")
      resFormat <- cbind(c("", rownames(object@Data[nrowShow, ncolShow])), rbind(colnames(object@Data[nrowShow, ncolShow]), apply(object@Data[nrowShow, ncolShow], 2, as.character)))
      resFormat <- rbind(resFormat, ".....")
      cat(paste("     ", apply(apply(matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)), 2, format, justify = "centre"), 1, paste, collapse = ""), "\n", sep = "", collapse = ""))
    }
    cat("\n  ~ LD :")
    if (length(object@LD) == 0) {
      cat(" NA")
    } else {
      cat("\n")
      if (length(object@LD) > 5) {
        tmpLD <- object@LD[seq_len(5)]
      } else {
        tmpLD <- object@LD
      }
      resFormat <- cbind(c("SNP1", "SNP2"), rbind(names(tmpLD), tmpLD))
      resFormat <- cbind(resFormat, "...")
      cat(paste("     ", apply(apply(matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)), 2, format, justify = "centre"), 1, paste, collapse = ""), "\n", sep = "", collapse = ""))
    }
    cat("\n  ~ eSNP :\n")
    .EnrichSNP.show(object@eSNP)
    cat("\n  ~ xSNP :\n")
    .EnrichSNP.show(object@xSNP)
    cat("\n")
    invisible(object)
  }
  cat("     ~~~ Class:", class(object), "~~~\n")
  .Chromosome.show(object)
  invisible(object)
})

methods::setMethod(f = "show", signature = "Enrichment", definition = function(object) {
  .Enrichment.show <- function(object) {
    .showArgs <- function(args) {
      for (iFuncArg in seq(args)) {
        if (is.null(unlist(args[iFuncArg], use.names = FALSE))) {
          cat(paste0("    ", names(args[iFuncArg]), "() : Not yet called."))
        } else {
          tmpArgs <- args[[names(args[iFuncArg])]]
          type <- lapply(tmpArgs, class)
          res <- NULL
          tmpArgs <- lapply(tmpArgs, function(li) if (length(li) > 1) deparse(li) else li)
          for (iArg in names(tmpArgs)) {
            res <- c(
              res,
              paste0(
                iArg,
                ifelse(type[[iArg]] == "character", '="', "="),
                ifelse(
                  test = type[[iArg]] == "NULL" | type[[iArg]] == "name",
                  yes = deparse(tmpArgs[[iArg]]),
                  no = tmpArgs[[iArg]]
                ),
                ifelse(type[[iArg]] == "character", '"', "")
              )
            )
          }
          cat(paste0(
            "    ", names(args[iFuncArg]), "(",
            paste(
              res,
              collapse = paste0(", \n", paste(rep(" ", nchar(names(args[iFuncArg])) + 5), collapse = ""))
            ),
            ")"
          ))
        }
      }
    }
    cat(" ~ Loss :", paste0("(", paste(dim(object@Loss), collapse = "x"), ")"))
    if (nrow(object@Loss) == 0) {
      cat(" NA")
    } else {
      cat("\n")
      resFormat <- cbind(
        c("", rownames(object@Loss[c(1, 2), ])),
        rbind(
          colnames(object@Loss[c(1, 2), ]),
          apply(round(object@Loss[c(1, 2), ], digits = 4), 2, as.character)
        )
      )
      resFormat <- rbind(resFormat, ".....")
      cat(paste(
        "     ",
        apply(
          X = apply(
            X = matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)),
            MARGIN = 2,
            FUN = format, justify = "centre"
          ),
          MARGIN = 1,
          FUN = paste, collapse = ""
        ),
        "\n",
        sep = "",
        collapse = ""
      ))
    }
    cat("\n ~ Call :\n")
    .showArgs(object@Call)
    cat("\n ~ eSNP :\n")
    .EnrichSNP.show(object@eSNP)
    cat("\n ~ xSNP :\n")
    .EnrichSNP.show(object@xSNP)
    nbVoid <- 0
    for (slot in names(object@Chromosomes)) {
      eval(parse(text = paste0(
        "nbVoid <- nbVoid + as.numeric(identical(chromosome(), ", paste0("object@Chromosomes$", slot), "))"
      )))
    }
    cat("\n ~ Chromosomes :", paste(nbVoid, length(names(object@Chromosomes)), sep = "/"), "empty Chromosomes")
    cat("\n")
    invisible(object)
  }
  cat("    ~~~ Class:", class(object), "~~~ \n")
  .Enrichment.show(object)
  invisible(object)
})
