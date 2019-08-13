#' Exclude SNPs from Enrichment analysis
#'
#' Remove a specify set of SNPs and compute a new enrichment analysis.
#'
#'
#' @name excludeSNP
#' @aliases excludeSNP excludeSNP-methods excludeSNP,Enrichment-method
#' excludeSNP,ANY-method
#' @docType methods
#' @param object [Enrichment]: an \code{\linkS4class{Enrichment}} object filled
#' by \code{\link{reSample}}.
#' @param excludeFile [vector(character)]: a list of SNPs to remove from a
#' previous enrichment analysis. A path to a file which the first column are
#' the SNPs.
#' @param mc.cores [numeric]: the number of cores to use (default is
#' \code{mc.cores=1}), i.e. at most how many child processes will be run
#' simultaneously.  Must be at least one, and parallelization requires at least
#' two cores.
#' @return Return the object given in argument where lists of SNPs are updated
#' by removing SNPs in \code{excludeFile}.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords excludeSNP methods
#' @examples
#'
#' \dontrun{data(toyEnrichment)
#' excludeFile <- c(
#'     "rs4376885", "rs17690928", "rs6460708", "rs13061537", "rs11769827",
#'     "rs12717054", "rs2907627", "rs1380109", "rs7024214", "rs7711972",
#'     "rs9658282", "rs11750720", "rs1793268", "rs774568", "rs6921786",
#'     "rs1699031", "rs6994771", "rs16926670", "rs465612", "rs3012084",
#'     "rs354850", "rs12803455", "rs13384873", "rs4364668", "rs8181047",
#'     "rs2179993", "rs12049335", "rs6079926", "rs2175144", "rs11564427",
#'     "rs7786389", "rs7005565", "rs17423335", "rs12474102", "rs191314",
#'     "rs10513168", "rs1711437", "rs1992620", "rs283115", "rs10754563",
#'     "rs10851727", "rs2173191", "rs7661353", "rs1342113", "rs7042073",
#'     "rs1567445", "rs10120375", "rs550060", "rs3761218", "rs4512977"
#' )
#' # OR
#' excludeFile <- system.file("extdata/Exclude/toyExclude.txt",
#'                            package = "snpEnrichment")
#'
#' toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)
#' toyEnrichment_exclude}
#'
methods::setGeneric(
  name = "excludeSNP",
  def = function(object, excludeFile, mc.cores = 1) standardGeneric("excludeSNP")
)
methods::setMethod(f = "excludeSNP", signature = "ANY", definition = function(object, excludeFile, mc.cores = 1) {
  if (!is.enrichment(object)) {
    stop('[Method:excludeSNP] not available for "', class(object), '" object.', call. = FALSE)
  }
})
methods::setMethod(f = "excludeSNP", signature = "Enrichment", definition = function(object, excludeFile, mc.cores = 1) {
  if (missing(excludeFile)) {
    stop('[Enrichment:excludeSNP] argument "excludeFile" is missing.', call. = FALSE)
  }
  cat("########## Exclude SNP list Start ##########\n")
  if (all(class(try(close(file(excludeFile)), silent = TRUE)) != "try-error")) {
    eSNPexclude <- utils::read.delim(file = excludeFile, header = FALSE, na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, stringsAsFactors = FALSE, sep = "\t")
  } else {
    eSNPexclude <- excludeFile
  }
  if (class(eSNPexclude) %in% c("matrix", "data.frame")) {
    if (ncol(eSNPexclude) > 1) {
      eSNPexclude <- eSNPexclude[, 1]
    }
  }
  eSNPexclude <- unlist(eSNPexclude, use.names = FALSE)
  callLD <- object@Call$readEnrichment$LD
  resParallel <- mclapply2(X = seq_len(22), mc.cores = min(22, mc.cores), FUN = function(iChr) {
    chrObject <- eval(parse(text = paste0("object@Chromosomes$Chrom", iChr)))
    temp <- chrObject@Data
    if (callLD) {
      xSNPexclude <- intersect(temp[, "SNP"], unique(c(eSNPexclude, chrObject@LD[names(chrObject@LD) %in% eSNPexclude])))
    } else {
      xSNPexclude <- eSNPexclude
    }
    temp[temp[, "SNP"] %in% xSNPexclude, "eSNP"] <- 0
    temp[temp[, "SNP"] %in% xSNPexclude, "xSNP"] <- 0
    chromosome(Data = temp, LD = chrObject@LD)
  })
  names(resParallel) <- paste0("Chrom", seq_len(22))
  result <- enrichment(
    Loss = object@Loss,
    Chromosomes = resParallel,
    Call = list(
      readEnrichment = object@Call$readEnrichment,
      reSample = list(
        object = NULL, nSample = NULL, empiricPvalue = NULL,
        MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL
      )
    )
  )
  rm(resParallel)
  GC()

  result <- computeER(object = result, sigThresh = object@Call$readEnrichment$sigThresh, mc.cores = mc.cores)
  cat("########### Update SNP list END ############\n")
  for (iType in c("eSNP", "xSNP")) {
    cat("   ", length(setdiff(object[iType]@List, result[iType]@List)), " SNPs are removed from", iType, "list.\n")
  }
  result@Loss <- cbind(result@Loss,
    exclude = c(
      result@Loss["Signal", "CIS"],
      length(result["List", seq_len(22)][["eSNP"]]),
      sapply(seq_len(22), function(iChr) {
        length(result["List", iChr][["eSNP"]])
      })
    )
  )

  cat("########### Exclude SNP list END ###########\n")
  result
})
