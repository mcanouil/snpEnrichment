#' ~ Overview: SNPs Enrichment Analysis With Sampling ~
#'
#' @details
#'   Implements classes and methods for large-scale SNP enrichment analysis
#'   (*e.g.*, SNPs associated with genes expression in a GWAS signal).
#'
#' @note Internal data management in [snpEnrichment] use RefSNP (rs) IDs.
#'
#' @seealso
#' + Overview: [snpEnrichment]
#' + Classes: [Enrichment-class]
#' + Methods:
#'     - [plot]
#'     - [reSample]
#'     - [getEnrichSNP]
#'     - [excludeSNP]
#'     - [compareEnrichment]
#'     - [enrichment]
#'     - [is.enrichment]
#' + Functions:
#'     - [initFiles]
#'     - [writeLD]
#'     - [readEnrichment]
#'
#' @docType package
#' @name snpEnrichment
#'
"_PACKAGE"
