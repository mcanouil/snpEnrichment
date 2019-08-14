#' Toy dataset with SNP data
#'
#' This data set gives an [Enrichment-class] object after, [initFiles] and [readEnrichment] is ran.
#' Compute LD for all SNPs in `snpListDir` files two by two. Genome Build 37.3 (hg19).
#'
#' @format See class [Enrichment-class] for details about the format.
"toyEnrichment"

#' Transcript information in order to check the CIS status for SNPs
#'
#' This dataset is used by [readEnrichment] and [compareEnrichment] in order to check the CIS status
#' for each SNP of signal. Genome Build 37.3 (hg19).
#'
#' @format See class [readEnrichment] and [compareEnrichment] for details about how to use this dataset.
#'
#' @name transcript
NULL
