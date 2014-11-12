

#' Class \code{"\linkS4class{EnrichSNP}"}
#' 
#' This class is defined to summarize the enrichment analysis. It's a part of
#' \code{\linkS4class{Chromosome}} and \code{\linkS4class{Enrichment}} classes.
#' 
#' 
#' @name EnrichSNP-class
#' @aliases EnrichSNP-class EnrichSNP [,EnrichSNP-method
#' [,EnrichSNP,ANY,ANY,ANY-method [<-,EnrichSNP-method
#' [<-,EnrichSNP,ANY,ANY,ANY-method show,EnrichSNP-method
#' print,EnrichSNP-method
#' @docType class
#' @note \code{\linkS4class{EnrichSNP}} object is not intended to be use
#' directly by user. It is a part of the \code{\linkS4class{Enrichment}} and
#' \code{\linkS4class{Chromosome}} object.
#' @section Slots: \describe{ \item{List}{[vector(character)]: a list of SNPs
#' used to compute enrichment (e.g. eSNP or xSNP).} \item{Table}{[matrix]:
#' Contingency table with SNPs (columns) and P-Values from signal (rows).}
#' \item{EnrichmentRatio}{[numeric]: Enrichment Ratio is computed on the
#' contingency table (\code{Table} slot).} \item{Z}{[numeric]: A statistic
#' computed from \code{EnrichmentRatio} and resampling results.}
#' \item{PValue}{[numeric]: P-Value associated with the statistic \code{Z}.}
#' \item{Resampling}{[matrix]: A matrix with by row, the contingency table and
#' the odds ratio for each resampling.} }
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords classes class enrichSNP
NULL





#' Class \code{\linkS4class{Chromosome}}
#' 
#' This class is defined to summarize the enrichment analysis about a
#' chromosome.
#' 
#' 
#' @name Chromosome-class
#' @aliases Chromosome-class Chromosome [,Chromosome-method
#' [,Chromosome,ANY,ANY,ANY-method [<-,Chromosome-method
#' [<-,Chromosome,ANY,ANY,ANY-method show,Chromosome-method chromosome
#' chromosome-methods chromosome,ANY-method
#' @docType class
#' @note \code{\linkS4class{Chromosome}} object is not intended to be used
#' alone on this version.\ It is a part of the \code{\linkS4class{Enrichment}}
#' object.
#' @section Objects from the Class: \code{\link{chromosome}} is defined to
#' build an object of class \code{\linkS4class{Chromosome}} in order to compute
#' an enrichment analysis.  A \code{\linkS4class{Chromosome}} object contains
#' the original data, a list of SNPs, some results and resampling data.\cr
#' 
#' When created, an \code{\linkS4class{Chromosome}} object is "empty".
#' \code{\link{readEnrichment}} initializes a \code{\linkS4class{Chromosome}}
#' object with value from PLINK computation and user's files.  In this step,
#' only the fields "Data", "LD", "SNP" are filled.  \code{reSample} fills the
#' fields: Table, EnrichmentRatio, Z, PValue and Resampling of a
#' \code{\linkS4class{Chromosome}}.
#' 
#' Note that if \code{\link{reSample}} is executed on an
#' \code{\linkS4class{Chromosome}} every new resampling is added to the
#' original ones, pre-existing statistics are erased and computed again with
#' the new resampling set.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords classes class chromosome
#' @examples
#' 
#' Data <- structure(
#'     list(
#'         SNP = c("rs4970420", "rs3766178",
#'                 "rs10910030", "rs10910036",
#'                 "rs2234167", "rs6661861"),
#'         PVALUE = c(0.9244, 0.167, 0.01177, 0.4267, 0.9728, 0.4063),
#'         CHR = c(1, 1, 1, 1, 1, 1),
#'         POS = c(1106473, 1478180, 2035684, 2183754, 2494330, 3043121),
#'         MAF = c(0.2149, 0.3117, 0.374, 0.3753, 0.1565, 0.06101),
#'         eSNP = c(0, 1, 1, 0, 0, 0),
#'         xSNP = c(0, 1, 1, 0, 0, 0)
#'     ),
#'     .Names = c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP", "xSNP"),
#'     row.names = c("rs4970420", "rs3766178",
#'                   "rs10910030", "rs10910036",
#'                   "rs2234167", "rs6661861"),
#' class = "data.frame")
#' 
#' toyChr <- chromosome(Data = Data)
#' show(toyChr)
#' toyChr
#' 
#' toyChr <- chromosome()
#' toyChr["Data"] <- Data
#' toyChr
#' 
#' is.chromosome(toyChr) # TRUE
#' 
NULL





#' Class \code{\linkS4class{Enrichment}}
#' 
#' This class is defined to summarize the enrichment analysis on each
#' chromosomes and the whole genome.
#' 
#' 
#' @name Enrichment-class
#' @aliases Enrichment-class Enrichment [,Enrichment-method
#' [,Enrichment,ANY,ANY,ANY-method [<-,Enrichment-method
#' [<-,Enrichment,ANY,ANY,ANY-method show,Enrichment-method enrichment
#' enrichment-methods enrichment,ANY-method
#' @docType class
#' @section Objects from the Class: \code{\link{enrichment}} is defined to
#' build an object of class \code{\linkS4class{Enrichment}} in order to compute
#' an enrichment analysis.  \code{Enrichment} is the object containing the
#' results for all \code{\linkS4class{Chromosome}} object and for the whole
#' genome.\cr
#' 
#' When an \code{\linkS4class{Enrichment}} object is created, it contains a
#' list of SNPs (e.g. eSNPs).  All the others slots are "empty". After
#' \code{\link{reSample}} is ran on an \code{\linkS4class{Enrichment}} object,
#' the slots: Table, EnrichmentRatio, Z, PValue and Resampling are filled.
#' 
#' Note that if \code{\link{reSample}} is executed on an
#' \code{\linkS4class{Enrichment}} every new resampling is added to the
#' original ones, pre-existing statistics are erased and computed again with
#' the new resampling set.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords classes class enrichment
#' @examples
#' 
#' data(toyEnrichment)
#' toyEnrich <- enrichment()
#' show(toyEnrich)
#' 
#' toyEnrich["Loss"] <- toyEnrichment["Loss"]
#' toyEnrich["Loss"]
#' 
#' toyEnrich <- enrichment(Loss = toyEnrichment["Loss"],
#'                         eSNP = toyEnrichment["eSNP"])
#' toyEnrich <- enrichment(Loss = toyEnrichment["Loss"])
#' 
#' \dontrun{reSample(object = toyEnrichment,
#'          nSample = 10,
#'          empiricPvalue = TRUE,
#'          MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
#'          mc.cores = 1,
#'          onlyGenome = TRUE)
#' print(toyEnrichment)
#' 
#' excludeFile <- c(
#' "rs7897180", "rs4725479", "rs315404", "rs17390391", "rs1650670",
#'     "rs6783390", "rs1642009", "rs4756586", "rs11995037", "rs4977345",
#'     "rs13136448", "rs4233536", "rs11151079", "rs2299657", "rs4833930",
#'     "rs1384", "rs7168184", "rs6909895", "rs7972667", "rs2293229",
#'     "rs918216", "rs6040608", "rs2817715", "rs13233541", "rs4486743",
#'     "rs2127806", "rs10912854", "rs1869052", "rs9853549", "rs448658",
#'     "rs2451583", "rs17483288", "rs10962314", "rs9612059", "rs1384182",
#'     "rs8049208", "rs12215176", "rs2980996", "rs1736976", "rs8089268",
#'     "rs10832329", "rs12446540", "rs7676237", "rs869922", "rs16823426",
#'     "rs1374393", "rs13268781", "rs11134505", "rs7325241", "rs7520109"
#' )
#' 
#' toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)
#' print(toyEnrichment_exclude)}
#' 
NULL





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
NULL





#' Get all eSNP/xSNP which are enriched
#' 
#' \code{\link{getEnrichSNP}} get all eSNP/xSNP in a
#' \code{\linkS4class{Enrichment}} object which are significant in the signal
#' according to \code{sigThresh} defined in \code{\link{readEnrichment}}.
#' 
#' 
#' @name getEnrichSNP
#' @aliases getEnrichSNP getEnrichSNP-methods getEnrichSNP,Enrichment-method
#' getEnrichSNP,ANY-method
#' @docType methods
#' @param object [Enrichment]: an object of class
#' \code{\linkS4class{Enrichment}}.
#' @param type [character]: extract \code{eSNP} or \code{xSNP} data.
#' @return Return a \code{data.frame} with eSNP/xSNP which are enriched in
#' signal given to \code{signalFile} in function \code{\link{readEnrichment}}.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords getEnrichSNP methods
#' @examples
#' 
#' \dontrun{data(toyEnrichment)
#' eSNPenriched <- getEnrichSNP(object = toyEnrichment, type = "eSNP")
#' head(eSNPenriched)}
#' 
NULL





#' Is an Chromosome object
#' 
#' 'is.chromosome' returns 'TRUE' if 'x' is an \code{\linkS4class{Chromosome}}
#' object and 'FALSE' otherwise.
#' 
#' 
#' @aliases is.chromosome is.chromosome-methods is.chromosome,ANY-method
#' @param object [ANY]: object to be tested.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords chromosome is is.chromosome
#' @examples
#' 
#' a <- chromosome()
#' c <- chromosome()
#' is.chromosome(list())                # FALSE
#' is.chromosome(1)                     # FALSE
#' is.chromosome(a)                     # TRUE
#' is.chromosome(c(a, c))               # TRUE TRUE
#' is.chromosome(list(a, b = "char"))   # TRUE FALSE
#' is.chromosome(c(a, b = list(12, c))) # TRUE FALSE TRUE
#' 
NULL





#' Is an Enrichment object
#' 
#' 'is.enrichment' returns 'TRUE' if 'x' is an \code{\linkS4class{Enrichment}}
#' object and 'FALSE' otherwise.
#' 
#' 
#' @aliases is.enrichment is.enrichment-methods is.enrichment,ANY-method
#' @param object [ANY]: object to be tested.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords enrichment is is.enrichment
#' @examples
#' 
#' a <- enrichment()
#' c <- enrichment()
#' is.enrichment(list())                # FALSE
#' is.enrichment(1)                     # FALSE
#' is.enrichment(a)                     # TRUE
#' is.enrichment(c(a, c))               # TRUE TRUE
#' is.enrichment(list(a, b = "char"))   # TRUE FALSE
#' is.enrichment(c(a, b = list(12, c))) # TRUE FALSE TRUE
#' 
NULL





#' Plot method (S4) for \code{\linkS4class{Enrichment}} object
#' 
#' \code{\link{plot}} is a generic function for plotting of R objects. The
#' function invokes particular \code{methods} which depend on the \code{class}
#' of the first argument.
#' 
#' 
#' @name plot-methods
#' @aliases plot plot-methods plot,Enrichment-method plot,Enrichment,ANY-method
#' @docType methods
#' @param x [Enrichment]: an object of class \code{\linkS4class{Enrichment}}
#' which the Z statistics or p-values have to be drawn.
#' @param what [character or vector(numeric)]: default \code{what="Genome"})
#' plot Z statistics or p-values for genome only (what must be: \code{"All"},
#' \code{"Genome"} or numeric vector).
#' @param type [vector(character)]: plot the selected analysis for
#' \code{"eSNP"} and/or \code{"xSNP"}.
#' @param ggplot [logical]: use ggplot (default \code{ggplot=FALSE}) instead of
#' classic plot method.
#' @param pvalue [logical]: if \code{TRUE}, p-value convergense is plotted.
#' Otherwise, Z statistic is plotted.
#' @param ... [any]: Arguments to be passed to methods, such as graphical
#' parameters (see \code{par})
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords plot methods
#' @examples
#' 
#' \dontrun{data(toyEnrichment)
#' reSample(toyEnrichment, 10)
#' plot(toyEnrichment)}
#' 
NULL





#' Print method (S4)
#' 
#' \code{\link{print}} is a generic function used to print results.
#' 
#' 
#' @name print-methods
#' @aliases print print-methods print,Chromosome-method print,Enrichment-method
#' @docType methods
#' @param x [Enrichment or Chromosome]: an object of class
#' \code{\linkS4class{Enrichment}} or \code{\linkS4class{Chromosome}}.
#' @param what [character or numeric]: \code{what="Genome"} (default) to print
#' results as a matrix. \code{what} could be \code{"All"}, \code{"Genome"} or a
#' numeric from 1 to 22 (numeric vector is allowed).
#' @param type [character]: select if results for \code{"eSNP"} and/or
#' \code{"xSNP"} should be print.
#' @return Return a \code{matrix} for classes \code{\linkS4class{Enrichment}}
#' and \code{\linkS4class{Chromosome}}.\cr
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords print methods
#' @examples
#' 
#' data(toyEnrichment)
#' print(toyEnrichment, "All", type = "eSNP")
#' print(toyEnrichment, "Genome")
#' print(toyEnrichment, 1)
#' 
NULL





#' Compare enrichment analysis between two SNPs list
#' 
#' Compare the enrichment analysis between two set of SNPs.
#' \code{\link{compareEnrichment}} compare two \code{\linkS4class{Enrichment}}
#' objects.
#' 
#' 
#' @name compareEnrichment
#' @aliases compareEnrichment compareEnrichment-methods
#' compareEnrichment,ANY-method
#' compareEnrichment,Enrichment,Enrichment,ANY-method
#' @docType methods
#' @param object.x,object.y [Enrichment]: an \code{\linkS4class{Enrichment}}
#' object fully filled (e.g. \code{\link{readEnrichment}}).
#' @param pattern [character]: character string containing a expression to be
#' matched with all chromosomes files (e.g."Chrom" for files which start by
#' "Chrom" followed by the chromosome number).
#' @param nSample [numeric]: the number of resampling done by
#' \code{\link{reSample}} for p-values computation (minimum is 100).
#' @param empiricPvalue [logical]: \code{empiricPvalue=TRUE} (default) compute
#' PValue based on the null distribution (resampling).  If
#' \code{empiricPvalue=TRUE}, the empirical p-values are computed instead.
#' @param mc.cores [numeric]: the number of cores to use (default is
#' \code{mc.cores=1}), i.e. at most how many child processes will be run
#' simultaneously.  Must be at least one, and parallelization requires at least
#' two cores.
#' @param onlyGenome [logical]: \code{onlyGenome=TRUE} (default) compute
#' resampling step for all chromosomes.
#' @return Return a \code{list} of three elements:
#' \item{object.xy}{\code{\linkS4class{Enrichment}} object from the comparison
#' between \code{object.x} and \code{object.y}.}
#' \item{object.x}{\code{\linkS4class{Enrichment}} object passed in
#' \code{object.x} with resampling data.}
#' \item{object.y}{\code{\linkS4class{Enrichment}} object passed in
#' \code{object.y} with resampling data.}
#' @note Still in development.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords compareEnrichment Enrichment
#' @examples
#' 
#' \dontrun{data(toyEnrichment)
#' 
#' reSample(object = toyEnrichment,
#'          nSample = 10,
#'          empiricPvalue = TRUE,
#'          MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
#'          mc.cores = 1,
#'          onlyGenome = TRUE)
#' 
#' excludeFile <- c(
#'     "rs7897180", "rs4725479", "rs315404", "rs17390391", "rs1650670",
#'     "rs6783390", "rs1642009", "rs4756586", "rs11995037", "rs4977345",
#'     "rs13136448", "rs4233536", "rs11151079", "rs2299657", "rs4833930",
#'     "rs1384", "rs7168184", "rs6909895", "rs7972667", "rs2293229",
#'     "rs918216", "rs6040608", "rs2817715", "rs13233541", "rs4486743",
#'     "rs2127806", "rs10912854", "rs1869052", "rs9853549", "rs448658",
#'     "rs2451583", "rs17483288", "rs10962314", "rs9612059", "rs1384182",
#'     "rs8049208", "rs12215176", "rs2980996", "rs1736976", "rs8089268",
#'     "rs10832329", "rs12446540", "rs7676237", "rs869922", "rs16823426",
#'     "rs1374393", "rs13268781", "rs11134505", "rs7325241", "rs7520109"
#' )
#' # OR
#' excludeFile <- system.file("extdata/Exclude/toyExclude.txt",
#'                            package = "snpEnrichment")
#' 
#' toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)
#' 
#' compareResults <- compareEnrichment(object.x = toyEnrichment,
#'                                     object.y = toyEnrichment_exclude,
#'                                     pattern = "Chrom",
#'                                     nSample = 10,
#'                                     empiricPvalue = FALSE,
#'                                     mc.cores = 1,
#'                                     onlyGenome = TRUE)}
#' 
NULL





#' ~ Internal: snpEnrichment objects and methods ~
#' 
#' Internal: snpEnrichment objects and methods.
#' 
#' These are not intended to be called by the user.\cr If you have a specific
#' need or for more details on snpEnrichment, feel free to contact the author.
#' 
#' @aliases .verbose .cleanTmpDir .showArgs .showArgs .writeSignal .writeFreq
#' .readSignal .readSNP .verbose .checkFilePath .showArgs .writeSignal
#' .writeFreq .readTranscript .readFreq .readLD .checkTranscript .readFiles
#' .enrichmentRatio .reSample .compareEnrich .EnrichSNP.show .Chromosome.show
#' .Enrichment.show .checkSnpInfoDir .checkSnpListDir .checkSignalFile
#' .splitByChrom enrichSNP enrichSNP-methods enrichSNP,ANY-method is.EnrichSNP
#' is.EnrichSNP-methods is.EnrichSNP,ANY-method reset reset-methods
#' reset,EnrichSNP-method reset,EnrichSNP,ANY-method reset,Chromosome-method
#' reset,Chromosome,ANY-method reset,Enrichment-method
#' reset,Enrichment,ANY-method reset,ANY-method computeER
#' computeER,Chromosome-method computeER,Enrichment-method computeER,ANY-method
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords snpEnrichment internal Enrichment
NULL





#' ~ Overview: SNPs enrichment analysis ~
#' 
#' Implements classes and methods for large-scale SNP enrichment analysis (e.g.
#' SNPs associated with genes expression in a GWAS signal).
#' 
#' \tabular{ll}{ Package: \tab snpEnrichment\cr Title: \tab SNPs enrichment
#' analysis\cr Author: \tab Mickael Canouil <mickael.canouil at good.ibl.fr>\cr
#' Contributor: \tab Loic Yengo\cr Maintainer: \tab Mickael Canouil
#' <mickael.canouil at good.ibl.fr>\cr % Description: \tab Implements classes
#' and methods for large-scale SNP enrichment analysis \cr (e.g. SNPs
#' associated with genes expression and a GWAS signal) \cr License: \tab GPL
#' (>= 2)\cr Depends: \tab R (>= 3.0.0), methods\cr Suggests: \tab grid,
#' ggplot2\cr Imports: \tab parallel, snpStats\cr URL: \tab
#' https://github.com/mcanouil/snpEnrichment\cr Encoding: \tab UTF-8\cr }
#' 
#' @name snpEnrichment-package
#' @aliases snpEnrichment-package snpEnrichment
#' @docType package
#' @note Internal data management in 'snpEnrichment' use RefSNP (rs) IDs.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords snpEnrichment package Enrichment
#' @examples
#' 
#' ###################
#' ### 1. Prepare data
#' \dontrun{snpInfoDir <- system.file("extdata/snpInfo",
#'                           package = "snpEnrichment")
#' signalFile <- system.file("extdata/Signal/toySignal.txt",
#'                           package = "snpEnrichment")
#' 
#' initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)
#' 
#' writeLD(pattern = "Chrom", snpInfoDir, signalFile,
#'         ldDir = NULL, ldThresh = 0.8, depth = 1000,
#'         mc.cores = 1)}
#' 
#' 
#' ################
#' ### 2. Read data
#' \dontrun{snpListDir <- system.file("extdata/List",
#'                           package = "snpEnrichment")
#' data(transcript)
#' transcriptFile <- transcript
#' 
#' toyData <- readEnrichment(pattern = "Chrom", signalFile,
#'                          transcriptFile, snpListDir,
#'                          snpInfoDir, distThresh = 1000,
#'                          sigThresh = 0.05, LD = TRUE,
#'                          ldDir = NULL, mc.cores = 1)
#' toyData}
#' 
#' 
#' ######################
#' ### 3. Compute results
#' \dontrun{reSample(object = toyData,
#'          nSample = 10,
#'          empiricPvalue = TRUE,
#'          MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
#'          mc.cores = 1,
#'          onlyGenome = TRUE)}
#' 
#' 
#' #######################
#' ### 4. Further analysis: Exclude SNP from original list.
#' \dontrun{excludeFile <- c(
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
#' toyData_exclude <- excludeSNP(toyData, excludeFile, mc.cores = 1)
#' 
#' # Warning: compareEnrichment is in development!!
#' compareResults <- compareEnrichment(object.x = toyData,
#'                                     object.y = toyData_exclude,
#'                                     pattern = "Chrom",
#'                                     nSample = 10,
#'                                     empiricPvalue = TRUE,
#'                                     mc.cores = 1,
#'                                     onlyGenome = TRUE)}
#' 
#' 
#' ####################
#' ### 5. Watch results
#' \dontrun{show(toyData)
#' print(toyData)
#' head(getEnrichSNP(toyData, type = "xSNP"))
#' 
#' show(toyData_exclude)
#' print(toyData_exclude)
#' head(getEnrichSNP(toyData_exclude, type = "eSNP"))}
#' 
NULL





#' Toy dataset with SNP data
#' 
#' This data set gives an \code{\linkS4class{Enrichment}} object after,
#' \code{\link{initFiles}} and \code{\link{readEnrichment}} is ran. Compute LD
#' for all SNPs in \code{snpListDir} files two by two. Genome Build 37.3
#' (hg19).
#' 
#' 
#' @name toyEnrichment-dataset
#' @aliases toyEnrichment toyEnrichment-dataset
#' @docType data
#' @format See class \code{\linkS4class{Enrichment}} for details about the
#' format.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords datasets toyEnrichment snpEnrichment
NULL





#' Transcript information in order to check the CIS status for SNPs
#' 
#' This dataset is used by \code{\link{readEnrichment}} and
#' \code{\link{compareEnrichment}} in order to check the CIS status for each
#' SNP of signal. Genome Build 37.3 (hg19).
#' 
#' 
#' @name transcript-dataset
#' @aliases transcript transcript-dataset
#' @docType data
#' @format See class \code{\link{readEnrichment}} and
#' \code{\link{compareEnrichment}} for details about how to use this dataset.
#' @author Mickael Canouil \email{mickael.canouil@@good.ibl.fr}
#' @seealso Overview : \code{\link{snpEnrichment-package}} \cr Classes :
#' \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}},
#' \code{\linkS4class{EnrichSNP}} \cr Methods : \code{\link{plot}},
#' \code{\link{reSample}}, \code{\link{getEnrichSNP}},
#' \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
#' \code{\link{enrichment}}, \code{\link{is.enrichment}},
#' \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr Functions :
#' \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
#' @keywords datasets snpEnrichment
NULL



