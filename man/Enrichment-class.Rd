\name{Enrichment-class}
\docType{class}
\alias{Enrichment-class}
\alias{Enrichment}
\alias{[,Enrichment-method}
\alias{[,Enrichment,ANY,ANY,ANY-method}
\alias{[<-,Enrichment-method}
\alias{[<-,Enrichment,ANY,ANY,ANY-method}
% \alias{computeER,Enrichment-method}
% \alias{excludeSNP,Enrichment-method}
% \alias{reSample,Enrichment-method}
% \alias{reset,Enrichment-method}
\alias{show,Enrichment-method}
% \alias{print,Enrichment-method}
\alias{enrichment}
\alias{enrichment-methods}
% \alias{enrichment,ANY,ANY,ANY,ANY,ANY-method}
\alias{enrichment,ANY-method}

\title{Class \code{\linkS4class{Enrichment}}}

\description{
This class is defined to summarize the enrichment analysis on each chromosomes and the whole genome.
}

\section{Objects from the Class}{
    \code{\link{enrichment}} is defined to build an object of class \code{\linkS4class{Enrichment}} in order to compute an enrichment analysis.
    \code{Enrichment} is the object containing the results for all \code{\linkS4class{Chromosome}} object and for the whole genome.\cr

    When an \code{\linkS4class{Enrichment}} object is created, it contains a list of SNPs (e.g. eSNPs).
    All the others slots are "empty". After \code{\link{reSample}} is ran on an \code{\linkS4class{Enrichment}} object,
    the slots: Table, EnrichmentRatio, Z, PValue and Resampling are filled.

    Note that if \code{\link{reSample}} is executed on an \code{\linkS4class{Enrichment}} every new resampling
    is added to the original ones, pre-existing statistics are erased and computed again with the new resampling set.
}

\section{Slots}{
    \describe{
        % \item{Signal}{[data.frame]: a three columns \code{data.frame}: "SNP", "PVALUE" and "IN" (e.g. GWAS). "IN" is computed during the reading step and gives informations about which SNPs are kept for the enrichment analysis.}
        \item{Loss}{[data.frame]: a four columns \code{data.frame}: "Rows", "Unique", "Intersect.Ref.Signal" and "CIS". This slot gives information on data losses.}
        \item{Call}{[list]: each parameters used for the reading or resampling step are stored in this slot.}
        \item{eSNP, xSNP}{[SNP]: contain a \code{\linkS4class{EnrichSNP}} object (whith slots: SNP, Table, EnrichmentRatio, Z, PValue and Resampling) for a list of SNPs (\code{eSNP}) and an extended one (\code{xSNP}).}
        % \item{xSNP}{[SNP]: contains a \code{\linkS4class{EnrichSNP}} object (whith slots: SNP, Table, EnrichmentRatio, Z, PValue and Resampling) for an extended list of SNPs.}
        \item{Chromosomes}{[list(Chromosome)]: a list of 22 \code{\linkS4class{Chromosome}} objects.}
    }
}

\section{Extends}{
Class \code{\linkS4class{Chromosome}}, directly.\cr
Class \code{\linkS4class{EnrichSNP}}, directly.\cr
}

\section{Methods}{
    \describe{
        \item{enrichment(Loss, Call, eSNP, xSNP, Chromosomes):}{Generate and initialize a new \code{\linkS4class{Enrichment}} object.}
        \item{object["slotName"]:}{Get the value of the field \code{slotName}.}
        \item{object["slotName"]<-value:}{Set \code{value} to the field \code{slotName}.}
        \item{show(object):}{Return the formatted values of \code{\linkS4class{Enrichment}} object.}
        \item{print(object, what, type):}{Return a summary table of an \code{\linkS4class{Enrichment}} object.}
        % \item{computeER(object, sigThresh):}{Compute Enrichment Ratio for \code{\linkS4class{Enrichment}} object (internal use only).}
        \item{reSample(object, nSample, MAFpool, mc.cores, onlyGenome):}{Compute P-Values based upon a resampling of SNPs (\code{eSNP} and \code{xSNP}) and update the \code{\linkS4class{Enrichment}} object.}
        \item{excludeSNP(object, excludeFile, mc.cores):}{Excludes SNPs given in \code{excludeFile} from the original list of eSNPs (xSNPs). Then a new enrichment analysis is computed.}
        \item{reset(object, "slotName"):}{Reset the field \code{slotName}.}
        \item{plot(object, what, type, ggplot, pvalue):}{Plot p-values or Z-statistics convergence.}
        \item{getEnrichSNP(object, type):}{Return eSNP/xSNP which are enriched as a \code{data.frame}.}
    }
}

% \note{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{EnrichSNP}} \cr
Methods : \code{\link{plot}}, \code{\link{reSample}}, \code{\link{getEnrichSNP}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{classes}
\keyword{class}
\keyword{enrichment}

\examples{
data(toyEnrichment)
toyEnrich <- enrichment()
show(toyEnrich)

toyEnrich["Loss"] <- toyEnrichment["Loss"]
toyEnrich["Loss"]

toyEnrich <- enrichment(Loss = toyEnrichment["Loss"],
                        eSNP = toyEnrichment["eSNP"])
toyEnrich <- enrichment(Loss = toyEnrichment["Loss"])

\dontrun{reSample(object = toyEnrichment,
         nSample = 10,
         empiricPvalue = TRUE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1,
         onlyGenome = TRUE)
print(toyEnrichment)

excludeFile <- c(
"rs7897180", "rs4725479", "rs315404", "rs17390391", "rs1650670",
    "rs6783390", "rs1642009", "rs4756586", "rs11995037", "rs4977345",
    "rs13136448", "rs4233536", "rs11151079", "rs2299657", "rs4833930",
    "rs1384", "rs7168184", "rs6909895", "rs7972667", "rs2293229",
    "rs918216", "rs6040608", "rs2817715", "rs13233541", "rs4486743",
    "rs2127806", "rs10912854", "rs1869052", "rs9853549", "rs448658",
    "rs2451583", "rs17483288", "rs10962314", "rs9612059", "rs1384182",
    "rs8049208", "rs12215176", "rs2980996", "rs1736976", "rs8089268",
    "rs10832329", "rs12446540", "rs7676237", "rs869922", "rs16823426",
    "rs1374393", "rs13268781", "rs11134505", "rs7325241", "rs7520109"
)

toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)
print(toyEnrichment_exclude)}
}
