\name{Chromosome-class}
\docType{class}
\alias{Chromosome-class}
\alias{Chromosome}
\alias{[,Chromosome-method}
\alias{[,Chromosome,ANY,ANY,ANY-method}
\alias{[<-,Chromosome-method}
\alias{[<-,Chromosome,ANY,ANY,ANY-method}
% \alias{computeER,Chromosome-method}
% \alias{doLDblock,Chromosome-method}
% \alias{reSample,Chromosome-method}
% \alias{reset,Chromosome-method}
\alias{show,Chromosome-method}
% \alias{print,Chromosome-method}
\alias{chromosome}
\alias{chromosome-methods}
% \alias{chromosome,ANY,ANY,ANY,ANY-method}
\alias{chromosome,ANY-method}

\title{Class \code{\linkS4class{Chromosome}}}

\description{
This class is defined to summarize the enrichment analysis about a chromosome.
}

\section{Objects from the Class}{
    \code{\link{chromosome}} is defined to build an object of class
    \code{\linkS4class{Chromosome}} in order to compute an enrichment analysis.
    A \code{\linkS4class{Chromosome}} object contains the original data,
    a list of SNPs, some results and resampling data.\cr

    When created, an \code{\linkS4class{Chromosome}} object is "empty".
    \code{\link{readEnrichment}} initializes a \code{\linkS4class{Chromosome}} object with value from PLINK computation and user's files.
    In this step, only the fields "Data", "LD", "SNP" are filled.
    \code{reSample} fills the fields: Table, EnrichmentRatio, Z, PValue and Resampling of a \code{\linkS4class{Chromosome}}.

    Note that if \code{\link{reSample}} is executed on an \code{\linkS4class{Chromosome}} every new resampling
    is added to the original ones, pre-existing statistics are erased and computed again with the new resampling set.
}

\section{Slots}{
    \describe{
        \item{Data}{[data.frame]: a \code{data.frame} with 6 columns ("SNP", "PVALUE", "CHR", "MAF", "eSNP", "xSNP"). Where "eSNP" and "xSNP" are logical columns defining the lists of SNPs (extended or not).}
        \item{LD}{[data.frame]: a \code{data.frame} which contains LD informations between SNPs (computed with \code{\link{writeLD}} or PLINK).}
        \item{eSNP, xSNP}{[SNP]: contain a \code{\linkS4class{EnrichSNP}} object (whith slots: SNP, Table, EnrichmentRatio, Z, PValue and Resampling) for a list of SNPs (\code{eSNP}) and an extended one (\code{xSNP}).}
        % \item{xSNP}{[SNP]: contains a \code{\linkS4class{EnrichSNP}} object (whith slots: SNP, Table, EnrichmentRatio, Z, PValue and Resampling) for an extended list of SNPs.}
    }
}

\section{Extends}{
Class \code{\linkS4class{EnrichSNP}}, directly.\cr
}

\section{Methods}{
    \describe{
        \item{chromosome(Data, LD, eSNP, xSNP):}{Generate and initialize a new \code{\linkS4class{Chromosome}} object.}
        \item{object["slotName"]:}{Get the value of the field \code{slotName}.}
        \item{object["slotName"]<-value:}{Set \code{value} to the field \code{slotName}.}
        \item{show(object):}{Return the formatted values of \code{\linkS4class{Chromosome}} object.}
        % \item{print(object):}{Return the formatted values of \code{\linkS4class{Chromosome}} object.}
        % \item{computeER(object, sigThresh, mc.cores):}{Compute an Enrichment Ratio for \code{\linkS4class{Chromosome}} object (internal use only).}
        % \item{doLDblock(object, mc.cores):}{Compute LD block and update \code{Data} field (internal use only).}
        % \item{reSample(object, nSample, sigThresh, MAFpool, extendMethod, mc.cores):}{Compute P-Value based upon a resampling of SNPs (\code{eSNP} and \code{xSNP}) and update an \code{\linkS4class{Chromosome}} object.}
        % \item{excludeSNP(object, excludeFile, nSample, MAFpool, mc.cores):}{Exclude SNPs from eSNPs list and xSNPs list and compute a new enrichment analysis.}
        % \item{reset(object, "slotName"):}{Reset the field \code{slotName}.}
    }
}

\note{
\code{\linkS4class{Chromosome}} object is not intended to be used alone on this version.\\
It is a part of the \code{\linkS4class{Enrichment}} object.
}


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
\keyword{chromosome}

\examples{
Data <- structure(
    list(
        SNP = c("rs4970420", "rs3766178",
                "rs10910030", "rs10910036",
                "rs2234167", "rs6661861"),
        PVALUE = c(0.9244, 0.167, 0.01177, 0.4267, 0.9728, 0.4063),
        CHR = c(1, 1, 1, 1, 1, 1),
        POS = c(1106473, 1478180, 2035684, 2183754, 2494330, 3043121),
        MAF = c(0.2149, 0.3117, 0.374, 0.3753, 0.1565, 0.06101),
        eSNP = c(0, 1, 1, 0, 0, 0),
        xSNP = c(0, 1, 1, 0, 0, 0)
    ),
    .Names = c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP", "xSNP"),
    row.names = c("rs4970420", "rs3766178",
                  "rs10910030", "rs10910036",
                  "rs2234167", "rs6661861"),
class = "data.frame")

toyChr <- chromosome(Data = Data)
show(toyChr)
toyChr

toyChr <- chromosome()
toyChr["Data"] <- Data
toyChr

is.chromosome(toyChr) # TRUE
}
