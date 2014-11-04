\name{EnrichSNP-class}
\docType{class}
\alias{EnrichSNP-class}
\alias{EnrichSNP}
\alias{[,EnrichSNP-method}
\alias{[,EnrichSNP,ANY,ANY,ANY-method}
\alias{[<-,EnrichSNP-method}
\alias{[<-,EnrichSNP,ANY,ANY,ANY-method}
\alias{show,EnrichSNP-method}
\alias{print,EnrichSNP-method}
% \alias{EnrichSNP}
% \alias{EnrichSNP-methods}
% \alias{EnrichSNP,ANY,ANY,ANY,ANY,ANY,ANY-method}

\title{Class \code{"\linkS4class{EnrichSNP}"}}

\description{
This class is defined to summarize the enrichment analysis.
It's a part of \code{\linkS4class{Chromosome}} and \code{\linkS4class{Enrichment}} classes.
}

% \section{Objects from the Class}{}

\section{Slots}{
    \describe{
        \item{List}{[vector(character)]: a list of SNPs used to compute enrichment (e.g. eSNP or xSNP).}
        \item{Table}{[matrix]: Contingency table with SNPs (columns) and P-Values from signal (rows).}
        \item{EnrichmentRatio}{[numeric]: Enrichment Ratio is computed on the contingency table (\code{Table} slot).}
        \item{Z}{[numeric]: A statistic computed from \code{EnrichmentRatio} and resampling results.}
        \item{PValue}{[numeric]: P-Value associated with the statistic \code{Z}.}
        \item{Resampling}{[matrix]: A matrix with by row, the contingency table and the odds ratio for each resampling.}
    }
}

\section{Methods}{
    \describe{
        % \item{enrichSNP(List, Table, EnrichmentRatio, Z, PValue, Resampling): }{Generate and initialize a new \code{\linkS4class{EnrichSNP}} object.}
        \item{object["slotName"]: }{Get the value of the field \code{slotName}.}
        \item{object["slotName"]<-value: }{Set \code{value} to the field \code{slotName}.}
        \item{show(object): }{Return the formatted values of \code{\linkS4class{EnrichSNP}} object.}
        % \item{print(x): }{Return the values of \code{\linkS4class{EnrichSNP}} object as a matrix.}
    }
}

\note{
\code{\linkS4class{EnrichSNP}} object is not intended to be use directly by user.
It is a part of the \code{\linkS4class{Enrichment}} and \code{\linkS4class{Chromosome}} object.
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
\keyword{enrichSNP}

% \examples{}
