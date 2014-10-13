\name{getEnrichSNP}
\alias{getEnrichSNP}
\docType{methods}
\alias{getEnrichSNP-methods}
\alias{getEnrichSNP,Enrichment-method}
% \alias{getEnrichSNP,Enrichment,ANY-method}
\alias{getEnrichSNP,ANY-method}

\title{
Get all eSNP/xSNP which are enriched
}

\description{
\code{\link{getEnrichSNP}} get all eSNP/xSNP in a \code{\linkS4class{Enrichment}} object which are significant in the signal according to \code{sigThresh} defined in \code{\link{readEnrichment}}.
}

\usage{
\S4method{getEnrichSNP}{Enrichment}(object, type = "eSNP")
}

\arguments{
    \item{object}{[Enrichment]: an object of class \code{\linkS4class{Enrichment}}.}
    \item{type}{[character]: extract \code{eSNP} or \code{xSNP} data.}
}

% \details{}

\value{
Return a \code{data.frame} with eSNP/xSNP which are enriched in signal given to \code{signalFile} in function \code{\link{readEnrichment}}.
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

\keyword{getEnrichSNP}
\keyword{methods}

\examples{
\dontrun{data(toyEnrichment)
eSNPenriched <- getEnrichSNP(object = toyEnrichment, type = "eSNP")
head(eSNPenriched)}
}
