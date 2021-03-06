% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getEnrichSNP.R
\name{getEnrichSNP}
\alias{getEnrichSNP}
\alias{getEnrichSNP,ANY-method}
\alias{getEnrichSNP,Enrichment-method}
\title{Get all eSNP/xSNP which are enriched}
\arguments{
\item{object}{An object of class \linkS4class{Enrichment}.}

\item{type}{A character definined the type of data to extract Extract, \emph{i.e.},  \code{"eSNP"} or \code{"xSNP"}.}
}
\value{
Return a \code{data.frame} with eSNP/xSNP which are enriched in signal
given to \code{signalFile} in function \link{readEnrichment}.
}
\description{
\link{getEnrichSNP} get all eSNP/xSNP in a \linkS4class{Enrichment} object which are significant in the signal
according to \code{sigThresh} defined in \link{readEnrichment}.
}
\examples{
if (interactive()) {
  data(toyEnrichment)
  eSNPenriched <- getEnrichSNP(object = toyEnrichment, type = "eSNP")
  head(eSNPenriched)
}

}
