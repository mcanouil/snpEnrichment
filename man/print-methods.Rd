\name{print-methods}
\docType{methods}
\alias{print}
\alias{print-methods}
\alias{print,Chromosome-method}
% \alias{print,Chromosome,ANY,ANY-method}
\alias{print,Enrichment-method}
% \alias{print,Enrichment,ANY,ANY-method}

\title{
Print method (S4)
}

\description{
\code{\link{print}} is a generic function used to print results.
}

\usage{
\S4method{print}{Enrichment}(x, what = "Genome", type = c("eSNP", "xSNP"))

\S4method{print}{Chromosome}(x, type = c("eSNP", "xSNP"))
}

\arguments{
    \item{x}{[Enrichment or Chromosome]: an object of class \code{\linkS4class{Enrichment}} or \code{\linkS4class{Chromosome}}.}
    \item{what}{[character or numeric]: \code{what="Genome"} (default) to print results as a matrix. \code{what} could be \code{"All"}, \code{"Genome"} or a numeric from 1 to 22 (numeric vector is allowed).}
    \item{type}{[character]: select if results for \code{"eSNP"} and/or \code{"xSNP"} should be print.}
}

% \details{}

\value{
Return a \code{matrix} for classes \code{\linkS4class{Enrichment}} and  \code{\linkS4class{Chromosome}}.\cr
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

\keyword{print}
\keyword{methods}

\examples{
data(toyEnrichment)
print(toyEnrichment, "All", type = "eSNP")
print(toyEnrichment, "Genome")
print(toyEnrichment, 1)
}
