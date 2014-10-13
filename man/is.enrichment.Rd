\name{is.enrichment}
\alias{is.enrichment}
% \docType{methods}
\alias{is.enrichment-methods}
\alias{is.enrichment,ANY-method}

\title{
Is an Enrichment object
}

\description{
'is.enrichment' returns 'TRUE' if 'x' is an \code{\linkS4class{Enrichment}} object and 'FALSE' otherwise.
}

\usage{
is.enrichment(object)
}

\arguments{
    \item{object}{[ANY]: object to be tested.}
}

% \details{}

% \value{}

% \note{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{EnrichSNP}} \cr
Methods : \code{\link{plot}}, \code{\link{reSample}}, \code{\link{getEnrichSNP}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{enrichment}
\keyword{is}
\keyword{is.enrichment}

\examples{
a <- enrichment()
c <- enrichment()
is.enrichment(list())                # FALSE
is.enrichment(1)                     # FALSE
is.enrichment(a)                     # TRUE
is.enrichment(c(a, c))               # TRUE TRUE
is.enrichment(list(a, b = "char"))   # TRUE FALSE
is.enrichment(c(a, b = list(12, c))) # TRUE FALSE TRUE
}
