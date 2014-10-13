\name{is.chromosome}
\alias{is.chromosome}
% \docType{methods}
\alias{is.chromosome-methods}
\alias{is.chromosome,ANY-method}

\title{
Is an Chromosome object
}

\description{
'is.chromosome' returns 'TRUE' if 'x' is an \code{\linkS4class{Chromosome}} object and 'FALSE' otherwise.
}

\usage{
is.chromosome(object)
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

\keyword{chromosome}
\keyword{is}
\keyword{is.chromosome}

\examples{
a <- chromosome()
c <- chromosome()
is.chromosome(list())                # FALSE
is.chromosome(1)                     # FALSE
is.chromosome(a)                     # TRUE
is.chromosome(c(a, c))               # TRUE TRUE
is.chromosome(list(a, b = "char"))   # TRUE FALSE
is.chromosome(c(a, b = list(12, c))) # TRUE FALSE TRUE
}
