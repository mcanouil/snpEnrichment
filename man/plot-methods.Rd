\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,Enrichment-method}
\alias{plot,Enrichment,ANY-method}

\title{
Plot method (S4) for \code{\linkS4class{Enrichment}} object
}

\description{
\code{\link{plot}} is a generic function for plotting of R objects.
The function invokes particular \code{methods} which depend on the \code{class} of the first argument.
}

\usage{
\S4method{plot}{Enrichment}(x, what = "Genome", type = c("eSNP", "xSNP"),
           ggplot = FALSE, pvalue = TRUE, ...)
}

\arguments{
    \item{x}{[Enrichment]: an object of class \code{\linkS4class{Enrichment}} which the Z statistics or p-values have to be drawn.}
    \item{what}{[character or vector(numeric)]: default \code{what="Genome"}) plot Z statistics or p-values for genome only (what must be: \code{"All"}, \code{"Genome"} or numeric vector).}
    \item{type}{[vector(character)]: plot the selected analysis for \code{"eSNP"} and/or \code{"xSNP"}.}
    \item{ggplot}{[logical]: use ggplot (default \code{ggplot=FALSE}) instead of classic plot method.}
    \item{pvalue}{[logical]: if \code{TRUE}, p-value convergense is plotted. Otherwise, Z statistic is plotted.}
    \item{...}{[any]: Arguments to be passed to methods, such as graphical parameters (see \code{par})}
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

\keyword{plot}
\keyword{methods}

\examples{
\dontrun{data(toyEnrichment)
reSample(toyEnrichment, 10)
plot(toyEnrichment)}
}
