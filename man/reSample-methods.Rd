\name{reSample}
\alias{reSample}
\docType{methods}
\alias{reSample-methods}
\alias{reSample,Chromosome-method}
% \alias{reSample,Chromosome,ANY,ANY,ANY,ANY,ANY-method}
% \alias{reSample,Chromosome,ANY-method}
\alias{reSample,Enrichment-method}
% \alias{reSample,Enrichment,ANY,ANY,ANY,ANY,ANY-method}
% \alias{reSample,Enrichment,ANY-method}
\alias{reSample,ANY-method}

\title{
Compute enrichment analysis on an \code{\linkS4class{Enrichment}} object
}

\description{
After \code{\link{initFiles}} and \code{\link{readEnrichment}} has been run.
\code{\link{reSample}} computes a statistic value and a p-value for each chromosomes and for the whole genome.
}

\usage{
\S4method{reSample}{Enrichment}(object, nSample = 100,
         empiricPvalue = TRUE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1, onlyGenome = TRUE)

\S4method{reSample}{Chromosome}(object, nSample = 100,
         empiricPvalue = TRUE, sigThresh = 0.05,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1)
}

\arguments{
    \item{object}{[Enrichment or Chromosome]: an object to be updated. It is intended, an object returned by the \code{\link{readEnrichment}} function.}
    \item{nSample}{[numeric]: the number of resampling done by \code{\link{reSample}} for p-values computation (minimum is 100).}
    \item{empiricPvalue}{[logical]: \code{empiricPvalue=TRUE} (default) compute PValue based on the null distribution (resampling).
    If \code{empiricPvalue=TRUE}, the empirical p-values are computed instead.}
    \item{sigThresh}{[numeric]: statistical threshold for signal (e.g. \code{sigThresh = 0.05} for a given GWAS signal) used to compute an Enrichment Ratio.}
    \item{MAFpool}{[vector(numeric)]: either a numeric vector giving the breaks points of intervals into which SNP's MAF (Minor Allele Frequency) is to be split.}
    \item{mc.cores}{[numeric]: the number of cores to use (default is \code{mc.cores=1}), i.e. at most how many child processes will be run simultaneously.
    Must be at least one, and parallelization requires at least two cores.}
    \item{onlyGenome}{[logical]: \code{onlyGenome=TRUE} (default) compute resampling step for all chromosomes.}
}

% \details{}

\value{
Return the object given in argument, updated by the resampling results.
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

\keyword{reSample}
\keyword{methods}

\examples{
\dontrun{data(toyEnrichment)
reSample(object = toyEnrichment,
         nSample = 10,
         empiricPvalue = TRUE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         onlyGenome = TRUE)
toyEnrichment}
}
