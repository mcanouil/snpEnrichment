\name{compareEnrichment}
\alias{compareEnrichment}
\docType{methods}
\alias{compareEnrichment-methods}
% \alias{compareEnrichment,ANY,ANY,ANY,ANY,ANY,ANY-method}
\alias{compareEnrichment,ANY-method}
\alias{compareEnrichment,Enrichment,Enrichment,ANY-method}

\title{
Compare enrichment analysis between two SNPs list
}

\description{
Compare the enrichment analysis between two set of SNPs.
\code{\link{compareEnrichment}} compare two \code{\linkS4class{Enrichment}} objects.
}

\usage{
compareEnrichment(object.x, object.y, pattern = "Chrom",
                  nSample = 100, empiricPvalue = FALSE,
                  mc.cores = 1, onlyGenome = FALSE)
}

\arguments{
    \item{object.x, object.y}{[Enrichment]: an \code{\linkS4class{Enrichment}} object fully filled (e.g. \code{\link{readEnrichment}}).}
    \item{pattern}{[character]: character string containing a expression to be matched with all chromosomes files (e.g."Chrom" for  files which start by "Chrom" followed by the chromosome number).}
    \item{nSample}{[numeric]: the number of resampling done by \code{\link{reSample}} for p-values computation (minimum is 100).}
    \item{empiricPvalue}{[logical]: \code{empiricPvalue=FALSE} (default) compute PValue based on the null distribution (resampling).
    If \code{empiricPvalue=TRUE}, the empirical p-values are computed instead.}
    \item{mc.cores}{[numeric]: the number of cores to use (default is \code{mc.cores=1}), i.e. at most how many child processes will be run simultaneously.
    Must be at least one, and parallelization requires at least two cores.}
    \item{onlyGenome}{[logical]: \code{onlyGenome=FALSE} (default) compute resampling step for all chromosomes.}
}

% \details{}

\value{
Return a \code{list} of three elements:
    \item{summary}{Two \code{data.frame} which summarize for \code{eSNP} and/or \code{xSNP} the enrichment analyses and the comparison between them.}
    \item{object1}{\code{\linkS4class{Enrichment}} object passed in \code{object.x} with resampling data.}
    \item{object2}{\code{\linkS4class{Enrichment}} object passed in \code{object.y} with resampling data.}
    \item{comparison}{\code{\linkS4class{Enrichment}} object from the comparison between \code{object.x} and \code{object.y}.}
}

\note{
Still in development.
}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{EnrichSNP}} \cr
Methods :  \code{\link{plot}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{compareEnrichment}
\keyword{Enrichment}

\examples{
\dontrun{
data(toyEnrichment)

reSample(object = toyEnrichment,
         nSample = 10,
         empiricPvalue = FALSE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1,
         onlyGenome = TRUE)

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
# OR
excludeFile <- system.file("extdata/Exclude/toyExclude.txt",
                           package = "snpEnrichment")

toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)

compareResults <- compareEnrichment(object.x = toyEnrichment,
                                    object.y = toyEnrichment_exclude,
                                    pattern = "Chrom",
                                    nSample = 10,
                                    empiricPvalue = FALSE,
                                    mc.cores = 1,
                                    onlyGenome = TRUE)}
}
