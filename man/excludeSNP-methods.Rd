\name{excludeSNP}
\alias{excludeSNP}
\docType{methods}
\alias{excludeSNP-methods}
\alias{excludeSNP,Enrichment-method}
\alias{excludeSNP,ANY-method}
% \alias{excludeSNP,Enrichment,ANY-method}
% \alias{excludeSNP,ANY,ANY,ANY,ANY,ANY,ANY,ANY-method}

\title{
Exclude SNPs from Enrichment analysis
}

\description{
Remove a specify set of SNPs and compute a new enrichment analysis.
}

\usage{
excludeSNP(object, excludeFile, mc.cores = 1)
}

\arguments{
    \item{object}{[Enrichment]: an \code{\linkS4class{Enrichment}} object filled by \code{\link{reSample}}.}
    \item{excludeFile}{[vector(character)]: a list of SNPs to remove from a previous enrichment analysis. A path to a file which the first column are the SNPs.}
    \item{mc.cores}{[numeric]: the number of cores to use (default is \code{mc.cores=1}), i.e. at most how many child processes will be run simultaneously.
    Must be at least one, and parallelization requires at least two cores.}
    % \item{onlyGenome}{[logical]: \code{onlyGenome=FALSE} (default) compute resampling step for all chromosomes.}
}

% \details{}

\value{
Return the object given in argument where lists of SNPs are updated by removing SNPs in \code{excludeFile}.
}

% \note{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{EnrichSNP}} \cr
Methods :  \code{\link{plot}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{excludeSNP}
\keyword{methods}

\examples{
\dontrun{
data(toyEnrichment)
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
toyEnrichment_exclude}
}
