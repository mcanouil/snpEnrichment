\name{transcript-dataset}
\docType{data}
\alias{transcript}
\alias{transcript-dataset}

\title{Transcript information in order to check the CIS status for SNPs}

\description{
This dataset is used by \code{\link{readEnrichment}} and \code{\link{compareEnrichment}}
in order to check the CIS status for each SNP of signal.
Genome Build 37.3 (hg19).
}

\usage{
data(transcript)
transcript
}

\format{
See class \code{\link{readEnrichment}} and \code{\link{compareEnrichment}} for details about how to use this dataset.
}

% \source{}

% \references{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{EnrichSNP}} \cr
Methods : \code{\link{plot}}, \code{\link{reSample}}, \code{\link{getEnrichSNP}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{datasets}
\keyword{snpEnrichment}