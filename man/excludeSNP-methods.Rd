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
}

% \details{}

\value{
Return the object given in argument where lists of SNPs are updated by removing SNPs in \code{excludeFile}.
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

\keyword{excludeSNP}
\keyword{methods}

\examples{
\dontrun{data(toyEnrichment)
excludeFile <- c(
    "rs4376885", "rs17690928", "rs6460708", "rs13061537", "rs11769827",
    "rs12717054", "rs2907627", "rs1380109", "rs7024214", "rs7711972",
    "rs9658282", "rs11750720", "rs1793268", "rs774568", "rs6921786",
    "rs1699031", "rs6994771", "rs16926670", "rs465612", "rs3012084",
    "rs354850", "rs12803455", "rs13384873", "rs4364668", "rs8181047",
    "rs2179993", "rs12049335", "rs6079926", "rs2175144", "rs11564427",
    "rs7786389", "rs7005565", "rs17423335", "rs12474102", "rs191314",
    "rs10513168", "rs1711437", "rs1992620", "rs283115", "rs10754563",
    "rs10851727", "rs2173191", "rs7661353", "rs1342113", "rs7042073",
    "rs1567445", "rs10120375", "rs550060", "rs3761218", "rs4512977"
)
# OR
excludeFile <- system.file("extdata/Exclude/toyExclude.txt",
                           package = "snpEnrichment")

toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)
toyEnrichment_exclude}
}
