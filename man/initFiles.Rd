\name{initFiles}
\alias{initFiles}

\title{
Initialize files for enrichment analysis
}

\description{
\code{\link{initFiles}} create several files needed to run \code{\link{readEnrichment}}.
".frq" and ".signal" are created with PLINK.
LD computation can be run with \code{\link{writeLD} or with PLINK}.
}

\usage{
initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)
}

\arguments{
    \item{pattern}{[character]: character string containing a expression to be matched with all chromosomes files (e.g."Chrom" for  files which start by "Chrom" followed by the chromosome number).}
    \item{snpInfoDir}{[character]: character string naming a directory containing the reference data in a PLINK format (*.bed, *.bim and *.fam).}
    \item{signalFile}{[character]: the name of the signal file which the data are to be read from (2 columns: "SNP" and "PVALUE").
    Each row of the table appears as one line of the file.
    If it does not contain an \code{_absolute_} path, the file name is \code{_relative_} to the current working directory, \code{getwd}.
    The fields separator character have to be a space \code{" "} or a tabulation \code{"\\t"}.}
    \item{mc.cores}{[numeric]: the number of cores to use (default is \code{mc.cores=1}), i.e. at most how many child processes will be run simultaneously.
    Must be at least one, and parallelization requires at least two cores.}
}

% \details{}

\value{
This function writes several files, in the temporary directory (defined in \code{R_SESSION_TMPDIR}), nothing else is returned.
These files are used to build an \code{\linkS4class{Enrichment}} object by \code{\link{readEnrichment}} in order to compute enrichment analysis (\code{\link{reSample}}).
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

\keyword{initFiles}
\keyword{initialize}

\examples{
\dontrun{snpInfoDir <- system.file("extdata/snpInfo",
                          package = "snpEnrichment")
signalFile <- system.file("extdata/Signal/toySignal.txt",
                          package = "snpEnrichment")
initFiles(pattern = "Chrom",
          snpInfoDir,
          signalFile,
          mc.cores = 1)}
}
