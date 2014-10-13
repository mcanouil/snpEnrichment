\name{writeLD}
\alias{writeLD}

\title{
Linkage Disequilibrium (LD) computation with PLINK
}

\description{
\code{\link{writeLD}} write a '.ld' file for each chromosomes which contains the LD (r^2).
}

\usage{
writeLD(pattern = "Chrom", snpInfoDir, signalFile,
        ldDir = NULL, ldThresh = 0.8, depth = 1000,
        mc.cores = 1)
}

\arguments{
    \item{pattern}{[character]: character string containing a expression to be matched with all chromosomes files (e.g."Chrom" for  files which start by "Chrom" followed by the chromosome number).}
    \item{snpInfoDir}{[character]: character string naming a directory containing the reference data in a PLINK format (*.bed, *.bim and *.fam).}
    % \item{snpListDir}{[character]: character string naming a directory containing a list of SNPs for each chromosomes.}
    \item{signalFile}{[character]: the name of the signal file which the data are to be read from (2 columns: "SNP" and "PVALUE").
    Each row of the table appears as one line of the file.
    If it does not contain an \code{_absolute_} path, the file name is \code{_relative_} to the current working directory, \code{getwd}.
    The fields separator character have to be a space \code{" "} or a tabulation \code{"\\t"}.}
    \item{ldDir}{[character]: character string naming a directory where the linkage Disequilibrium files should be written (default \code{ldDir=NULL} is in temporary directory).}
    \item{ldThresh}{[numeric]: threshold  value for LD calculation.}
    \item{depth}{[numeric]: this parameter is mandatory and controls the maximum lag between SNPs considered.}
    \item{mc.cores}{[numeric]: the number of cores to use (default is \code{mc.cores=1}), i.e. at most how many child processes will be run simultaneously.
    Must be at least one, and parallelization requires at least two cores.}
}

% \details{}

\value{
One ".ld" file per chromosome is returned by \code{\link{writeLD}} in \code{snpInfoDir} directory.
}

\note{
The LD computation can take a long time depending on number of SNPs in \code{signalFile}. It is recommended to save LD results in a directory (\code{ldDir}) which is not a temporary directory.
}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{EnrichSNP}} \cr
Methods : \code{\link{plot}}, \code{\link{reSample}}, \code{\link{getEnrichSNP}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{initFiles}
\keyword{Enrichment}

\examples{
\dontrun{signalFile <- system.file("extdata/Signal/toySignal.txt",
                          package = "snpEnrichment")
snpInfoDir <- system.file("extdata/snpInfo",
                          package = "snpEnrichment")
writeLD(pattern = "Chrom", snpInfoDir, signalFile,
        ldDir = NULL, ldThresh = 0.8, mc.cores = 1)}
}
