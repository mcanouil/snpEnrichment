\name{snpEnrichment-package}
\alias{snpEnrichment-package}
\docType{package}
\alias{snpEnrichment}

\title{~ Overview: SNPs enrichment analysis ~}

\description{
Implements classes and methods for large-scale SNP enrichment analysis (e.g. SNPs associated with genes expression in a GWAS signal).
}

\details{
    \tabular{ll}{
    Package: \tab snpEnrichment\cr
    Title: \tab SNPs enrichment analysis\cr
    Author: \tab Mickael Canouil <mickael.canouil at good.ibl.fr>\cr
    Contributor: \tab Loic Yengo\cr
    Maintainer: \tab Mickael Canouil <mickael.canouil at good.ibl.fr>\cr
    % Description: \tab Implements classes and methods for large-scale SNP enrichment analysis \cr (e.g. SNPs associated with genes expression and a GWAS signal) \cr
    License: \tab GPL (>= 2)\cr
    Depends: \tab R (>= 2.15.0), methods\cr
    Suggests: \tab grid, ggplot2\cr
    Imports: \tab parallel, snpStats\cr
    URL: \tab https://github.com/mcanouil/snpEnrichment\cr
    Encoding: \tab UTF-8\cr
    }
}

\note{
Internal data management in 'snpEnrichment' use RefSNP (rs) IDs.
}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{EnrichSNP}} \cr
Methods :  \code{\link{plot}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{snpEnrichment}
\keyword{package}
\keyword{Enrichment}

\examples{
#######################
### 1. Data preparation
\dontrun{
snpInfoDir <- system.file("extdata/snpInfo",
                          package = "snpEnrichment")
signalFile <- system.file("extdata/Signal/toySignal.txt",
                          package = "snpEnrichment")
initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)
writeLD(pattern = "Chrom", snpInfoDir, signalFile,
        ldDir = NULL, ldThresh = 0.8, depth = 1000,
        mc.cores = 1)}

###################
### 2. Reading data
\dontrun{
snpListDir <- system.file("extdata/List",
                          package = "snpEnrichment")
data(transcript)
transcriptFile <- transcript

toy_M1 <- readEnrichment(pattern = "Chrom", signalFile,
                         transcriptFile, snpListDir,
                         snpInfoDir, distThresh = 1000,
                         sigThresh = 0.05, LD = FALSE,
                         ldDir = NULL, mc.cores = 1)
toy_M1}


########################
### 3. Computing results
\dontrun{
reSample(object = toy_M1,
         nSample = 10,
         empiricPvalue = FALSE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1,
         onlyGenome = FALSE)
# OR
data(toyEnrichment)
reSample(object = toyEnrichment,
         nSample = 10,
         empiricPvalue = FALSE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1,
         onlyGenome = FALSE)
}


#######################
### 4. Further analysis: Exclude SNP from original list.
\dontrun{
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


####################
### 5. Watch results
\dontrun{
show(toy_M1)
summary(toy_M1)

show(toyEnrichment)
summary(toyEnrichment)

show(toyEnrichment_exclude)
summary(toyEnrichment_exclude)}
}
