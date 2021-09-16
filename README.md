
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SNPs Enrichment Analysis With Sampling <!--<img src="man/figures/snpEnrichment_hex.png" align="right" width="120" /> -->

<!-- badges: start -->

[![Lifecycle:
retired](https://img.shields.io/badge/lifecycle-retired-orange.svg)](https://www.tidyverse.org/lifecycle/#retired)
[![GitHub
tag](https://img.shields.io/github/tag/mcanouil/snpEnrichment.svg?label=latest%20tag)](https://github.com/mcanouil/snpEnrichment)
[![Travis-CI Build
Status](https://travis-ci.org/mcanouil/snpEnrichment.svg?branch=master)](https://travis-ci.org/mcanouil/snpEnrichment)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/mcanouil/snpEnrichment?branch=master&svg=true)](https://ci.appveyor.com/project/mcanouil/snpEnrichment)
[![Coverage Status
(codecov)](https://codecov.io/gh/mcanouil/snpEnrichment/branch/master/graph/badge.svg)](https://codecov.io/gh/mcanouil/snpEnrichment)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version-ago/snpEnrichment)](https://cran.r-project.org/package=snpEnrichment)
[![cran
checks\_worst](https://cranchecks.info/badges/worst/snpEnrichment)](https://cran.r-project.org/web/checks/check_results_snpEnrichment.html)
[![CRAN\_Download\_total](http://cranlogs.r-pkg.org/badges/grand-total/snpEnrichment)](https://cran.r-project.org/package=snpEnrichment)
<!--[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/2719/badge)](https://bestpractices.coreinfrastructure.org/projects/2719) -->
<!--[![cran checks_summary](https://cranchecks.info/badges/summary/snpEnrichment)](https://cran.r-project.org/web/checks/check_results_snpEnrichment.html) -->
<!--[![CRAN_Download_month](http://cranlogs.r-pkg.org/badges/snpEnrichment?color=brightgreen)](https://cran.r-project.org/package=snpEnrichment) -->
<!--[![Coverage Status (coveralls)](https://coveralls.io/repos/github/mcanouil/snpEnrichment/badge.svg?branch=master)](https://coveralls.io/github/mcanouil/snpEnrichment?branch=master) -->
<!-- badges: end -->

Implements classes and methods for large scale SNP enrichment analysis
(*e.g.*, SNPs associated with genes expression in a GWAS signal).

**Not maintained anymore.**

## Installation

``` r
# Install snpEnrichment from CRAN:
install.packages("snpEnrichment")

# Or the the development version from GitHub:
# install.packages("remotes")
remotes::install_github("mcanouil/snpEnrichment")
```

## Overview

**Note**: Internal data management in ‘snpEnrichment’ use RefSNP (rs)
IDs.

### Load snpEnrichment

``` r
library(snpEnrichment)
```

### Prepare data

``` r
snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment")
signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment")

initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)

writeLD(
  pattern = "Chrom", 
  snpInfoDir, 
  signalFile, 
  ldDir = NULL,
  ldThresh = 0.8, 
  depth = 1000, 
  mc.cores = 1
)
```

### Read data

``` r
snpListDir <- system.file("extdata/List", package = "snpEnrichment")
data(transcript)
transcriptFile <- transcript

toyData <- readEnrichment(
    pattern = "Chrom", 
    signalFile, 
    transcriptFile,
    snpListDir, 
    snpInfoDir,
    distThresh = 1000, 
    sigThresh = 0.05, 
    LD = TRUE, 
    ldDir = NULL,
    mc.cores = 1
)
toyData
```

### Compute results

``` r
reSample(
  object = toyData,
  nSample = 10,
  empiricPvalue = TRUE,
  MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
  mc.cores = 1,
  onlyGenome = TRUE
)
```

### 5\. Further analysis: Exclude SNP from original list.

``` r
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

# or
excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment")
```

``` r
toyData_exclude <- excludeSNP(toyData, excludeFile, mc.cores = 1)
```

**Warning:** `compareEnrichment` is in development\!\!

``` r
compareResults <- compareEnrichment(
  object.x = toyData,
  object.y = toyData_exclude,
  pattern = "Chrom",
  nSample = 10,
  empiricPvalue = TRUE,
  mc.cores = 1,
  onlyGenome = TRUE
)
```

### Watch results

``` r
show(toyData)
print(toyData)
head(getEnrichSNP(toyData, type = "xSNP"))

show(toyData_exclude)
print(toyData_exclude)
head(getEnrichSNP(toyData_exclude, type = "eSNP"))
```

## Getting help

If you encounter a clear bug, please file a minimal reproducible example
on [github](https://github.com/mcanouil/CARoT/issues). For questions and
other discussion, please contact the package maintainer.

-----

Please note that this project is released with a [Contributor Code of
Conduct](.github/CODE_OF_CONDUCT.md). By participating in this project
you agree to abide by its terms.
