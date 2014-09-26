snpEnrichment
=============

Implements classes and methods for large scale SNP enrichment analysis.
(e.g. SNPs associated with genes expression in a GWAS signal)

## Note
Internal data management in 'snpEnrichment' use RefSNP (rs) IDs.



## Example
### 1. Install and load snpEnrichment
To install the stable version from CRAN, simply run the following from an R console:
```r
install.packages("snpEnrichment")
```
To install the latest development builds directly from GitHub, run this instead:
```r
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("mcanouil/snpEnrichment")
```
Load *snpEnrichment*
```r
library(snpEnrichment)
```

### 2. Data preparation
```r
snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment")
signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment")
initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)
writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL,
        ldThresh = 0.8, depth = 1000, mc.cores = 1)
```

### 3. Reading data
```r
snpListDir <- system.file("extdata/List", package = "snpEnrichment")
data(transcript)
transcriptFile <- transcript

toy_M1 <- readEnrichment(
    pattern = "Chrom", signalFile, transcriptFile, snpListDir, snpInfoDir,
    distThresh = 1000, sigThresh = 0.05, LD = FALSE, ldDir = NULL,
    mc.cores = 1
)
toy_M1
```

### 4. Computing results
```r
reSample(object = toy_M1,
         nSample = 10,
         empiricPvalue = FALSE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1,
         onlyGenome = FALSE)
```
OR
```r
data(toyEnrichment)
reSample(object = toyEnrichment,
         nSample = 10,
         empiricPvalue = FALSE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1,
         onlyGenome = FALSE)
```

### 5. Further analysis: Exclude SNP from original list.
```r
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
```
OR
```r
excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment")
```

```r
toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)

compareResults <- compareEnrichment(object.x = toyEnrichment,
                                    object.y = toyEnrichment_exclude,
                                    pattern = "Chrom",
                                    nSample = 10,
                                    empiricPvalue = FALSE,
                                    mc.cores = 1,
                                    onlyGenome = TRUE)
```

### 6. Watch results
```r
show(toy_M1)
summary(toy_M1)

show(toyEnrichment)
summary(toyEnrichment)

show(toyEnrichment_exclude)
summary(toyEnrichment_exclude)
```
