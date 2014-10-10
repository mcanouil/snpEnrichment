## Example
### 1. Install and load snpEnrichment
# To install the stable version from CRAN, simply run the following from an R console:
install.packages("snpEnrichment")

# To install the latest development builds directly from GitHub, run this instead:
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("mcanouil/snpEnrichment")

# Load *snpEnrichment*:
library(snpEnrichment)


### 2. Prepare data
snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment")
signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment")

initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)

writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL,
        ldThresh = 0.8, depth = 1000, mc.cores = 1)


### 3. Read data
snpListDir <- system.file("extdata/List", package = "snpEnrichment")
data(transcript)
transcriptFile <- transcript

toyData <- readEnrichment(
    pattern = "Chrom", signalFile, transcriptFile, snpListDir, snpInfoDir,
    distThresh = 1000, sigThresh = 0.05, LD = TRUE, ldDir = NULL,
    mc.cores = 1
)
toyData


### 4. Compute results
reSample(object = toyData,
         nSample = 10,
         empiricPvalue = TRUE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1,
         onlyGenome = TRUE)


### 5. Further analysis: Exclude SNP from original list.
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
excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment")

toyData_exclude <- excludeSNP(toyData, excludeFile, mc.cores = 1)

# Warning, compareEnrichment is in development!!
compareResults <- compareEnrichment(object.x = toyData,
                                    object.y = toyData_exclude,
                                    pattern = "Chrom",
                                    nSample = 10,
                                    empiricPvalue = TRUE,
                                    mc.cores = 1,
                                    onlyGenome = TRUE)


### 6. Watch results
show(toyData)
print(toyData)

show(toyData_exclude)
print(toyData_exclude)
