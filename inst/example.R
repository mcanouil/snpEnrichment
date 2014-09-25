#######################
### 1. Data Preparation
snpInfoDir <- system.file("extdata/snpInfo",
                          package = "snpEnrichment")
snpListDir <- system.file("extdata/List",
                          package = "snpEnrichment")
signalFile <- system.file("extdata/Signal/toySignal.txt",
                          package = "snpEnrichment")
initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)
writeLD(pattern = "Chrom", snpInfoDir, signalFile,
        ldDir = NULL, ldThresh = 0.8, depth = 1000,
        mc.cores = 1)

###################
### 2. Reading data
data(transcript)
transcriptFile <- transcript

initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)
toy_M1 <- readEnrichment(pattern = "Chrom", signalFile,
                         transcriptFile, snpListDir,
                         snpInfoDir, distThresh = 1000,
                         sigThresh = 0.05, LD = FALSE,
                         ldDir = NULL, mc.cores = 1)
toy_M1


########################
### 3. Computing results
reSample(object = toy_M1,
         nSample = 10,
         empiricPvalue = FALSE,
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         mc.cores = 1,
         onlyGenome = FALSE)


#######################
### 5. Further analysis
## Exclude SNP from original list.
excludeFile <- system.file("extdata/Exclude/toyExclude.txt",
                           package = "snpEnrichment")

data(toyEnrichment)
toyEnrichment_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)

compareResults <- compareEnrichment(object.x = toyEnrichment, # Cannot compute if LD is missing, i'm working on it.
                                    object.y = toyEnrichment_exclude,
                                    pattern = "Chrom",
                                    nSample = 10,
                                    empiricPvalue = FALSE,
                                    mc.cores = 1,
                                    onlyGenome = TRUE)


####################
### 4. Watch results
show(toyEnrichment)
print(toyEnrichment)

show(toyEnrichment_exclude)
print(toyEnrichment_exclude)}
