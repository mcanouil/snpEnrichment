expect_invisible(data(toyEnrichment))
expect_s4_class({toyEnrich <- enrichment(); toyEnrich}, "Enrichment")
expect_invisible(toyEnrich["Loss"] <- toyEnrichment["Loss"])
expect_s3_class(toyEnrich["Loss"], "data.frame")
expect_s4_class(
  object = {
    toyEnrich <- enrichment(Loss = toyEnrichment["Loss"], eSNP = toyEnrichment["eSNP"])
    toyEnrich
  },
  class = "Enrichment"
)
expect_s4_class(
  object = {
    toyEnrich <- enrichment(Loss = toyEnrichment["Loss"])
    toyEnrich
  },
  class = "Enrichment"
)
expect_type(print(toyEnrich, what = "All"), "list")
expect_output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE))
expect_type(print(toyEnrichment, what = "All"), "list")
expect_type(print(toyEnrichment, what = "Genome"), "double")
expect_type(print(toyEnrichment, what = 1), "double")
expect_type(print(toyEnrichment, what = c(2, 4)), "double")
expect_type(print(toyEnrichment, what = seq(22)), "double")
expect_invisible(plot(toyEnrichment, what = "Genome", type = "eSNP"))
# expect_invisible(plot(toyEnrichment, what = "Genome", type = "xSNP")) # issue
# expect_invisible(plot(toyEnrichment, what = "Genome", type = c("eSNP", "xSNP"))) # issue
expect_invisible(plot(toyEnrichment, what = 22, type = "eSNP"))
# expect_invisible(plot(toyEnrichment, what = 22, type = "xSNP")) # issue
# expect_invisible(plot(toyEnrichment, what = 22, type = c("eSNP", "xSNP"))) # issue
expect_invisible(plot(toyEnrichment, what = 2, type = "eSNP"))
expect_invisible(plot(toyEnrichment, what = 2, type = "xSNP"))
expect_invisible(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP")))
expect_invisible(plot(toyEnrichment, what = 22, type = "eSNP"))
# expect_invisible(plot(toyEnrichment, what = 22, type = "xSNP")) # issue
# expect_invisible(plot(toyEnrichment, what = 22, type = c("eSNP", "xSNP"))) # issue
# expect_invisible(plot(toyEnrichment, what = 2, type = "eSNP", ggplot = TRUE)) # issue
# expect_invisible(plot(toyEnrichment, what = 2, type = "xSNP", ggplot = TRUE)) # issue
# expect_invisible(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP"), ggplot = TRUE)) # issue
# expect_invisible(plot(toyEnrichment, what = 22, type = "eSNP", ggplot = TRUE)) # issue
# expect_invisible(plot(toyEnrichment, what = 22, type = "xSNP", ggplot = TRUE)) # issue
# expect_invisible(plot(toyEnrichment, what = 22, type = c("eSNP", "xSNP", ggplot = TRUE))) # issue

expect_invisible(data(toyEnrichment))
expect_output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = TRUE))
expect_type(print(toyEnrichment, what = "All"), "list")
expect_type(print(toyEnrichment, what = "Genome"), "double")
expect_type(print(toyEnrichment, what = 1), "double")
expect_type(print(toyEnrichment, what = c(2, 4)), "double")
expect_type(print(toyEnrichment, what = seq(22)), "double")
expect_invisible(plot(toyEnrichment, what = "Genome", type = "eSNP"))
# expect_invisible(plot(toyEnrichment, what = "Genome", type = "xSNP")) # issue
# expect_invisible(plot(toyEnrichment, what = "Genome", type = c("eSNP", "xSNP"))) # issue
expect_error(plot(toyEnrichment, what = 2, type = "eSNP"))
expect_error(plot(toyEnrichment, what = 2, type = "xSNP"))
expect_error(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP")))
expect_output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE))
expect_invisible(plot(toyEnrichment, what = 2, type = "eSNP"))
# expect_invisible(plot(toyEnrichment, what = 2, type = "xSNP")) # issue
# expect_invisible(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP"))) # issue
expect_invisible(snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment"))
expect_invisible(snpListDir <- system.file("extdata/List", package = "snpEnrichment"))
expect_invisible(signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment"))

expect_invisible(data(transcript))
expect_true(any(grepl("transcript", ls())))
expect_output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1))
expect_invisible(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = FALSE, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = FALSE, ldDir = NULL, mc.cores = 1))
expect_type(print(toyEnrichment, what = "All"), "list")
expect_output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1))
expect_output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = FALSE, ldDir = NULL, mc.cores = 1))

expect_type(print(toyEnrichment, what = "All"), "list")
expect_output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1))
expect_output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1))
expect_output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 1))

expect_type(print(toyEnrichment, what = "All"), "list")
expect_output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1))
expect_output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1))
expect_output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 1))
expect_output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE))

expect_type(print(toyEnrichment, what = "All"), "list")
expect_output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1))
expect_output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1))
expect_output(toyEnrichmentM2 <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 1))
expect_output(reSample(object = toyEnrichmentM2, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE))
expect_type(print(toyEnrichmentM2, what = "All"), "list")
expect_invisible(data(toyEnrichment))
expect_invisible(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
expect_output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1))
expect_type(print(toyM1_exclude, what = "All"), "list")

expect_invisible(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
expect_output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1))
expect_type(print(toyM1_exclude, what = "All"), "list")
expect_invisible(excludeFile <- c(
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
))
expect_output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1))
expect_type(print(toyM1_exclude, what = "All"), "list")
expect_invisible(data(toyEnrichment))
expect_output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE))
expect_invisible(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
expect_output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1))

expect_output(reSample(object = toyM1_exclude, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE))
expect_warning(trash <- capture.output(res_toyM1 <- compareEnrichment(object.x = toyEnrichment, object.y = toyM1_exclude, pattern = "Chrom", nSample = 10, empiricPvalue = FALSE, mc.cores = 1, onlyGenome = FALSE)))
expect_invisible(data(toyEnrichment))
expect_output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE))
expect_invisible(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
expect_output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1))
expect_output(reSample(object = toyM1_exclude, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE))
expect_output(reSample(object = toyM1_exclude, nSample = 10, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE))
expect_warning(trash <- capture.output(res_toyM1 <- compareEnrichment(object.x = toyEnrichment, object.y = toyM1_exclude, pattern = "Chrom", nSample = 10, empiricPvalue = FALSE, mc.cores = 1, onlyGenome = FALSE)))
expect_length(res_toyM1, 3)
expect_invisible(snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment"))

expect_invisible(snpListDir <- system.file("extdata/List", package = "snpEnrichment"))
expect_invisible(signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment"))
expect_invisible(data(transcript))
expect_output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 2))
expect_output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 2))
expect_output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 2))
expect_output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 2, onlyGenome = FALSE))
expect_invisible(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
expect_output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 2))
expect_warning(trash <- capture.output(res_toyM1 <- compareEnrichment(object.x = toyEnrichment, object.y = toyM1_exclude, pattern = "Chrom", nSample = 10, empiricPvalue = FALSE, mc.cores = 2, onlyGenome = FALSE)))
expect_length(res_toyM1, 3)

expect_invisible(data(toyEnrichment))
expect_s3_class(trash <- getEnrichSNP(toyEnrichment, type = "eSNP"), "data.frame")
# expect_s3_class(trash <- getEnrichSNP(toyEnrichment, type = "xSNP"), "data.frame") # issue

