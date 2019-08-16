# expect_true(object = {data(toyEnrichment); grepl("toyEnrichment", ls())})
# expect_s4_class({toyEnrich <- enrichment(); toyEnrich}, "Enrichment")
# expect_invisible(toyEnrich["Loss"] <- toyEnrichment["Loss"])
# expect_s3_class(toyEnrich["Loss"], "data.frame")
# expect_s4_class(
#   object = {
#     toyEnrich <- enrichment(Loss = toyEnrichment["Loss"], eSNP = toyEnrichment["eSNP"])
#     toyEnrich
#   },
#   class = "Enrichment"
# )
# expect_s4_class(
#   object = {
#     toyEnrich <- enrichment(Loss = toyEnrichment["Loss"])
#     toyEnrich
#   },
#   class = "Enrichment"
# )
# expect_type(print(toyEnrich, what = "All"), "list")
# expect_(data(toyEnrichment))
# expect_(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)))
# expect_(print(toyEnrichment, what = "All"))
# expect_(print(toyEnrichment, what = "Genome"))
# expect_(print(toyEnrichment, what = 1))
# expect_(print(toyEnrichment, what = c(2, 4)))
# expect_(print(toyEnrichment, what = seq(22)))
# expect_(plot(toyEnrichment, what = "Genome", type = c("eSNP", "xSNP")))
# expect_(plot(toyEnrichment, what = "Genome", type = "xSNP"))
# expect_(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP")))
# expect_(plot(toyEnrichment, what = 22, type = "eSNP"))
# expect_(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP"), ggplot = TRUE))
# expect_(plot(toyEnrichment, what = 22, type = "eSNP", ggplot = TRUE))
#
# expect_(data(toyEnrichment))
# expect_(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = TRUE)))
# expect_(print(toyEnrichment, what = "All"))
# expect_(print(toyEnrichment, what = "Genome"))
# expect_(print(toyEnrichment, what = 1))
# expect_(print(toyEnrichment, what = c(2, 4)))
# expect_(print(toyEnrichment, what = seq(22)))
# expect_(plot(toyEnrichment, what = "Genome", type = c("eSNP", "xSNP")))
# expect_(plot(toyEnrichment, what = "Genome", type = "xSNP"))
# expect_(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP")))
#
# expect_(plot(toyEnrichment, what = 22, type = "eSNP"))
# expect_(a <- capture.output(dev.off()))
# expect_(snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment"))
# expect_(snpListDir <- system.file("extdata/List", package = "snpEnrichment"))
# expect_(signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment"))
# expect_(data(transcript))
# expect_(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)))
# expect_(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = FALSE, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = FALSE, ldDir = NULL, mc.cores = 1)))
# expect_(print(toyEnrichment, what = "All"))
# expect_(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)))
#
# expect_(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = FALSE, ldDir = NULL, mc.cores = 1)))
# expect_(print(toyEnrichment, what = "All"))
# expect_(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)))
# expect_(trash <- capture.output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1)))
# expect_(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 1)))
# expect_(print(toyEnrichment, what = "All"))
# expect_(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)))
# expect_(trash <- capture.output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1)))
# expect_(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 1)))
# expect_(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)))
#
# expect_(print(toyEnrichment, what = "All"))
# expect_(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)))
# expect_(trash <- capture.output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1)))
# expect_(trash <- capture.output(toyEnrichmentM2 <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 1)))
# expect_(trash <- capture.output(reSample(object = toyEnrichmentM2, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)))
# expect_(print(toyEnrichmentM2, what = "All"))
# expect_(data(toyEnrichment))
# expect_(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
# expect_(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)))
# expect_(print(toyM1_exclude, what = "All"))
#
# expect_(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
# expect_(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)))
# expect_(print(toyM1_exclude, what = "All"))
# expect_(excludeFile <- c(
#   "rs4376885", "rs17690928", "rs6460708", "rs13061537", "rs11769827",
#   "rs12717054", "rs2907627", "rs1380109", "rs7024214", "rs7711972",
#   "rs9658282", "rs11750720", "rs1793268", "rs774568", "rs6921786",
#   "rs1699031", "rs6994771", "rs16926670", "rs465612", "rs3012084",
#   "rs354850", "rs12803455", "rs13384873", "rs4364668", "rs8181047",
#   "rs2179993", "rs12049335", "rs6079926", "rs2175144", "rs11564427",
#   "rs7786389", "rs7005565", "rs17423335", "rs12474102", "rs191314",
#   "rs10513168", "rs1711437", "rs1992620", "rs283115", "rs10754563",
#   "rs10851727", "rs2173191", "rs7661353", "rs1342113", "rs7042073",
#   "rs1567445", "rs10120375", "rs550060", "rs3761218", "rs4512977"
# ))
# expect_(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)))
# expect_(print(toyM1_exclude, what = "All"))
# expect_(data(toyEnrichment))
# expect_(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)))
# expect_(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
# expect_(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)))
#
# expect_(trash <- capture.output(reSample(object = toyM1_exclude, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)))
# expect_(trash <- capture.output(res_toyM1 <- compareEnrichment(object.x = toyEnrichment, object.y = toyM1_exclude, pattern = "Chrom", nSample = 10, empiricPvalue = FALSE, mc.cores = 1, onlyGenome = FALSE)))
# expect_(data(toyEnrichment))
# expect_(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)))
# expect_(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
# expect_(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)))
# expect_(trash <- capture.output(reSample(object = toyM1_exclude, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)))
# expect_(trash <- capture.output(reSample(object = toyM1_exclude, nSample = 10, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)))
# expect_(trash <- capture.output(res_toyM1 <- compareEnrichment(object.x = toyEnrichment, object.y = toyM1_exclude, pattern = "Chrom", nSample = 10, empiricPvalue = FALSE, mc.cores = 1, onlyGenome = FALSE)))
# expect_(snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment"))
#
# expect_(snpListDir <- system.file("extdata/List", package = "snpEnrichment"))
# expect_(signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment"))
# expect_(data(transcript))
# expect_(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 2)))
# expect_(trash <- capture.output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 2)))
# expect_(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 2)))
# expect_(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 2, onlyGenome = FALSE)))
# expect_(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"))
# expect_(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 2)))
# expect_(trash <- capture.output(res_toyM1 <- compareEnrichment(object.x = toyEnrichment, object.y = toyM1_exclude, pattern = "Chrom", nSample = 10, empiricPvalue = FALSE, mc.cores = 2, onlyGenome = FALSE)))
#
# expect_(data(toyEnrichment))
# expect_(trash <- getEnrichSNP(toyEnrichment, type = "eSNP"))
# expect_(trash <- getEnrichSNP(toyEnrichment, type = "xSNP"))
