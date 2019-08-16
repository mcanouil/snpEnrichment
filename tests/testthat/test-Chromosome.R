expect_s4_class(snpEnrichment:::chromosome(), "Chromosome")
expect_false(snpEnrichment:::is.chromosome(list()))
expect_false(snpEnrichment:::is.chromosome(1))
expect_true(snpEnrichment:::is.chromosome(snpEnrichment:::chromosome()))
expect_true(all(snpEnrichment:::is.chromosome(c(snpEnrichment:::chromosome(), snpEnrichment:::chromosome()))))
expect_equivalent(snpEnrichment:::is.chromosome(list(snpEnrichment:::chromosome(), b = "char")), c(TRUE, FALSE))
expect_equivalent(snpEnrichment:::is.chromosome(c(snpEnrichment:::chromosome(), b = list(12, snpEnrichment:::chromosome()))), c(TRUE, FALSE, TRUE))
expect_equivalent(
  print(snpEnrichment:::chromosome(), type = "eSNP"),
  structure(
    c(NA, NA, NA, 0, 0, 0),
    .Dim = c(1L, 6L),
    .Dimnames = list("Chrom:eSNP", c("EnrichmentRatio", "Z", "PValue", "nSample", "TotalSNP", "eSNP"))
  )
)
