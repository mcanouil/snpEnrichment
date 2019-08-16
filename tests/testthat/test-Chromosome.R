expect_s4_class(chromosome(), "Chromosome")
expect_false(is.chromosome(list()))
expect_false(is.chromosome(1))
expect_true(is.chromosome(chromosome()))
expect_true(all(is.chromosome(c(chromosome(), chromosome()))))
expect_equivalent(is.chromosome(list(chromosome(), b = "char")), c(TRUE, FALSE))
expect_equivalent(is.chromosome(c(chromosome(), b = list(12, chromosome()))), c(TRUE, FALSE, TRUE))
expect_equivalent(
  print(chromosome(), type = "eSNP"),
  structure(
    c(NA, NA, NA, 0, 0, 0),
    .Dim = c(1L, 6L),
    .Dimnames = list("Chrom:eSNP", c("EnrichmentRatio", "Z", "PValue", "nSample", "TotalSNP", "eSNP"))
  )
)
