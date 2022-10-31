test_that("summarise_region works", {
  x <- rbind(rep(NA, 10), sample(100, 10), sample(100, 10))
  x[,3] <- NA
  expect_type(summarise_region(x, 2, "median"), "double")
  expect_true(all(is.na(summarise_region(x, 3, "mean"))))
})
test_that("reduce_to_regions works", {
  expect_error(reduce_to_regions(tumor_toy_data, illumina27k_hg38, cpg_islands), "No shared chromosomes")
  reduced_data <- reduce_to_regions(bs_seq_toy_matrix, bs_seq_toy_sites, cpg_islands)
  expect_is(reduced_data, "matrix")
  expect_equal(nrow(reduced_data), nrow(cpg_islands))
})
