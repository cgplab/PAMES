context("AUC")  ##################################################
test_that("'compute_AUC' errors and warnings", {
  expect_error(compute_AUC(tumor_toy_data[1:10,]/100, control_toy_data[1:10,]), "tumor table")
  expect_error(compute_AUC(tumor_toy_data[1:10,]/100, control_toy_data[1:10,]), "tumor table")
  expect_error(compute_AUC(tumor_toy_data[1:10,], control_toy_data[1:10,]/100), "control table")
})
test_that("'compute_AUC' works", {
  auc <- compute_AUC(tumor_toy_data, control_toy_data)
  expect_type(auc, "double")
})

context("selection of sites") ##################################################
test_that("'selection_of_sites' errors and warnings", {
  expect_error(select_informative_sites(tumor_toy_data,
    runif(nrow(tumor_toy_data)), platform="27", max_sites=19), "even")
})
test_that("'selection_of_sites' works", {
  auc <- runif(27578)
  site_list <- select_informative_sites(tumor_toy_data, auc, platform="27")
  expect_type(site_list, "list")
  expect_length(site_list, 2)
})
test_that("'cluster_reduction' works", {
 res <- cluster_reduction(c(4722, 14309, 4522, 19915), 4, 1e6, illumina27k_hg19)
 expect_equal(res, c(4722))

 res <- cluster_reduction(c(4722, 14309, 4522, 19915), 2, 1e4, illumina27k_hg19)
 expect_equal(res, c(4722, 4522))

 res <- cluster_reduction(c(4722, 14309, 4522, 19915), 2, 1e3, illumina27k_hg19)
 expect_equal(res, c(4722, 14309))
})
test_that("'too_close' works", {
  expect_true(is_too_close(4722, 14309, 1e6, illumina27k_hg19))
  expect_false(is_too_close(1, 2:10, 1e6, illumina27k_hg19))
})

context("CpG regions") ##################################################
test_that("median of regions", {
  x <- rbind(rep(NA, 10), sample(100, 10), sample(100, 10))
  x[,3] <- NA
  expect_type(median_of_region(x, 2), "double")
  expect_true(is.na(median_of_region(x, 3)))
})
test_that("reduce_to_regions works", {
  expect_error(reduce_to_regions(tumor_toy_data, illumina27k_hg19, cpg_islands), "No shared chromosomes")

  reduced_tumor <- reduce_to_regions(bs_toy_matrix, bs_toy_sites, cpg_islands)
  expect_is(reduced_tumor, "matrix")
})
