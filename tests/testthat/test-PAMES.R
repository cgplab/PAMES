context("AUC")
test_that("'compute_AUC' errors and warnings", {
  skip("AUC")
  expect_error(compute_AUC(tumor_toy_data[1:10,]/100, control_toy_data[1:10,]),
    "computation efficiency")
  expect_error(compute_AUC(tumor_toy_data[1:10,], control_toy_data[1:10,]/100),
    "computation efficiency")
  expect_warning(compute_AUC(tumor=tumor_toy_data[1:10,],
      control=control_toy_data[1:10,]), "deprecated")
})
test_that("'compute_AUC' works", {
  skip("AUC")
  auc <- compute_AUC(tumor_toy_data, control_toy_data, ncores=4)
  expect_is(auc, "numeric")
})

context("selection of sites")
test_that("'selection_of_sites' errors and warnings", {
  skip("sel of sites")
  expect_warning(select_informative_sites(tumor=tumor_toy_data,
    auc=runif(nrow(tumor_toy_data)), platform="27"), "deprecated")
  expect_error(select_informative_sites(tumor_toy_data,
    runif(nrow(tumor_toy_data)), platform="27", max_sites=19), "even")
})
test_that("'selection_of_sites' works", {
  skip("sel of sites")
  auc <- round(runif(27578))
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
  expect_true(too_close(4722, 14309, 1e6, illumina27k_hg19))
  expect_false(too_close(1, 2:10, 1e6, illumina27k_hg19))
})

context("Bisulphite sequencing")
test_that("reduce_to_island works", {
  idx <- unlist(bs_toy_indexes[1:10])
  reduced_tumor <-
    reduce_to_islands(bs_toy_matrix[idx, 1:10], bs_toy_indexes[1:10], 3)
  expect_equal(nrow(reduced_tumor), 10)
  expect_is(reduced_tumor, "matrix")
})
test_that("compute_islands_indexes works", {
  cpg_indexes <- compute_islands_indexes(bs_toy_sites, head(cpg_islands_df, 10))
  expect_is(cpg_indexes, "list")
  expect_equal(length(cpg_indexes), 10)
})

context("Compute purity")
test_that("purity sites", {
  compute_gg
})
test_that("purity islands", {
})
