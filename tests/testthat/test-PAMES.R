context("Test PAMES")

test_that("'compute_AUC' works", {
    auc <- compute_AUC(tumor_toy_data, control_toy_data)
    expect_is(auc, "numeric")
})

test_that("'selection_of_sites' returns a list with 2 vectors", {
    x <- as.matrix(runif(27578))
    y <- as.matrix(runif(27578))
    site_list <- select_informative_sites(x, y, platform="27k")
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

test_that("reduce_island works", {
  res <- reduce_island(bs_tumor_toy_data[bs_toy_indexes[[1]], ], 3)
  expect_length(res, 10)
  expect_is(res, "numeric")
})

test_that("reduce_to_islands", {
  reduced_tumor <- reduce_to_islands(bs_tumor_toy_data, bs_toy_indexes)
  expect_is(reduced_tumor, "matrix")
  expect_equal(nrow(reduced_tumor), length(bs_toy_indexes))
})
