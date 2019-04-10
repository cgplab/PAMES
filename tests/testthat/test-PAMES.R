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
  set.seed(252)
  auc <- runif(nrow(tumor_toy_data))
  expect_error(select_informative_sites(tumor_toy_data, auc, max_sites = 19, platform="27"), "'method' is set to 'even'")
  site_list <- select_informative_sites(tumor_toy_data, auc, platform="27", method = "even")
  expect_type(site_list, "list")
  site_list <- select_informative_sites(tumor_toy_data, auc, platform="27", method = "top")
  expect_type(site_list, "list")
  site_list <- select_informative_sites(tumor_toy_data, auc, platform="27", method = "hyper")
  expect_length(site_list$hyper, 18)
  site_list <- select_informative_sites(tumor_toy_data, auc, platform="27", method = "hypo")
  expect_length(site_list$hypo, 14)
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
test_that("select_informative_regions works", {
  reduced_tumor <- reduce_to_regions(bs_toy_matrix, bs_toy_sites, cpg_islands)
  set.seed(252)
  auc <- runif(nrow(reduced_tumor))
  region_list <- select_informative_regions(reduced_tumor, auc)
  expect_type(region_list, "list")
  region_list <- select_informative_regions(reduced_tumor, auc, method="top")
  expect_type(region_list, "list")
  region_list <- select_informative_regions(reduced_tumor, auc, method="hyper")
  expect_length(region_list$hyper, 20)
  region_list <- select_informative_regions(reduced_tumor, auc, method="hypo")
  expect_length(region_list$hypo, 20)
})

context("Compute purity") ##################################################
test_that("compute_purity works", {
  set.seed(252)
  auc <- runif(nrow(tumor_toy_data))
  site_list <- select_informative_sites(tumor_toy_data, auc, platform="27")
  expect_error(compute_purity(tumor_toy_data*10, site_list), "Unexpected range")
  purity <- compute_purity(tumor_toy_data, site_list)
  expect_length(purity, ncol(tumor_toy_data))
  purity <- compute_purity(tumor_toy_data/100, site_list)
  expect_true(all(purity < 1))
})
