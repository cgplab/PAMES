context("AUC")  ##################################################
test_that("'compute_AUC' errors and warnings", {
  expect_error(compute_AUC(tumor_toy_data[1:10,]/100, control_toy_data[1:10,]), "tumor table")
  expect_error(compute_AUC(tumor_toy_data[1:10,]/100, control_toy_data[1:10,]), "tumor table")
  expect_error(compute_AUC(tumor_toy_data[1:10,], control_toy_data[1:10,]/100), "control table")
})
test_that("'compute_AUC' works", {
  auc <- compute_AUC(tumor_toy_data, control_toy_data)
  expect_type(auc, "double")
  # x <- matrix(NA, 10,10); x[,1:2] <- c(100,90)
  # y <- matrix(NA, 10,10); y[,1:2] <- c(0,10)
  # auc <- compute_AUC(x, y, 1, 1)
})

context("selection of sites") ##################################################
test_that("'selection_of_sites' errors and warnings", {
  auc <- runif(nrow(tumor_toy_data))
  expect_error(select_informative_sites(tumor_toy_data, auc, max_sites = 19, illumina27k_hg19[3:4]), "method is set to 'even'")
  expect_error(select_informative_sites(tumor_toy_data, auc, percentiles = c(0,1000), illumina27k_hg19[3:4]), "not true")
  expect_error(select_informative_sites_ext(tumor_toy_data, control_toy_data, auc, illumina27k_hg19), "not a numeric or integer vector")
})
test_that("'selection_of_sites' works", {
  set.seed(252)
  auc <- runif(nrow(tumor_toy_data))
  site_list <- select_informative_sites(tumor_toy_data, auc, illumina27k_hg19[3:4], method = "even")
  expect_type(site_list, "list")
  site_list <- select_informative_sites(tumor_toy_data, auc, illumina27k_hg19[3:4], method = "top")
  expect_type(site_list, "list")
  site_list <- select_informative_sites(tumor_toy_data, auc, illumina27k_hg19[3:4], method = "hyper")
  expect_length(site_list$hyper, 18)
  site_list <- select_informative_sites(tumor_toy_data, auc, illumina27k_hg19[3:4], method = "hypo")
  expect_length(site_list$hypo, 14)
  site_list <- select_informative_sites(tumor_toy_data, auc, illumina27k_hg19[3:4], percentiles=c(5,95))
  expect_length(site_list$hyper, 10)
  site_list <- select_informative_sites_ext(tumor_toy_data, control_toy_data, auc, illumina27k_hg19[3:4], return_info = TRUE)
  expect_length(site_list, 3)
  expect_length(site_list$hypo, 1)
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
  reduced_data <- reduce_to_regions(bs_toy_matrix, bs_toy_sites, cpg_islands)
  expect_is(reduced_data, "matrix")
  expect_equal(nrow(reduced_data), nrow(cpg_islands))
})
test_that("select_informative_regions works", {
  reduced_data <- round(reduce_to_regions(bs_toy_matrix, bs_toy_sites, cpg_islands))
  set.seed(252)
  auc <- runif(nrow(reduced_data))
  # auc <- compute_AUC(reduced_tumor, reduced_control)
  region_list <- select_informative_regions(reduced_data, auc)
  expect_type(region_list, "list")
  region_list <- select_informative_regions(reduced_data, auc, method="top")
  expect_type(region_list, "list")
  region_list <- select_informative_regions(reduced_data, auc, method="hyper")
  # expect_length(region_list$hyper, 20)
  expect_type(region_list, "list")
  region_list <- select_informative_regions(reduced_data, auc, method="hypo")
  # expect_length(region_list$hypo, 20)
  expect_type(region_list, "list")
  region_list <- select_informative_regions(reduced_data, auc, percentiles=c(5,95))
  # expect_length(region_list$hypo, 10)
  expect_type(region_list, "list")
  reduced_tumor <- reduced_data[,1:10]
  reduced_control <- reduced_data[,11:20]
  region_list <- select_informative_regions_ext(reduced_tumor, reduced_control, auc, percentiles=c(5,95))
  # expect_length(region_list$hypo, 10)
  expect_type(region_list, "list")
})

context("Compute purity") ##################################################
test_that("compute_purity works", {
  set.seed(252)
  auc <- runif(nrow(tumor_toy_data))
  site_list <- select_informative_sites(tumor_toy_data, auc, illumina27k_hg19[3:4])
  expect_error(compute_purity(tumor_toy_data*10, site_list, illumina27k_hg19[3:4]), "Unexpected range of beta")
  purity <- compute_purity(tumor_toy_data, site_list, illumina27k_hg19[3:4])
  expect_length(purity, ncol(tumor_toy_data))
})
