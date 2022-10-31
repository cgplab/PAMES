test_that("compute_purity works", {
  auc <- get_AUC(tumor_toy_data, control_toy_data, cores=2)
  site_list <- find_informative_sites(tumor_toy_data, control_toy_data, auc, illumina27k_hg19, cores=2)
  purity <- get_purity(tumor_toy_data, site_list)
  expect_type(purity, "list")
})
