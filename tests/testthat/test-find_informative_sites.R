test_that("'find_informative_sites' works", {
  set.seed(252)
  auc1 <- get_AUC(tumor_toy_data, control_toy_data, cores = 2)
  info_sites <- find_informative_sites(tumor_toy_data, control_toy_data, auc1, illumina27k_hg19, cores=2, percentiles=c(5,95), full_info = FALSE)
  expect_type(info_sites, "list")
})
