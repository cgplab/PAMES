test_that("select_informative_regions works", {
  reduced_data <- round(reduce_to_regions(bs_seq_toy_matrix, bs_seq_toy_sites, cpg_islands))
  reduced_tumor <- reduced_data[,1:10]
  reduced_control <- reduced_data[,11:20]
  auc_bs_seq <- get_AUC(reduced_tumor, reduced_control, cores = 2)
  region_list <- find_informative_regions(reduced_tumor, reduced_control, auc_bs_seq, cores=2, percentiles=c(5,95))
  expect_type(region_list, "list")
})
