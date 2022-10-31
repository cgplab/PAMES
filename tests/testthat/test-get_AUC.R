test_that("'get_AUC' works", {
  auc <- get_AUC(tumor_toy_data, control_toy_data, cores=2)
  expect_type(auc, "list")
  auc <- get_AUC(tumor_toy_data, control_toy_data, cores=2, full_info=TRUE)
  expect_length(auc, 4)
})
