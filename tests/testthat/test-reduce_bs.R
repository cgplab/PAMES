context("Verify BS data reducing")

test_that("'remap_bs_data works", {
    n <- 1000
    remapped <- remap_bs_data(bs_control_toy_data[, ],
                              bs_toy_sites[, ],
                              islands_df=cpg_islands[seq_len(n), ],
                              100)
    expect_equal(nrow(remapped), n)
    remapped <- remap_bs_data(bs_control_toy_data[, ],
                              bs_toy_sites[, ])
    expect_equal(nrow(remapped), nrow(cpg_islands))
})
