context("Verify selection functions")

test_that("'selection_of_sites' return a list", {
    x <- as.matrix(runif(27578))
    y <- as.matrix(runif(27578))
    site_list <- select_informative_sites(x, y, platform="27k")
    expect_type(site_list, "list")
})
