#' Compute Purity of Tumor Samples
#'
#' Estimate the proportion of cancer cells in the admixture of cells
#' forming tumor microenvironment.
#'
#' @param tumor_table A matrix of beta-values from tumor samples.
#' @param list_of_sites A list of indexes generated by
#' \code{\link{select_informative_sites}} or
#' \code{\link{select_informative_regions}}.
#' @param ref_table Reference table used to find sites.
#' @return A vector of purity estimates.
#' @importFrom stats median
#' @export
#' @examples
#' purity <- compute_purity(tumor_toy_data,
#'     list_of_sites=list(hyper=c(1, 10, 20), hypo=c(15,30,45),
#'     ref_table=illumina27k_hg19[3:4]))
compute_purity <- function(tumor_table, list_of_sites, ref_table) {
    message(sprintf("[%s] # Compute purity #", Sys.time()))
    # check parameters
    assertthat::assert_that(nrow(tumor_table) == nrow(ref_table))
    assertthat::assert_that(is.list(list_of_sites))
    assertthat::assert_that(any(c("hyper", "hypo") %in% names(list_of_sites)))
    message(sprintf("- Using %i hyper- and %i hypo-methylated sites",
                    length(list_of_sites[["hyper"]]),
                    length(list_of_sites[["hypo"]])))

    beta_values <- rbind(tumor_table[list_of_sites[["hyper"]],],
                         100 - tumor_table[list_of_sites[["hypo"]],])
    purity <- apply(beta_values, 2, median, na.rm = TRUE)
    message(sprintf("[%s] Done",  Sys.time()))
    return(purity)
}
