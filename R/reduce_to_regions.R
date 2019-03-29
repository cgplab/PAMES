#' Reduce beta values from CpG sites to genomic regions
#'
#' Reduce several beta values from different CpG sites to single beta values
#' associated to one genomic region (such as CpG islands) where CpG sites are located.
#' Different technologies retrieve DNA methylation levels of
#' different CpG sites. This function makes a direct comparison possible.
#'
#' @param beta_table A matrix of beta-values (percentage).
#' @param cpg_sites A data.frame reporting the genomic location of CpG sites.
#' @param cpg_regions A data.frame reporting the genomic location of genomic regions.
#' @param min_CpGs An integer (default to 3). Minimum number of CpG sites
#' within a single genomic region required to compute the reduced beta value
#' (return NA otherwise).
#' @return A matrix of beta values (nrow == length(cpg_indexes)).
#' @importFrom stats median
#' @export
#' @examples
#' reduced_data <- reduce_to_regions(bs_toy_matrix, bs_toy_sites, cpg_islands[1:10,])
reduce_to_regions <- function(beta_table, cpg_sites, cpg_regions, min_CpGs = 3){
  message(sprintf("[%s] # Reduce to regions #", Sys.time()))
  # check parameters
  min_CpGs <- as.integer(min_CpGs)
  assertthat::assert_that(min_CpGs > 0)
  assertthat::assert_that(ncol(cpg_sites) >= 2)
  assertthat::assert_that(ncol(cpg_regions) >= 3)
  assertthat::assert_that(length(intersect(cpg_regions[[1]], cpg_sites[[1]])) > 0,
    msg="No shared chromosomes between cpg_sites and cpg_regions. Check chromosome names.")
  diff_range <- diff(range(beta_table, na.rm = TRUE))
  assertthat::assert_that(diff_range > 1, diff_range <= 100,
    msg=paste("For computation efficiency convert table to percentage values."))
  beta_table <- as.matrix(beta_table)
  beta_table <- round(beta_table)
  storage.mode(beta_table) <- "integer"

  if (ncol(cpg_regions) < 4) {
    cpg_regions[[4]] <- paste0("CpG_", seq_len(nrow(cpg_regions)))
  }

  message(sprintf("[%s] Finding overlaps...",  Sys.time()))
  regions_range <- GenomicRanges::GRanges(seqnames = cpg_regions[[1]],
                ranges = IRanges::IRanges(start = cpg_regions[[2]],
                                          end = cpg_regions[[3]],
                                          names = cpg_regions[[4]]))
  sites_range <- GenomicRanges::GRanges(cpg_sites[[1]], ranges = IRanges::IRanges(start=cpg_sites[[2]], width=1))
  overlaps <- GenomicRanges::findOverlaps(regions_range, sites_range)

  message(sprintf("[%s] Reducing beta values...",  Sys.time()))
  reduced_table <- do.call(rbind, tapply(overlaps@to, overlaps@from,
          function(idx) {median_of_region(beta_table[idx, , drop = FALSE], min_CpGs)}))

  rownames(reduced_table) <- names(regions_range)[unique(overlaps@from)]

  message(sprintf("[%s] Done",  Sys.time()))
  return(reduced_table)
}

#' Transform CpG sites to one CpG region
#'
#' If the number of sites is sufficient take the median value else return NA.
#' @param x A subset matrix.
#' @param n Minimum required number of sites per region (return NA otherwise).
#' @return A vector
#' @keywords internal
median_of_region <- function(x, n) {
    # remove sites non reported for all samples
    valid_sites <- which(!rowSums(is.na(x)) == ncol(x))
    x <- x[valid_sites, , drop=FALSE]
    if (nrow(x) < n) {
      return(NA)
    } else {
      return(apply(x, 2, median, na.rm = TRUE))
    }
}
