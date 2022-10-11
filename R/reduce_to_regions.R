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
#' @param method Take `median` or `mean`of CpG sites.
#' @return A matrix of beta values (nrow == length(cpg_indexes)).
#' @importFrom stats median
#' @export
#' @examples
#' reduced_data <- reduce_to_regions(bs_toy_matrix, bs_toy_sites, cpg_islands[1:10,])
reduce_to_regions <- function(beta_table, cpg_sites, cpg_regions, min_CpGs = 3, method=c("median", "mean")){

    message(sprintf("[%s] # Reduce to regions #", Sys.time()))
    # check parameters
    method <- match.arg(method)
    min_CpGs <- as.integer(min_CpGs)
    assertthat::assert_that(min_CpGs > 0)
    assertthat::assert_that(ncol(cpg_sites) >= 2)
    assertthat::assert_that(ncol(cpg_regions) >= 3)
    assertthat::assert_that(length(intersect(cpg_regions[[1]], cpg_sites[[1]])) > 0,
                            msg="No shared chromosomes between cpg_sites and cpg_regions. Check chromosome names.")
    beta_table <- round(as.matrix(beta_table))
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
    if (length(overlaps) == 0){
        message(sprintf("[%s] No overlaps found",  Sys.time()))
        return(NULL)
    }

    message(sprintf("[%s] Reducing beta values...",  Sys.time()))
    reduced_table <- matrix(NA, length(regions_range), ncol(beta_table),
                            dimnames = list(names(regions_range), colnames(beta_table)))
    # insert reduced beta values at appropriate positions (leave uncovered regions to NA)

    subject_hits <- S4Vectors::subjectHits(overlaps)
    query_hits <- S4Vectors::queryHits(overlaps)
    idx_list <- tapply(subject_hits, query_hits, function(idx) return(idx))
    above_thr_regions <- which(sapply(idx_list, length) >= min_CpGs) # skip regions with not enough sites

    if (length(above_thr_regions) == 0){
        message(sprintf("[%s] No region overlaps enough sites",  Sys.time()))
        return(reduced_table)
    }

    pb <- utils::txtProgressBar(max=length(idx_list), style=3, width=80)
    reduced_data <- lapply(seq_along(idx_list[above_thr_regions]), function(i) {
            utils::setTxtProgressBar(pb, i)
            idx <- idx_list[[i]]
            summarise_region(beta_table[idx,,drop = FALSE], min_CpGs, method)
    })
    close(pb)

    reduced_table[unique(query_hits)[above_thr_regions],] <- do.call(rbind, reduced_data)

    message(sprintf("[%s] Done",  Sys.time()))
    return(reduced_table)
}

#' Reduce many CpG sites to one CpG region
#'
#' If the number of sites is sufficient, take the median/mean value else return NA.
#' @param x A subset matrix.
#' @param n Minimum required number of sites per region (return NA otherwise).
#' @param method Either `median` or `mean`.
#' @return A vector
#' @keywords internal
summarise_region <- function(x, n, method) {
    # remove fully NA sites
    valid_sites <- which(rowSums(is.na(x)) != ncol(x))
    x <- x[valid_sites,,drop=FALSE]
    if (nrow(x) < n) {
        return(rep(NA, ncol(x)))
    } else if (method=="median"){
        return(apply(x, 2, median, na.rm = TRUE))
    } else if (method=="mean"){
        return(apply(x, 2, mean, na.rm = TRUE))
    }
}
