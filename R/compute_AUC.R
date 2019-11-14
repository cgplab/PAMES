#' Compute Area Under Curve for a matrix of samples
#'
#' This function computes the Area Under Curve used to define the segregation
#' between tumor and control samples accordingly to their methylation (beta)
#' values.
#'
#' @param tumor_table A matrix of beta-values (percentage) from tumor samples.
#' @param control_table A matrix of beta-values (percentage) from normal/control samples.
#' @param ncores Number of parallel processes to use for parallel computing.
#' @param min_samples_frac Fraction of samples (independently in tumor and
#' control samples) that are not NA required to analyze a site (default=100).
#' @param simplify If TRUE (default) return a vector else compute all AUC and return a data.frame
#' reporting fraction of NAs in tumor and control tables.
#' @param na_threshold (DEPRECATED) Fraction of NAs (considered independently in tumor and
#' control samples) above which a site will not be selected (default=0).
#' @return A vector of AUC scores (NA if not analyzed) or a data.frame with AUC scores
#' and the fraction of non-NA samples in tumor and contol tables.
#' @examples
#' auc_data <- compute_AUC(tumor_toy_table, control_toy_table)
#' @export
compute_AUC <- function(tumor_table, control_table, ncores=1, na_threshold, simplify=TRUE, min_samples_frac=1) {
    message(sprintf("[%s] # Compute AUC #", Sys.time()))
    # check parameters
    if (!missing(na_threshold)){
        warning("'na_threshold' is deprecated, use 'min_samples_frac' instead.")
        assertthat::assert_that(is.numeric(na_threshold))
        assertthat::assert_that(dplyr::between(na_threshold, 0, 1))
        min_samples_frac <- 1-na_threshold
    }

    ncores <- min(max(as.integer(ncores), 1), parallel::detectCores())
    assertthat::assert_that(is.numeric(min_samples_frac))
    assertthat::assert_that(is.logical(simplify))
    assertthat::assert_that(dplyr::between(min_samples_frac, 0, 1))
    assertthat::assert_that(nrow(tumor_table) == nrow(control_table))

    diff_range_t <- diff(range(tumor_table, na.rm = TRUE))
    diff_range_c <- diff(range(control_table, na.rm = TRUE))
    assertthat::assert_that(diff_range_t > 1, diff_range_t <= 100, msg="For computation efficiency, convert tumor table to percentage values.")
    assertthat::assert_that(diff_range_c > 1, diff_range_c <= 100, msg="For computation efficiency, convert control table to percentage values.")

    beta_table <- as.matrix(cbind(tumor_table, control_table))
    beta_table <- round(beta_table)
    storage.mode(beta_table) <- "integer"

    is_tumor <- c(rep(TRUE, ncol(tumor_table)), rep(FALSE, ncol(control_table)))

    cl <- parallel::makeCluster(ncores)
    if (isTRUE(simplify)){
        # select rows by NAs
        message(sprintf("[%s] Selecting sites with fraction of valid beta-scores greater than or equal to %.2f...", Sys.time(), min_samples_frac))
        tumor_valid_sites <- which(rowSums(!is.na(beta_table[,is_tumor]))/sum(is_tumor) >= min_samples_frac)
        control_valid_sites <- which(rowSums(!is.na(beta_table[,!is_tumor]))/sum(!is_tumor) >= min_samples_frac)
        valid_sites <- intersect(tumor_valid_sites, control_valid_sites)

        message(sprintf("[%s] Computing...", Sys.time()))
        auc <- setNames(rep(NA_real_, nrow(beta_table)), rownames(beta_table))
        auc[valid_sites] <- parallel::parApply(cl, beta_table[valid_sites,,drop=FALSE], 1, single_AUC, is_tumor=is_tumor)

    } else {
        message(sprintf("[%s] Computing...", Sys.time()))
        auc_scores <- parallel::parApply(cl, beta_table, 1, single_AUC, is_tumor = is_tumor)
        auc <- data.frame(auc = auc_scores,
                          tumor_nonNA_frac = rowSums(!is.na(beta_table[,is_tumor]))/sum(is_tumor),
                          control_nonNA_frac = rowSums(!is.na(beta_table[,!is_tumor]))/sum(!is_tumor))
    }
    parallel::stopCluster(cl)
    message(sprintf("[%s] Done",  Sys.time()))
    return(auc)
}

#' Compute AUC
#'
#' Use Wilcoxon method to compute AUC.
#'
#' @param scores integer vector (range 1-100)
#' @param is_tumor logical vector (class labels)
#' @keywords internal
#' \href{http://blog.revolutionanalytics.com/2017/03/auc-meets-u-stat.html}{http://blog.revolutionanalytics.com}
single_AUC <- function(scores, is_tumor) {
    assertthat::assert_that(is.integer(scores))
    assertthat::assert_that(is.logical(is_tumor))
    na_idx <- is.na(scores)
    scores <- scores[!na_idx]
    is_tumor <- is_tumor[!na_idx]
    n1 <- sum(is_tumor)
    n2 <- sum(!is_tumor)
    R1 <- sum(rank(scores)[is_tumor])
    U1 <- R1 - n1*(n1+1)/2
    return(U1/(n1*n2))
}
