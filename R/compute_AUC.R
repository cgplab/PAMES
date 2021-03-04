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
#' control samples) that are not NA required to analyze a site (range=0-1, default=1).
#' @param return_info If TRUE (default) return a vector else compute all AUC and return a data.frame
#' reporting fraction of NAs in tumor and control tables.
#' @param na_threshold (DEPRECATED) Fraction of NAs (considered independently in tumor and
#' control samples) above which a site will not be selected (default=0).
#' @return A vector of AUC scores (NA if not analyzed) or a data.frame with AUC scores
#' and the fraction of non-NA samples in tumor and control tables.
#' @examples
#' auc_data <- compute_AUC(tumor_toy_data, control_toy_data)
#' @importFrom stats setNames
#' @export
compute_AUC <- function(tumor_table, control_table, ncores=1, na_threshold, return_info=FALSE, min_samples_frac=1) {
    message(sprintf("[%s] # Compute AUC #", Sys.time()))
    # check parameters
    if (!missing(na_threshold)){
        warning("'na_threshold' is deprecated, use 'min_samples_frac' instead.")
        assertthat::assert_that(is.numeric(na_threshold))
        assertthat::assert_that(dplyr::between(na_threshold, 0, 1))
        min_samples_frac <- 1-na_threshold
    }

    assertthat::assert_that(nrow(tumor_table) == nrow(control_table))
    if (!is.null(rownames(tumor_table)) & !is.null(rownames(control_table))){
        if (any(rownames(tumor_table) != rownames(control_table))){
            warning("tumor_table and control_table have different rownames")
        }
    }

    assertthat::assert_that(is.numeric(ncores))
    ncores <- min(max(ncores, 1), parallel::detectCores())
    assertthat::assert_that(is.numeric(min_samples_frac), dplyr::between(min_samples_frac, 0, 1))
    assertthat::assert_that(is.logical(return_info))

    beta_table <- as.matrix(cbind(tumor_table, control_table))
    beta_table <- round(beta_table)
    storage.mode(beta_table) <- "integer"

    is_tumor <- c(rep(TRUE, ncol(tumor_table)), rep(FALSE, ncol(control_table)))

    cl <- parallel::makeCluster(ncores)
    if (isFALSE(return_info)){
        # select rows by NAs
        message(sprintf("[%s] Filter sites with fraction of available beta-scores greater than or equal to %.2f...", Sys.time(), min_samples_frac))
        tumor_available_sites <- which(rowSums(!is.na(beta_table[,is_tumor]))/sum(is_tumor) >= min_samples_frac)
        control_available_sites <- which(rowSums(!is.na(beta_table[,!is_tumor]))/sum(!is_tumor) >= min_samples_frac)
        available_sites <- intersect(tumor_available_sites, control_available_sites)

        message(sprintf("[%s] Computing...", Sys.time()))
        auc <- setNames(rep(NA_real_, nrow(beta_table)), rownames(beta_table))
        auc[available_sites] <- parallel::parApply(cl, beta_table[available_sites,,drop=FALSE], 1, single_AUC, is_tumor=is_tumor)

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
