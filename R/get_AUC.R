#' Calculate Area Under Curve comparing tumor and control samples
#'
#' This function calculates the Area Under Curve that define the segregation
#' between tumor and control samples accordingly to their methylation (beta)
#' values.
#'
#' @param tumor_table A matrix of beta-values from tumor samples.
#' @param control_table A matrix of beta-values from normal/control samples.
#' @param cores Number of cores for parallel computing.
#' @param max_NA_data At most this fraction of NAs (in tumor and control data
#' independently) is permitted to calculate AUC of a site (default=1, i.e. only
#' sites without NAs).
#' @param full_info If TRUE (default=FALSE) returns also sites infos (fraction
#' of NAs in tumor and control tables) without skipping probes with missing values.
#' @return A data.frame with AUC scores and probes names
#' @examples
#' auc_data <- get_AUC(tumor_toy_data, control_toy_data)
#' @importFrom stats setNames
#' @export
get_AUC <- function(tumor_table, control_table, max_NA_data=0, cores=1, full_info=FALSE) {
    message(sprintf("[%s] # Calculate AUC #", Sys.time()))

    # check parameters
    assertthat::assert_that(is.matrix(tumor_table))
    assertthat::assert_that(is.matrix(control_table))
    assertthat::assert_that(!is.null(rownames(tumor_table)))
    assertthat::assert_that(identical(rownames(tumor_table), rownames(control_table)))
    assertthat::assert_that(is.numeric(max_NA_data))
    assertthat::assert_that(is.numeric(cores))
    assertthat::assert_that(is.logical(full_info))
    assertthat::assert_that(max_NA_data>=0 & max_NA_data<=1)

    cores <- min(round(cores), parallel::detectCores())

    beta_table <- cbind(tumor_table, control_table)

    is_tumor <- c(rep(TRUE, ncol(tumor_table)), rep(FALSE, ncol(control_table)))

    cl <- parallel::makeCluster(cores)
    if (isFALSE(full_info)){
        # identify non-valid sites
        message(sprintf("[%s] Skip sites with proportion of missing beta-scores greater than %.2f...", Sys.time(), max_NA_data))
        tumor_available_sites   <- which(rowSums(is.na(tumor_table))/ncol(tumor_table) >= max_NA_data)
        control_available_sites <- which(rowSums(is.na(control_table))/ncol(control_table) >= max_NA_data)
        available_sites <- intersect(tumor_available_sites, control_available_sites)

        message(sprintf("[%s] Calculating...", Sys.time()))
        auc_scores <- rep(NA_real_, nrow(beta_table))
        auc_scores[available_sites] <- parallel::parApply(cl, beta_table[available_sites,,drop=FALSE], 1, single_AUC, is_tumor=is_tumor)

        AUC_df <- data.frame(Probe=rownames(tumor_table),
                             AUC = auc_scores)

    } else {
        message(sprintf("[%s] Calculating...", Sys.time()))
        auc_scores <- parallel::parApply(cl, beta_table, 1, single_AUC, is_tumor = is_tumor)
        AUC_df <- data.frame(Probe=rownames(tumor_table),
                             AUC = auc_scores,
                             tumor_NA_fraction   = rowSums(is.na(tumor_table))/ncol(tumor_table),
                             control_NA_fraction = rowSums(is.na(control_table))/ncol(control_table))
    }
    parallel::stopCluster(cl)
    message(sprintf("[%s] Done",  Sys.time()))
    return(AUC_df)
}

#' Calculate AUC
#'
#' Use Wilcoxon method to calculate AUC.
#'
#' @param scores integer vector (range 1-100)
#' @param is_tumor logical vector (class labels)
#' @keywords internal
#' \href{http://blog.revolutionanalytics.com/2017/03/auc-meets-u-stat.html}{http://blog.revolutionanalytics.com}
single_AUC <- function(scores, is_tumor) {
    assertthat::assert_that(is.numeric(scores))
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
