#' Compute Area Under Curve for a matrix of samples
#'
#' This function computes the Area Under Curve used to define the
#' segregation between tumor and control samples accordingly to their
#' methylation (beta) values.
#'
#' @param tumor_table A matrix of beta-values (percentage) from tumor samples.
#' @param control_table A matrix of beta-values (percentage) from normal/control samples.
#' @param ncores Number of parallel processes to use for parallel computing.
#' @param na_threshold Fraction of NAs (considered independently in tumor and
#' control samples) above which a site will not be selected (default=0).
#' @return A vector of AUC scores.
#' @examples
#' auc_data <- compute_AUC(tumor_toy_data, control_toy_data)
#' @export
compute_AUC <- function(tumor_table, control_table, ncores = 1, na_threshold = 0) {
  message(sprintf("[%s] # Compute AUC #", Sys.time()))
  # check parameters
  ncores <- as.integer(ncores)
  system_cores <- parallel::detectCores()
  assertthat::assert_that(ncores < system_cores)

  na_threshold <- as.numeric(na_threshold)
  assertthat::assert_that(na_threshold >= 0, na_threshold < 1)

  assertthat::assert_that(nrow(tumor_table) == nrow(control_table))
  diff_range_t <- diff(range(tumor_table, na.rm = TRUE))
  assertthat::assert_that(diff_range_t > 1, diff_range_t <= 100, msg="For computation efficiency, convert tumor table to percentage values.")
  diff_range_c <- diff(range(control_table, na.rm = TRUE))
  assertthat::assert_that(diff_range_c > 1, diff_range_c <= 100, msg="For computation efficiency, convert control table to percentage values.")

  beta_table <- as.matrix(cbind(tumor_table, control_table))
  beta_table <- round(beta_table)
  storage.mode(beta_table) <- "integer"

  is_tumor <- c(rep(TRUE, ncol(tumor_table)), rep(FALSE, ncol(control_table)))

  # select rows by NAs
  message(sprintf("[%s] Selecting sites with fraction of NAs \u2264 %.2f...", Sys.time(), na_threshold))
  tumor_valid_sites <- rowSums(is.na(beta_table[,is_tumor]))/sum(is_tumor) <= na_threshold
  control_valid_sites <- rowSums(is.na(beta_table[,!is_tumor]))/sum(!is_tumor) <= na_threshold
  valid_sites_idx <- which(tumor_valid_sites & control_valid_sites)

  message(sprintf("[%s] Computing...", Sys.time()))
  auc <- setNames(rep(NA_real_, nrow(beta_table)), rownames(beta_table))
  cl <- parallel::makeCluster(ncores)
  auc[valid_sites_idx] <- parallel::parApply(cl, beta_table[valid_sites_idx,,drop=FALSE], 1, single_AUC, is_tumor = is_tumor)
  parallel::stopCluster(cl)

  message(sprintf("[%s] Done",  Sys.time()))
  return(auc)
}

#' Compute AUC a single vector
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
