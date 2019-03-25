#' Compute Area Under Curve for a matrix of samples
#'
#' This function computes the Area Under Curve that is then used to define the
#' segregation between tumor and control samples accordingly to their
#' methylation (beta) values.
#'
#' @param tumor_table A matrix of beta-values (percentage) from tumor samples.
#' @param control_table A matrix of beta-values (percentage) from
#' normal/control samples.
#' @param ncores Number of parallel processes to use for parallel computing.
#' @param na_threshold Fraction of NAs (considered independently in tumor and
#' control samples) above which a site will not be selected (default=0).
#' @return A vector of AUC scores.
#' @export
compute_AUC <- function(tumor_table, control_table, ncores = 1, na_threshold = 0) {
  # check parameters
  ncores <- as.integer(ncores)
  system_cores <- parallel::detectCores()
  assertthat::assert_that(ncores < system_cores)

  na_threshold <- as.numeric(na_threshold)
  assertthat::assert_that(na_threshold >= 0, na_threshold < 1)

  beta_table <- as.matrix(cbind(tumor_table, control_table))
  diff_range <- diff(range(beta_table, na.rm = TRUE))
  if (diff_range <= 1 || diff_range > 100) {
    stop(paste("For computation efficiency convert tumor and control",
               "tables to percentage values."))
  } else {
    beta_table <- round(beta_table)
    storage.mode(beta_table) <- "integer"
  }
  sample_state <- c(rep(TRUE, ncol(tumor_table)), rep(FALSE, ncol(control_table)))

  # select rows by NAs
  message(sprintf("[%s] Filter NA rows", Sys.time()))
  tumor_valid_sites <- rowSums(is.na(beta_table[,sample_state]))/sum(sample_state) <= na_threshold
  control_valid_sites <- rowSums(is.na(beta_table[,!sample_state]))/sum(!sample_state) <= na_threshold

  below_threshold_idx <- which(tumor_valid_sites & control_valid_sites)

  auc <- rep(NA_real_, nrow(beta_table))
  message(sprintf("[%s] Compute AUC", Sys.time()))
  cl <- parallel::makeCluster(ncores)
  auc[below_threshold_idx] <- parallel::parApply(cl, beta_table[below_threshold_idx,], 1, single_AUC, states = sample_state)
  parallel::stopCluster(cl)
  message(sprintf("[%s] Done",  Sys.time()))
  return(auc)
}

#' Compute AUC a single vector
#'
#' Use Wilcoxon method to compute AUC.
#'
#' @param scores integer vector (range 1-100)
#' @param states logical vector (class labels)
#' @keywords internal
#' \href{http://blog.revolutionanalytics.com/2017/03/auc-meets-u-stat.html}{http://blog.revolutionanalytics.com}
single_AUC <- function(scores, states) {
  assertthat::assert_that(is.integer(scores))
  assertthat::assert_that(is.logical(states))
  na_idx <- is.na(scores)
  scores <- scores[!na_idx]
  states <- states[!na_idx]
  n1 <- sum(states)
  n2 <- sum(!states)
  R1 <- sum(rank(scores)[states])
  U1 <- R1 - n1*(n1+1)/2
  return(U1/(n1*n2))
}
