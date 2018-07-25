#' Compute Area Under Curve for a matrix of samples
#'
#' This function computes the Area Under Curve that is then used to define the
#' segregation between tumor and control samples accordingly to their
#' methylation (beta) values.
#'
#' @param tumor_table A matrix of beta-values (percentage) from tumor samples.
#' @param control_table A matrix of beta-values (percentage) from
#' normal/control samples.
#' @param ncores Number of parallel processes to use for parallel computing
#' @param na_threshold Fraction of NAs (considered independently in tumor and
#' control samples) above which a site will not be selected (default=0).
#' @return A vector of AUC scores.
#' @export
#' @examples
#' auc_data <- compute_AUC(tumor_toy_data, control_toy_data)
#' auc_data_bs <- compute_AUC(bs_toy_matrix[,1:10], bs_toy_matrix[,11:20])
compute_AUC <- function(tumor_table, control_table, ncores=1, na_threshold=0,
  max_NAs_frac, tumor, control){
  # check parameters
  if (!missing(max_NAs_frac)){
    warning("'max_NAs_frac' is deprecated. Use 'na_threshold' instead.")
    na_threshold <- max_NAs_frac
  }
  if (!missing(tumor) || !missing(control)){
    warning(paste0("'tumor' and 'control' are deprecated.
      Use 'tumor_table' and 'control_table' instead."))
      tumor_table <- tumor
      control_table <- control
  }
  ncores <- as.integer(ncores)
  system_cores <- parallel::detectCores()
  assertthat::assert_that(ncores < system_cores)

  na_threshold <- as.numeric(na_threshold)
  assertthat::assert_that(na_threshold >= 0 || na_threshold < 1)

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
  tumor_NAs <- apply(beta_table[,sample_state], 1, function(x) {
    sum(is.na(x))
  })/sum(sample_state) <= na_threshold
  control_NAs <- apply(beta_table[,!sample_state], 1, function(x) {
    sum(is.na(x))
  })/sum(!sample_state) <= na_threshold

  NAs_idx <- which(tumor_NAs & control_NAs)

  auc <- rep(NA_real_, nrow(beta_table))
  message(sprintf("[%s] Compute AUC", Sys.time()))
  cl <- parallel::makeCluster(ncores)
  auc[NAs_idx] <- parallel::parApply(cl, beta_table[NAs_idx,], 1,
                                     single_AUC, states = sample_state)
  parallel::stopCluster(cl)
  message(sprintf("[%s] Done",  Sys.time()))
  return(auc)
}

#' Compute AUC a single vector
#'
#' Use Wilcoxon method to compute AUC.
#' Return NA if NA samples are more than threshold
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
