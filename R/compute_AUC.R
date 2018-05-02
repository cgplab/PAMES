#' Compute Area Under Curve for a matrix of samples
#'
#' This function computes the Area Under Curve that is then used to define the
#' segregation between tumor and control samples accordingly to their
#' methylation (beta) values.
#'
#' @param tumor_table A matrix of beta-values (percentage) from tumor samples.
#' @param control_table A matrix of beta-values (percentage) from
#' normal/control samples.
#' @param na_threshold Fraction of NAs (considered independently in tumor and
#' control samples) above which a site will not be selected (default=0).
#' @param ncores Number of parallel processes to use for parallel computing
#' @return A vector of AUC scores.
#' @export
#' @examples
#' auc_data <- compute_AUC(tumor_toy_data, control_toy_data)
#' auc_data_bs <- compute_AUC(bs_toy_matrix[,1:10], bs_toy_matrix[,11:20])
compute_AUC <- function(tumor_table, control_table, na_threshold=1, ncores=1,
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
  assertthat::assert_that(na_threshold >= 0 || na_threshold <= 1)

  diff_range_t <- diff(range(tumor_table, na.rm = TRUE))
  diff_range_c <- diff(range(control_table, na.rm = TRUE))
  assertthat::assert_that(diff_range_t > 1 && diff_range_t <= 100 &&
    diff_range_c > 1 && diff_range_c <= 100,
    msg=paste("For computation efficiency, convert tumor and control tables",
        "to percentage values."))
  beta_table <- round(as.matrix(cbind(tumor_table, control_table)))
  storage.mode(beta_table) <- "integer"
  sample_state <- c(rep(TRUE, ncol(tumor_table)), rep(FALSE, ncol(control_table)))

  message(sprintf("[%s] Computing AUC ", Sys.time()))
  cl <- parallel::makeCluster(ncores)
  auc <- parallel::parApply(cl, beta_table, 1, single_AUC,
    state = sample_state, na_threshold = na_threshold)
  parallel::stopCluster(cl)
  message(sprintf("[%s] Done",  Sys.time()))
  return(auc)
}

#' Compute AUC a single vector
#'
#' Return NA if NA samples are more than threshold
#'
#' @param x integer vector (range 0-1)
#' @param state logical vector
#' @param na_threshold numeric
#' @keywords internal
single_AUC <- function(x, state, na_threshold) {
  stopifnot(is.numeric(x))
  all_NA <- all(is.na(x))
  tumor_NA <- sum(is.na(x[state])) / length(x[state]) > na_threshold
  control_NA <- sum(is.na(x[!state])) / length(x[!state]) > na_threshold
  if (all_NA || tumor_NA || control_NA ) {
    return(NA)
  } else {
    roc <- ROC::rocdemo.sca(state[!is.na(x)], x[!is.na(x)],
      cutpts = seq(0, 100, 1))
    return(ROC::AUC(roc))
  }
}
