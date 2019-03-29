#' PAMES: Purity Assessment from clonal MEthylation Sites
#'
#' The PAMES package provides a set of functions to estimate
#' the level of purity of tumor samples.
#'
#' The basic workflow of PAMES requires to \code{\link{compute_AUC}} (to evaluate tumor-control methylation differences),
#' \code{\link{select_informative_sites}} (to retreive sites of interest),
#' and \code{\link{compute_purity}} of tumor samples.
#' When working with methylation data obtained with other technologies (such as Bisulphite Sequencing),
#' users should must map their set of CpG sites to differentially methylated regions
#' (such as CpG islands) using data \code{\link{reduce_to_regions}}, then
#' \code{\link{compute_AUC}}, \code{\link{select_informative_regions}} and finally
#' \code{\link{compute_purity}}.
#
#' @docType package
#' @name PAMES
NULL
