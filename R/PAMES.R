#' PAMES: Purity Assessment from clonal MEthylation Sites
#'
#' The PAMES package provides a set of functions to estimate
#' the level of purity or the tumor content of tumor samples.
#'
#' The basic workflow of PAMES requires to \code{\link{get_AUC}} (to evaluate tumor-control methylation differences),
#' \code{\link{find_informative_sites}} (to retrieve sites of interest),
#' and \code{\link{get_purity}} of tumor samples.
#' When working with methylation data obtained with other technologies (such as Bisulphite Sequencing),
#' users should must map their set of CpG sites to differentially methylated regions
#' (such as CpG islands) using data \code{\link{reduce_to_regions}}, then
#' \code{\link{get_AUC}}, \code{\link{find_informative_regions}} and finally
#' \code{\link{get_purity}}.
#
#' @docType package
#' @name PAMES
NULL
