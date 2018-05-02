#' Associate CpG islands to corresponding CpG sites
#'
#' @param cpg_sites A data.frame reporting the location of CpG sites from Bisulphite
#' Sequencing data. It requires two columns: chromosome and 1-based coordinate of CpG sites.
#' @param cpg_islands A data.frame reporting the location of the CpG islands. By default it
#' loads the `cpg_islands_df` data.frame that comes with the package (based on data
#' downloaded from [UCSC Genome Browser](https://genome.ucsc.edu/) on date 2016-03-01.
#' @export
#' @examples
#' cpg_indexes_short <- compute_islands_indexes(bs_toy_sites, head(cpg_islands_df, 10))
#' \dontrun{
#' cpg_indexes <- compute_islands_indexes(bs_toy_sites)
#' }
compute_islands_indexes <- function(cpg_sites, cpg_islands = cpg_islands_df){
  assertthat::assert_that(is.data.frame(cpg_sites))
  assertthat::assert_that(is.data.frame(cpg_islands))
  assertthat::assert_that(is.character(cpg_sites[[1]]))
  assertthat::assert_that(is.numeric(cpg_sites[[2]]))
  assertthat::assert_that(all(startsWith(cpg_sites[[1]], "chr")))

  assertthat::assert_that(is.character(cpg_islands[[1]]))
  assertthat::assert_that(is.numeric(cpg_islands[[2]]))
  assertthat::assert_that(is.numeric(cpg_islands[[3]]))

  message(sprintf("[%s] Computing indexes for CpG island mapping...",
      Sys.time()))
  cpg_indexes <- vector("list", nrow(cpg_islands))
  for (i in seq_along(cpg_islands)){
    cpg_island <- cpg_islands[i,]
    same_chromosome  <- cpg_sites[[1]] == cpg_island[["chr"]]
    is_between <- cpg_sites[[2]] >= cpg_island[["start"]] & cpg_sites[[2]] <= cpg_island[["end"]]
    cpg_indexes[[i]] <- which(same_chromosome & is_between)
  }
  message(sprintf("[%s] Done",  Sys.time()))
  return(cpg_indexes)
}
