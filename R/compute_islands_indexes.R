#' Associate CpG islands to corresponding CpG sites
#'
#' @param cpg_sites a data.frame reporting the location of CpG sites from Bisulphite
#' Sequencing data and requiring two columns: chromosome and 1-based coordinate of CpG sites.
#' @param cpg_islands_df a data.frame reporting the location of the CpG islands. By default it
#' loads the `cpg_islands` data.frame that comes with the package (based on data
#' downloaded from [UCSC Genome Browser](https://genome.ucsc.edu/) on date 2016-03-01.
#' @export
#' @examples
#' cpg_indexes_short <- compute_islands_indexes(bs_toy_sites, head(cpg_islands,10))
#' \dontrun{
#' cpg_indexes <- compute_islands_indexes(bs_toy_sites)
#' }
compute_islands_indexes <- function(cpg_sites, cpg_islands_df = cpg_islands){
  if (!is.data.frame(cpg_sites) | !is.data.frame(cpg_islands_df))
    stop("'cpg_sites' and 'cpg_islands_df' must be data.frames")
  if (!is.character(cpg_sites[[1]]) | !is.numeric(cpg_sites[[2]]))
    stop("First two columns of 'cpg_sites' must be 'character' and 'numeric', respectively")
  if (any(!startsWith(cpg_sites[[1]], "chr")))
    stop("Chromosome symbols in 'cpg_sites' must have the format 'chrN'")
  if (!is.character(cpg_islands_df[[1]]) | !is.numeric(cpg_islands_df[[2]])| !is.numeric(cpg_islands_df[[3]]))
    stop("First column of 'cpg_islands_df' must be 'character', second and third 'numeric'")

  message(sprintf("[%s] Computing indexes for 'CpG sites to CpG islands' mapping...",  Sys.time()))
  cpg_indexes <- vector("list", nrow(cpg_islands_df))
  for (i in seq_along(cpg_indexes)){
    cpg_island <- cpg_islands_df[i,]
    same_chromosome  <- cpg_sites[[1]] == cpg_island[["chr"]]
    is_between <- cpg_sites[[2]] >= cpg_island[["start"]] & cpg_sites[[2]] <= cpg_island[["end"]]
    cpg_indexes[[i]] <- which(same_chromosome & is_between)
  }
  message(sprintf("[%s] Done",  Sys.time()))
  return(cpg_indexes)
}
