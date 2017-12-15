library(PAMES)
library(dplyr)
load("/projects/data/methdata/lin2013_beltran2016/Beta_table_Lin2013_Beltran2015_BeltranRapidAutopsies2016.RData")
compute_indexes <- function(bs_sites, cpg_islands_df = cpg_islands){
  if (!is.data.frame(bs_sites) | !is.data.frame(cpg_islands_df))
    stop("'bs_sites' and 'cpg_islands_df' must be data.frames")
  if (!is.character(bs_sites[[1]]) | !is.numeric(bs_sites[[2]]))
    stop("First to colums of 'bs_sites' must be 'character' and 'numeric', respectively")
  if (any(startsWith(bs_sites[[1]], "chr")))
    stop("Please be sure that chromosome symbols in 'bs_sites' don't contain 'chr' string")
  if (!is.character(cpg_islands_df[[1]]) | !is.numeric(cpg_islands_df[[2]])| !is.numeric(cpg_islands_df[[3]]))
    stop("First column of 'cpg_islands_df' must be 'character', second and third 'numeric'")

  message(sprintf("[%s] Computing indexes for 'CpG sites to CpG islands' mapping...",  Sys.time()))
  cpg_indexes <- apply(cpg_islands_df, 1, function(x) {
    same_chromosome  <- bs_sites[[1]] == x[["chr"]]
    after_start      <- bs_sites[[2]] >= x[["start"]]
    before_end       <- bs_sites[[2]] <= x[["end"]]
    return(which(same_chromosome & after_start & before_end))
  })
}
x = rownames(beta_table)[1:10] %>% 
  strsplit("\\.") %>% 
  do.call(what= "rbind") %>% 
  as.data.frame() %>% 
  mutate(V1 = sub("^chr", "", V1))

indexes = compute_indexes(x[1:1e3,])

# load("/projects/data/methdata/lin2013_beltran2016/Beta_table_Lin2013_Beltran2015_BeltranRapidAutopsies2016.RData")
# time estimates:
# 1e4 bisulphite sequencing; full cpg_islands; 95 seconds elapsed
# 1e5 bisulphite sequencing; full cpg_islands; 151 seconds elapsed
# 6e6 bisulphite sequencing; full cpg_islands; 1h10m seconds elapsed
