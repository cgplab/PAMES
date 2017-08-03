#' Convert beta values from Bisulphite Sequencing to single CpG island beta values
#'
#' Remap beta values from different CpG sites to single beta values
#' associated to one CpG island where CpG sites are located.
#'
#' Bisulphite Sequencing is used to retrieve level of DNA methylation at
#' base resolution. Different experiments though may retrieve
#' different CpG sites. This function makes a direct comparison possible.
#'
#' @param bs_data a matrix of beta values from Bisulphite Sequencing data.
#' @param sites_df a data.frame reporting the location of CpG sites from Bisulphite
#' Sequencing data. It requires two columns, (1) chromosome (format: `chrN`)
#' and (2) 1-based coordinate of CpG site.
#' @param islands_df a data.frame reporting the location of CpG islands. By default it
#' loads the `cpg_islands` data.frame that comes with the package (based on data
#' downloaded from [UCSC Genome Browser](https://genome.ucsc.edu/) in 2016-03-01.
#' @param min_CpGs an integer (default to 3). Minimum number of CpG sites within
#' a single CpG island
#' required to compute the mapped beta value (return NA if below that number).
#' @return a matrix of beta values (one row per CpG island).
#' @examples
#' \donttest{
#' mapped <- map_to_islands(bs_control_toy_data, bs_toy_sites)
#' }
#' \dontrun{
#' mapped <- map_to_islands(bs_tumor_toy_data, bs_toy_data_sites,
#'                islands_df = cpg_islands_new, min_CpGs=5)
#' }
#' @export
map_to_islands <- function(bs_data, sites_df, islands_df=cpg_islands, min_CpGs=3){
    bs_data <- as.matrix(bs_data)
    if (!is.data.frame(sites_df) | !is.data.frame(islands_df))
        stop("'sites_df' and 'islands_df' must be data.frames")
    if (!is.character(sites_df[[1]]) | !is.numeric(sites_df[[2]]))
        stop("First to colums of 'sites_df' must be 'character' and 'numeric', respectively")
    if (!is.character(islands_df[[1]]) | !is.numeric(islands_df[[2]])| !is.numeric(islands_df[[3]]))
        stop("First column of 'islands_df' must be 'character', second and third 'numeric'")

    message(sprintf("[%s] Mapping CpG sites to CpG islands...",  Sys.time()))
    cpg_indexes <- vector("list", nrow(islands_df))
    #### insert time bar
    pb <- txtProgressBar(min=0, max=length(cpg_indexes), style=3)
    for (i in seq_len(nrow(islands_df))){
        cpg_indexes[[i]] <- compute_indexes(sites_df, islands_df[i, ])
        setTxtProgressBar(pb, i)
    }
    close(pb)

    # running time magnitude of following step should be in terms of minutes
    message(sprintf("[%s] Reducing beta values...",  Sys.time()))
    mapped <- do.call("rbind", lapply(cpg_indexes, function(idxs){
        if (length(idxs) < min_CpGs){
            beta_values <- rep(NA, ncol(bs_data))
        } else {
            beta_values <- apply(bs_data[idxs, , drop=F], 2, median, na.rm=T)
        }
        return(beta_values)
    }))
    message(sprintf("[%s] Done",  Sys.time()))
    return(mapped)
}

compute_indexes <- function(sites_df, one_island){
    same_chromosome  <- sites_df[[1]] == one_island[[1]]
    after_start      <- sites_df[[2]] >= one_island[[2]]
    before_end       <- sites_df[[2]] <= one_island[[3]]
    return(which(same_chromosome & after_start & before_end))
}

# load("/projects/data/methdata/lin2013_beltran2016/Beta_table_Lin2013_Beltran2015_BeltranRapidAutopsies2016.RData")
# time estimates:
# 1e4 bisulphite sequencing; full cpg_islands; 95 seconds elapsed
# 1e5 bisulphite sequencing; full cpg_islands; 151 seconds elapsed
# 6e6 bisulphite sequencing; full cpg_islands; 1h10m seconds elapsed
