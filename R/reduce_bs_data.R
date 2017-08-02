#' Reduce beta values from Bisulphite Sequencing to single CpG island scores
#'
#' Substitute beta values from a variable and unknown a priori number of CpG sites
#' to a fixed number of beta values equal to the number of curated CpG islands
#' (retireved from UCSC browser)
#'
#' Data from different BS experiment not necessarly retrieve DNA methylation
#' levels of the same CpG sites making impossible a direct comparison between
#' different experiments. To convert beta values from
#' different CpG sites located within a single CpG island to one single beta value
#' makes the comparison "doable".
#'
#' @param bs_data A matrix of beta values from Bisulphite Sequencing data
#' @param site_df A data.frame with the location of CpG sites from Bisulphite
#' Sequencing data. Requires two columns, `chr` and `pos` (respectively,
#' chromosome - in the form `chrN` and 1-based coordinate of CpG site)
#' @param island_df A data.frame with location of CpG islands. By default it
#' loads the `cpg_islands` data.frame that comes with the package, which can be
#' used as reference to create a customized data.frame
#' @param min_CpGs An integer (default to 3). Minimum number of CpGs within a single CpG island
#' required to compute the reduced beta value (return )
#' @return A matrix of beta values (one row per CpG island provided by
#' `island_df`)
#' \donottest{
#' @examples reduce_bs_data(bs_tumor_toy_data, bs_toy_data_sites)
#' }
#' \donotrun{
#' @examples reduce_bs_data(brca_bs_matrix, sites_df, cpg_islands_v3)
#' }
#' @export
reduce_bs_data <- function(bs_data, sites_df, islands_df=cpg_islands, min_CpGs=3){
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
    reduced <- do.call("rbind", lapply(cpg_indexes, function(idxs){
        if (length(idxs) < min_CpGs){
            beta_values <- rep(NA, ncol(bs_data))
        } else {
            beta_values <- apply(bs_data[idxs, , drop=F], 2, median, na.rm=T)
        }
        return(beta_values)
    }))
    message(sprintf("[%s] Done",  Sys.time()))
    return(reduced)
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
