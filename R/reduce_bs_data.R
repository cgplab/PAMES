#' Reduce beta values from Bisulphite Sequencing to single CpG island scores
#'
#' @param bs_data a matrix of beta values from Bisulphite Sequencing data
#' @param cpg_site_df a data.frame with the location of CpG sites from Bisulphite
#' Sequencing data. Requires two columns, "chr" and "pos" (respectively,
#' chromosome - in the form `chr1` and 1-based coordinate of CpG site).
#' @param cpg_island_df a data.frame with location of CpG islands. Set to "default"
#' to use default one or provide a different one. Must contain "chr", "start",
#' "end" (1-based coordinates).
#' @return a reduced matrix of beta values
#' @examples reduce_bs_data(brca_bs_matrix, sites_df)
#' @examples reduce_bs_data(brca_bs_matrix, sites_df, islands_df=cpg_islands_v3)
reduce_bs_data <- function(bs_data, sites_df, islands_df=cpg_islands){
    bs_data <- as.matrix(bs_data)
    if (!is.data.frame(sites_df) | !is.data.frame(islands_df))
        stop("'sites_df' and 'islands_df' must be data.frames")
    if (!is.character(sites_df[[1]]) | !is.numeric(sites_df[[2]]))
        stop("First to colums of 'sites_df' must be 'character' and 'numeric', respectively")
    if (!is.character(islands_df[[1]]) | !is.numeric(islands_df[[2]])| !is.numeric(islands_df[[3]]))
        stop("First column of 'islands_df' must be 'character', second and third 'numeric'")

    cat(sprintf("[%s] Mapping CpG sites to CpG islands...\n",  Sys.time()))
    cpg_indexes <- vector("list", nrow(islands_df))
    #### insert time bar
    pb <- txtProgressBar(min=0, max=length(cpg_indexes), style=3)
    for (i in seq_along(cpg_indexes)){
        same_chromosome  <- sites_df[[1]] == unlist(islands_df[i, 1])  # unlist if using a tbl_df
        after_start      <- sites_df[[2]] >= unlist(islands_df[i, 2])  # unlist if using a tbl_df
        before_end       <- sites_df[[2]] <= unlist(islands_df[i, 3])  # unlist if using a tbl_df
        cpg_indexes[[i]] <- which(same_chromosome & after_start & before_end)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    # running time magnitude of following step should be in terms of minutes
    cat(sprintf("[%s] Reducing beta values...\n",  Sys.time()))
    reduced <- do.call("rbind", lapply(cpg_indexes, function(idxs){
        if (length(idxs) == 0){
            beta_values <- rep(NA, ncol(bs_data))
        } else {
            beta_values <- apply(bs_data[idxs, , drop=F], 2, median, na.rm=T)
        }
        return(beta_values)
    }))
    cat(sprintf("[%s] Done\n",  Sys.time()))
    return(reduced)
}
# load("/projects/data/methdata/lin2013_beltran2016/Beta_table_Lin2013_Beltran2015_BeltranRapidAutopsies2016.RData")
# time estimates:
# 1e4 bisulphite sequencing; full cpg_islands; 95 seconds elapsed
# 1e5 bisulphite sequencing; full cpg_islands; 151 seconds elapsed
# 6e6 bisulphite sequencing; full cpg_islands; 1h10m seconds elapsed
