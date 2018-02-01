library(readr)
library(dplyr)
library(devtools)

# illumina27k_hg19 <- read_tsv("data-raw/illumina27k_hg19.txt.gz") %>%
#     select(-Beta_value) %>%
#     mutate(Genomic_Coordinate=ifelse(Genomic_Coordinate == 0, NA, Genomic_Coordinate))
# illumina27k_hg38 <- read_tsv("data-raw/illumina27k_hg38.txt.gz") %>%
#     select(-Beta_value)
# illumina450k_hg19 <- read_tsv("data-raw/illumina450k_hg19.txt.gz") %>%
#     select(-Beta_value) %>% filter(startsWith(`Composite Element REF`, "cg"))
# illumina450k_hg38 <- read_tsv("data-raw/illumina450k_hg38.txt.gz") %>%
#     select(-Beta_value) %>% filter(startsWith(`Composite Element REF`, "cg"))
#
# use_data(illumina27k_hg19, overwrite=T)
# use_data(illumina450k_hg19, overwrite=T)
# use_data(illumina27k_hg38, overwrite=T)
# use_data(illumina450k_hg38, overwrite=T)

N <- 1000 # firt N islands
bs_tumor_toy_data <- c()
bs_control_toy_data <- c()
cpg_sites <- c()
nsamples <- 10 # columns
cpg_islands_names <- c()
i <- 1
for (i in seq_len(N)){
  island = cpg_islands[i,]
  cpg_islands_names <- c(cpg_islands_names, paste(island, collapse = "-"))
  min_x = island$start
  max_x = island$end
  set.seed(i)
  n = sample(3:6, 1) # number of rows within an island
  set.seed(i)
  pos = sample(seq(min_x, max_x), n)
  set.seed(i)
  m = sample(c(0.1, 0.2, 0.5, 0.8, 0.9), 1) # mean score
  set.seed(i)
  tumors <- matrix(abs(rnorm(nsamples*n, m, 0.05)), byrow = T, nrow = n)
  set.seed(i*-1)
  m = sample(c(0.1, 0.2, 0.5, 0.8, 0.9), 1)
  set.seed(i*-1)
  controls <- matrix(abs(rnorm(nsamples*n, m, 0.05)), nrow = n)
  bs_tumor_toy_data <- rbind(bs_tumor_toy_data, tumors)
  bs_control_toy_data <- rbind(bs_control_toy_data, controls)
  cpg_sites <- c(cpg_sites, pos)
}
bs_toy_coordinates <- data_frame(chr = "chr1", pos = cpg_sites)
dimnames(bs_tumor_toy_data) <- list(paste0("chr1_", cpg_sites), paste0("tumor", seq_len(nsamples)))
dimnames(bs_control_toy_data) <- list(paste0("chr1_", cpg_sites), paste0("control", seq_len(nsamples)))
bs_toy_matrix <- cbind(bs_tumor_toy_data, bs_control_toy_data)
bs_toy_indexes <- compute_islands_indexes(bs_toy_coordinates, head(cpg_islands, N))
names(bs_toy_indexes) <- cpg_islands_names
use_data(bs_toy_matrix, bs_toy_coordinates, bs_toy_indexes, overwrite = T)
