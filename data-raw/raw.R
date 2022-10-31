library(readr)
library(dplyr)
library(magrittr)
library(devtools)

illumina27k_hg19 <- read_tsv("illumina27k_hg19.txt.gz")
illumina27k_hg19 %<>% select(Chromosome, Genomic_Coordinate, Probe=`Composite Element REF`, Gene_Symbol)
illumina27k_hg19$Genomic_Coordinate[illumina27k_hg19$Genomic_Coordinate==0] <- NA
illumina27k_hg19 <- illumina27k_hg19[startsWith(illumina27k_hg19$Probe, "cg"),]

illumina450k_hg19 <- read_tsv("illumina450k_hg19.txt.gz")
illumina450k_hg19 %<>% select(Chromosome, Genomic_Coordinate, Probe=`Composite Element REF`, Gene_Symbol)
illumina450k_hg19$Genomic_Coordinate[illumina450k_hg19$Genomic_Coordinate==0] <- NA
illumina450k_hg19 <- illumina450k_hg19[startsWith(illumina450k_hg19$Probe, "cg"),]

illumina27k_hg38 <- read_tsv("illumina27k_hg38.txt.gz")
illumina27k_hg38 %<>% select(Chromosome, Genomic_Coordinate=Start, Probe=`Composite Element REF`, Gene_Symbol:Feature_Type)
illumina27k_hg38$Genomic_Coordinate[illumina27k_hg38$Genomic_Coordinate==0] <- NA
illumina27k_hg38 <- illumina27k_hg38[startsWith(illumina27k_hg38$Probe, "cg"),]

illumina450k_hg38 <- read_tsv("illumina450k_hg38.txt.gz")
illumina450k_hg38 %<>% select(Chromosome, Genomic_Coordinate=Start, Probe=`Composite Element REF`, Gene_Symbol:Feature_Type)
illumina450k_hg38$Genomic_Coordinate[illumina450k_hg38$Genomic_Coordinate==0] <- NA
illumina450k_hg38 <- illumina450k_hg38[startsWith(illumina450k_hg38$Probe, "cg"),]

use_data(illumina27k_hg19, overwrite=T)
use_data(illumina450k_hg19, overwrite=T)
use_data(illumina27k_hg38, overwrite=T)
use_data(illumina450k_hg38, overwrite=T)

load("../data/cpg_islands.rda")
N <- 1000 # first N islands
bs_seq_tumor_toy_matrix <- c()
bs_seq_control_toy_matrix <- c()
cpg_sites <- c()
nsamples <- 10 # columns
cpg_islands_names <- c()
i <- 1
for (i in seq_len(N)){
  island <- cpg_islands_df[i,]
  cpg_islands_names <- c(cpg_islands_names, paste(island, collapse = "-"))
  min_x <- island$start
  max_x <- island$end
  set.seed(i)
  nrows <- sample(3:6, 1) # number of rows within an island
  set.seed(i)
  pos <- sample(seq(min_x, max_x), nrows)
  set.seed(i)
  m <- sample(c(0.1, 0.2, 0.5, 0.8, 0.9), 1) # mean score
  set.seed(i)
  tumors <- matrix(abs(runif(nsamples*nrows, m-0.05, m+0.05)), byrow = TRUE, nrow = nrows)
  set.seed(i*-1)
  m <- sample(c(0.1, 0.2, 0.5, 0.8, 0.9), 1)
  set.seed(i*-1)
  controls <- matrix(abs(runif(nsamples*nrows, m-0.05, m+0.05)), nrow = nrows)
  bs_seq_tumor_toy_matrix <- rbind(bs_seq_tumor_toy_matrix, tumors)
  bs_seq_control_toy_matrix <- rbind(bs_seq_control_toy_matrix, controls)
  cpg_sites <- c(cpg_sites, pos)
}

bs_seq_toy_sites <- tibble(chr = "1", pos = cpg_sites)
dimnames(bs_seq_tumor_toy_matrix) <- list(paste0("chr1_", cpg_sites), paste0("tumor", seq_len(nsamples)))
dimnames(bs_seq_control_toy_matrix) <- list(paste0("chr1_", cpg_sites), paste0("control", seq_len(nsamples)))
bs_seq_toy_matrix <- round(cbind(bs_seq_tumor_toy_matrix, bs_seq_control_toy_matrix)*100)
use_data(bs_seq_toy_matrix, bs_seq_toy_sites, overwrite = TRUE)
