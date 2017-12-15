library(readr)
library(dplyr)
library(devtools)

illumina27k_hg19 <- read_tsv("data-raw/illumina27k_hg19.txt.gz") %>%
    select(-Beta_value) %>%
    mutate(Genomic_Coordinate=ifelse(Genomic_Coordinate == 0, NA, Genomic_Coordinate))
illumina27k_hg38 <- read_tsv("data-raw/illumina27k_hg38.txt.gz") %>%
    select(-Beta_value) %>%
    mutate(Genomic_Coordinate=ifelse(Genomic_Coordinate == 0, NA, Genomic_Coordinate))
illumina450k_hg19 <- read_tsv("data-raw/illumina450k_hg19.txt.gz") %>%
    select(-Beta_value) %>% filter(startsWith(`Composite Element REF`, "cg"))
illumina450k_hg38 <- read_tsv("data-raw/illumina450k_hg38.txt.gz") %>%
    select(-Beta_value) %>% filter(startsWith(`Composite Element REF`, "cg"))

use_data(illumina27k_hg19, overwrite=T)
use_data(illumina450k_hg19, overwrite=T)
use_data(illumina27k_hg38, overwrite=T)
use_data(illumina450k_hg38, overwrite=T)
