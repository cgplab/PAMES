library(readr)
library(dplyr)

illumina27k_hg19 <- read_tsv("data-raw/illumina27k_hg19.txt.gz") %>%
    select(-Beta_value)
illumina27k_hg38 <- read_tsv("data-raw/illumina27k_hg38.txt.gz") %>%
    select(-Beta_value)
illumina450k_hg19 <- read_tsv("data-raw/illumina450k_hg19.txt.gz") %>%
    select(-Beta_value) %>% filter(startsWith(`Composite Element REF`, "cg"))
illumina450k_hg38 <- read_tsv("data-raw/illumina450k_hg38.txt.gz") %>%
    select(-Beta_value) %>% filter(startsWith(`Composite Element REF`, "cg"))
load("data-raw/cpg_islands.rda")

devtools::use_data(illumina27k_hg19,  illumina27k_hg38,  illumina450k_hg19,
                   illumina450k_hg38, cpg_islands, internal=T, overwrite=T)

