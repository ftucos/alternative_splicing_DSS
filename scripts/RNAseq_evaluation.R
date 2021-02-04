setwd("/Volumes/TucosHDD/Bioinformatics/workspace/alternative_splicing_DSS")

library("tidyverse")
library("tidybulk")
library("tidyHeatmap")
options(scipen = 999)
# Import reads ------------------------------------
## list files to import
files <- list.files(path = "output/read_counts/", pattern = ".tab", full.names = TRUE)

## Define function to import reads
parse_counts <- function(x) {
  read_delim(x, "\t", skip = 4,
             col_names = c("transcript", "count", "count_fw_stranded", "count_rev_stranded")) %>%
    ## extract sample name
    mutate(sample = gsub(x, pattern = ".*/", replacement = "") %>% gsub(x, pattern = ".ReadsPerGene.out.tab", replacement = ""),
           ## remove version from ENS gene ID
           transcript = str_extract(transcript, "^[^.]+")) %>%
    select(sample, transcript, count)
}

## Import all reads
df <- map(files, parse_counts) %>%
  do.call(rbind, .)

## Import metadata
metadata <- read_csv("data/SraRunTable.txt")

selected_info <- metadata %>%
  select(Run, Age, agent, rna_isolation_date)

# TidyBulk ------------------------------------
df.1 <- df %>%
  left_join(selected_info, by=c("sample" = "Run")) %>%
  tidybulk(sample, transcript, count) %>%
# Normalize for read depths with "Trimmed Mean of M-values” method from edgeR
  scale_abundance(method = "TMM",)

df.2 <- df.1 %>%
  # Remove outlier sample
#  filter(sample != "SRR7529788") %>%
  keep_variable(top = 1000, log_transform = TRUE)

# automatically computes SD
# Based on complex_heatmap https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#heatmap-titles
df.2 %>%
  heatmap(transcript, sample, count_scaled, palette_value = c("red", "white", "blue"), show_row_dend = FALSE) %>%
  add_tile(rna_isolation_date) %>%
  add_tile(Age) %>%
  add_tile(agent) 
  
  