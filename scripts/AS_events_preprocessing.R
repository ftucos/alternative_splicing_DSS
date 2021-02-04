#=================================
# To do:
# - Import metadata
# - filter significant event
# - Plot gene expression to group samples and exclude the faulty one
#================================

library(tidyverse)
setwd("/Volumes/TucosHDD/Bioinformatics/workspace/alternative_splicing_DSS")

# Functions ------------------------------------------------------------
# Function required for parsing .miso_summary files
import_AS <- function(x, path) {
  read_delim(paste0(path, x), delim ="\t") %>%
    mutate(file_name = gsub(x, pattern=".miso_summary", replacement = ""),
           # Extract the SRA accession ID from the first part of the file name
           sample = str_extract(file_name, "^[^_]+"),
           event_type = str_extract(file_name, "(?<=_).*")) %>%
    select(sample, event_type, -file_name, everything())
}

# Import files ----------------------------------
# Import DSS AS events
## List DSS analysis files
DSS_path <- "processed/miso/summary/"
DSS_files <- list.files(path=DSS_path, pattern = "*.miso_summary")
DSS_metadata <- read_csv("data/SraRunTable.txt")

DSS_AS <- map(DSS_files, ~import_AS(.x, DSS_path)) %>%
  do.call(rbind, .)

# Import Rbm3 analysis files
Rbm3_path <- "/Volumes/TucosHDD/Bioinformatics/workspace/alternative_splicing_siRbm3/processed/miso/summary/"
Rbm3_files <- list.files(path=Rbm3_path, pattern = "*.miso_summary")
Rbm3_metadata <- read_csv("/Volumes/TucosHDD/Bioinformatics/workspace/alternative_splicing_siRbm3/data/SraRunTable.txt")

Rbm3_AS <- map(Rbm3_files, ~import_AS(.x, Rbm3_path)) %>%
  do.call(rbind, .)


# Sandbox ================================================
# Evaluate total number of possible events

events.files <- list.files(path = "data/miso_indices", pattern = "mm10_exon_index_*/*.gff")
  
 





DSS_AS %>%
  ggplot(aes(x=sample, fill=event_type)) +
    geom_bar()

DSS_AS %>% group_by(event_name) %>% summarize(count = n()) %>% pull(count) %>% hist()

a<-DSS_AS %>% group_by(event_name) %>% summarize(count = n()) 

# List siRbm3 


Rbm3_AS %>%
  filter(event_type != "isoform") %>%
  ggplot(aes(x=sample, fill=event_type)) +
  geom_bar()
