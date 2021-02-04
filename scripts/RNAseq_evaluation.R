setwd("/Volumes/TucosHDD/Bioinformatics/workspace/alternative_splicing_DSS")

library("tidyverse")
library("tidybulk")
library("tidyHeatmap")
library("ggrepel")
set.seed(42)
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
  filter(lowly_abundant == FALSE) %>%
  # sort by variance of log transformed values
  keep_variable(top = 1000, log_transform = TRUE)

# automatically computes SD
# Based on complex_heatmap https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#heatmap-titles
df.2 %>%
  heatmap(transcript, sample, count_scaled, palette_value = c("red", "white", "blue"), show_row_dend = FALSE) %>%
  add_tile(rna_isolation_date) %>%
  add_tile(Age) %>%
  add_tile(agent) 


# Cluster samples with the list of DEG from Mathijs list   -------------
gene_list <- read_csv("data/genes_DSS_absLFC_0.5.txt", col_names = "symbol")

# Remap ENSG to gene symbol
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ensembl2symbol <- getBM(values = df.1$transcript,
                       filters = "ensembl_gene_id",
                       mart = mart,
                       attributes = c("external_gene_name", "ensembl_gene_id"))

df.symbol <- df.1 %>%
  left_join(ensembl2symbol, by = c("transcript" = "ensembl_gene_id"))

df.symbol %>%
  filter(external_gene_name %in% gene_list$symbol) %>%
  heatmap(external_gene_name, sample, count_scaled, palette_value = c("red", "white", "blue"), show_row_dend = FALSE) %>%
  add_tile(rna_isolation_date) %>%
  add_tile(Age) %>%
  add_tile(agent) %>%
  save_pdf("output/plots/heatmap_genes_of_interest.pdf", width = 20, height = 20, units = "cm")


# PCA plot ---------------------------
df.PCA <- df.1 %>%
  filter(lowly_abundant == FALSE) %>%
  reduce_dimensions(method="PCA", .dims = 2, scale = TRUE, log_transform = TRUE) %>%
  group_by(sample, agent) %>%
  # actually you could just pick the value from the first row,
  # it's the same in all the rows from the same sample
  summarize(PC1 = median(PC1), PC2 = median(PC2)) %>%
  ungroup()

ggplot(df.PCA, aes(x=PC1, y=PC2, group=sample, color = agent)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(aes(label = sample), size = 3, point.padding = unit(5, "cm"), show.legend = F) +
  ggtitle("PCA: no lowly expressed genes")+
  theme_bw()
