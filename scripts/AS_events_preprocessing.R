#=================================
# To do:
# - filter significant event
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


fill_in_missing_events <- function(x) {
  data.frame(event_name = rep(unique(x$event_name),
                              times = length(unique(x$sample))),
             sample = rep(unique(x$sample),
                          each = length(unique(x$event_name)))) %>%
    # Add event type field
    left_join(x %>%
                select(event_type, event_name) %>%
                distinct())
}


# Import files ----------------------------------
# Import DSS AS events
## List DSS analysis files
DSS_path <- "processed/miso/summary/"
DSS_files <- list.files(path=DSS_path, pattern = "*.miso_summary")
DSS_metadata <- read_csv("data/SraRunTable.txt")

DSS_AS <- map(DSS_files, ~import_AS(.x, DSS_path)) %>%
  do.call(rbind, .) %>%
  filter(event_type != "isoform") %>%
  mutate(miso_posterior_mean = as.numeric(miso_posterior_mean),
         ci_low = as.numeric(ci_low),
         ci_high = as.numeric(ci_high)) %>%
  # Add sample metadata
  left_join(select(DSS_metadata, Run, Treatment=agent), by=c("sample" = "Run")) %>%
  mutate(Run = sample,
         sample = paste0(sample, "-", Treatment)) %>%
  # Add 0-counts_events
  full_join(fill_in_missing_events(.)) 

# Import Rbm3 analysis files
Rbm3_path <- "/Volumes/TucosHDD/Bioinformatics/workspace/alternative_splicing_siRbm3/processed/miso/summary/"
Rbm3_files <- list.files(path=Rbm3_path, pattern = "*.miso_summary")
Rbm3_metadata <- read_csv("/Volumes/TucosHDD/Bioinformatics/workspace/alternative_splicing_siRbm3/data/SraRunTable.txt")

Rbm3_AS <- map(Rbm3_files, ~import_AS(.x, Rbm3_path)) %>%
  do.call(rbind, .) %>%
  filter(event_type != "isoform") %>%
  mutate(miso_posterior_mean = as.numeric(miso_posterior_mean),
         ci_low = as.numeric(ci_low),
         ci_high = as.numeric(ci_high)) %>%
  left_join(select(Rbm3_metadata, Run, Treatment), by=c("sample" = "Run")) %>%
  mutate(Run = sample,
         sample = paste0(Run, "-", str_extract(Treatment, "(mock|siRbm3)"))) %>%
  full_join(fill_in_missing_events(.))

# X + Y >= N, Y >= 1
# in at least one of the RNA-Seq samples processed by MISO, where X, Y are the number of reads in the (1,0) and (0,1) read classes,
# respectively, and N is an arbitrary but sizeable number (e.g. 10 or 20).
# This filter requires that the sum of inclusion and exclusion-supporting reads be greater than or equal to N,
# and that the read class supporting exclusion is non-zero in at least one of the samples.
# This guarantees that at least some of the reads are informative about the inclusive or exclusive isoform.
# Note that this filter does not guarantee that there will be junction evidence for the inclusion isoform,
# but does guarantee that there will be a junction evidence for skipping of the exon.

AS = rbind(DSS_AS %>% mutate(experiment = "DSS"),
           Rbm3_AS %>% mutate(experiment = "siRbm3"))

# Add low quality event flag
N = 20

AS.1 <- AS %>%
  mutate(X =  as.numeric(str_extract(counts, "(?<=\\(0,1\\):)[0-9]+")),
         X = replace_na(X, 0),
         Y =  as.numeric(str_extract(counts, "(?<=\\(1,0\\):)[0-9]+")),
         Y = replace_na(Y, 0),
         Status = ifelse(X + Y >= N, "Non-missing", "Missing")    
  )
  
AS.1 %>%
  ggplot(aes(x=sample, fill=Status)) +
  geom_bar() +
  theme_bw() +
  ylab("Events Number") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  facet_wrap(~experiment,scales = "free_x") +
  labs(title = "Alternative Splicing events",
       subtitle = "",
       caption = "DSS: read length = 50 bp, read depth ~ 32 M\nsiRbm3: read length = 75 bp, read depth ~ 40 M\nThreshold: noAS_reads + AS_reads >= 20")

AS.2 <- AS.1 %>%
  group_by(experiment, event_name) %>%
  # remove events missing in more than 1 sample
  filter(sum(Status == "Non-missing", na.rm = TRUE) >= length(unique(sample))-1)

AS.2 %>%
  ggplot(aes(x=sample, fill=Status)) +
  geom_bar() +
  ylab("Events Number") + 
  theme_bw() +
  ylab("Events Number") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  facet_wrap(~experiment, scales = "free_x") +
  labs(title = "Alternative Splicing events filtered",
       subtitle = "Filter: event missing in max one sample",
       caption = "DSS: read length = 50 bp, read depth ~ 32 M\nsiRbm3: read length = 75 bp, read depth ~ 40 M\nThreshold: noAS_reads + AS_reads >= 20")

                    
# Differential events analysis --------------------------------
DE_DSS <- AS.2 %>%
  filter(experiment == "DSS",
         Status == "Non-missing") %>%
  group_by(event_name, event_type, Treatment) %>%
  summarize(mean_PSI = mean(miso_posterior_mean),
            min_CI_low = min(ci_low),
            max_CI_high = max(ci_high)) %>%
  ungroup()
  
DE_DSS.1 <- full_join(
  DE_DSS %>%
    filter(Treatment == "Water") %>%
    select(-Treatment) %>%
    rename(mean_PSI.water = mean_PSI,
           min_CI_low.water = min_CI_low,
           max_CI_high.water = max_CI_high),
  DE_DSS %>%
    filter(Treatment == "DSS") %>%
    select(-Treatment) %>%
    rename(mean_PSI.dss = mean_PSI,
           min_CI_low.dss = min_CI_low,
           max_CI_high.dss = max_CI_high),
) %>%
  # Identify Differential Events
  mutate(PSI_diff = mean_PSI.dss - mean_PSI.water,
         direction = case_when(
           PSI_diff >= 0.1 & min_CI_low.dss > max_CI_high.water ~ "UP",
           PSI_diff <- 0.1  & max_CI_high.dss < min_CI_low.water ~ "DOWN",
           TRUE ~ "NON SIGNIFICANT")
         )


DE_Rbm3 <- AS.2 %>%
  filter(experiment == "siRbm3",
         Status == "Non-missing") %>%
  group_by(event_name, event_type, Treatment) %>%
  summarize(mean_PSI = mean(miso_posterior_mean),
            min_CI_low = min(ci_low),
            max_CI_high = max(ci_high)) %>%
  ungroup()

DE_Rbm3.1 <- full_join(
  DE_Rbm3 %>%
    filter(Treatment == "mock transfected at 37C") %>%
    select(-Treatment) %>%
    rename(mean_PSI.mock = mean_PSI,
           min_CI_low.mock = min_CI_low,
           max_CI_high.mock = max_CI_high),
  DE_Rbm3 %>%
    filter(Treatment == "siRbm3 transfected at 37C") %>%
    select(-Treatment) %>%
    rename(mean_PSI.si = mean_PSI,
           min_CI_low.si = min_CI_low,
           max_CI_high.si = max_CI_high),
) %>%
  # Identify Differential Events
  mutate(PSI_diff = mean_PSI.si - mean_PSI.mock,
         direction = case_when(
           PSI_diff >= 0.1 & min_CI_low.si > max_CI_high.mock ~ "UP",
           PSI_diff <= -0.1  & max_CI_high.si < min_CI_low.mock ~ "DOWN",
           TRUE ~ "NON SIGNIFICANT")
  )


# Plot
rbind(
  DE_DSS.1 %>% select(event_name, event_type, PSI_diff, direction) %>% mutate(experiment = "DSS"),
  DE_Rbm3.1 %>% select(event_name, event_type, PSI_diff, direction) %>% mutate(experiment = "siRbm3")) %>%
  group_by(experiment,  event_type, direction) %>%
  summarize(count = n()) %>%
  rbind(data.frame(
    experiment = c("DSS", "DSS"),
    event_type = c(NA, NA),
    direction = c("UP", "DOWN"),
    count = c(0, 0)
  )) %>%
  filter(direction != "NON SIGNIFICANT") %>%
  ggplot(aes(x=as.factor(experiment), y=count, fill = event_type)) +
  geom_col(position="stack") +
  theme_bw() +
  facet_wrap(~direction) +
  scale_fill_discrete(na.translate=FALSE, name = "AS event type") +
  ggtitle("Significantly altered AS events") +
  xlab("")

# Find overlapping events
a <- DE_Rbm3.1 %>%
  group_by(event_name) %>%
  summarize(count = n())
