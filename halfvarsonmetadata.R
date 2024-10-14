library(tidyverse)
library(dplyr)

# Load data
metadata <- read_delim("halfvarson_metadata.tsv", delim="\t")

# Select columns of interest
metadata_select <- select(metadata, "sample-id", "cd_resection", "collection_timestamp", "diagnosis_full", 
                          "ibd_subtype", "patientnumber", "perianal_disease", "timepoint")


# Filter out LC and CC ibd subtypes
metadata_filter <- metadata_select %>%
  filter(ibd_subtype != "LC") %>%
  filter(ibd_subtype != "CC")

# Obtain timepoints of interest (1 and 2)
cd_timepoint_1_2 <- metadata_filter %>%
  filter(diagnosis_full == "CD" & timepoint %in% c(1, 2)) %>%
  group_by(timepoint) %>%
  summarize(count = n())

patients_with_timepoint_1_2 <- metadata_filter %>%
  filter(timepoint %in% c(1, 2)) %>%
  group_by(patientnumber) %>%
  filter(all(c(1, 2) %in% timepoint)) %>%
  ungroup()

# Check counts of patients
count_HC <- patients_with_timepoint_1_2 %>%
  filter(diagnosis_full == "HC") %>%
  nrow()/2

count_UC <- patients_with_timepoint_1_2 %>%
  filter(diagnosis_full == "UC") %>%
  nrow()/2

# Combine IBD subtype column into one column
metadata_final <- patients_with_timepoint_1_2 %>%
  mutate(diagnosis = case_when(
    diagnosis_full %in% c("UC", "CD") ~ "IBD",
    diagnosis_full == "HC" ~ "HC",
    TRUE ~ diagnosis_full
  ))

# Select necessary columns
metadata_final <- metadata_final %>%
  select("sample-id", "cd_resection", "ibd_subtype", "patientnumber", "timepoint", "diagnosis")

# Save metadata as a .tsv file
write.table(metadata_final, file="halfvarson_metadata_filtered.tsv", sep="\t", quote = T, row.names = F, col.names = T)
