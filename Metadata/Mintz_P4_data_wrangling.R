library(tidyverse)

# Load data
metadata <- read_delim("data/fmt_p4_metadata.txt", delim="\t")

# Remove CDI 
metadata <- metadata %>%
  filter(condition != "CDI")

# Rename columns to match Halfvarson metadata
metadata <- metadata %>%
  rename(ibd_subtype = condition) %>%
  mutate(ibd_subtype = case_when(
    grepl("Donor$", ibd_subtype) ~ "Healthy_control",
    TRUE ~ ibd_subtype))

# Filter names to match Halfverson paper
metadata <- metadata %>%
  mutate(treatment_type = case_when(
    grepl("Donor$", Group) ~ "Healthy_control",
    grepl("3MoPST$", Group) ~ "IBD_fmt_long",
    grepl("PRE$", Group) ~ "IBD_no_fmt",
    grepl("1WkPST$", Group) ~ "IBD_fmt_short")) %>%
  select(-c("Group", "timepoint"))

# Adding diagnosis column
metadata_final <- metadata %>%
  mutate(diagnosis = case_when(ibd_subtype %in% c("CDI+UC", "UC") ~ "IBD",
                               ibd_subtype == "Healthy_control" ~ "Healthy_control"))

# Save metadata as a .tsv file
write.table(metadata_final, file="data/Mintz_P4_filtered.tsv", sep="\t", quote = T, row.names = F, col.names = T)
