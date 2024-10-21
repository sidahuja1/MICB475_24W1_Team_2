library(tidyverse)

# Load data
metadata <- read_delim("data/IBD_halfvarson_metadata.tsv", delim="\t")

# Select columns of interest
metadata <- select(metadata, "sample-id", "cd_resection", "diagnosis_full", 
                   "ibd_subtype", "patientnumber", "timepoint")


# Filter out LC and CC ibd subtypes
metadata <- metadata %>%
  filter(ibd_subtype != "LC" & ibd_subtype != "CC") 

# Obtain timepoints of interest (1 and 2) and only keep patient ids with 
# both timepoint 1 and 2
metadata <- metadata %>%
  filter(timepoint %in% c(1, 2)) %>%
  group_by(patientnumber) %>%
  filter(n() == 2) %>%
  ungroup()

# Check counts of patients
count_HC <- metadata %>%
  filter(diagnosis_full == "HC") %>%
  nrow()/2 # 9

count_UC <- metadata %>%
  filter(diagnosis_full == "UC") %>%
  nrow()/2 # 48

count_CD <- metadata %>%
  filter(diagnosis_full == "CD") %>%
  nrow()/2 # 38

count_surgery <- metadata %>%
  filter(cd_resection == "yes") %>%
  nrow()/2 # 17

# Filtering for resection only and adding treatment type column
metadata_surg <- metadata %>%
  filter(cd_resection == "yes") %>%
  mutate(treatment_type = case_when(timepoint == 1 ~ "IBD_surg_short",
                               timepoint == 2 ~ "IBD_surg_long"))

# Filtering for healthy controls and only keeping timepoint 1 and adding 
# treatment type column
metadata_healthy <- metadata %>%
  filter(diagnosis_full == "HC" & timepoint == 1)  %>%
  mutate(treatment_type = "HC")

# Filtering for no resection IBD patients
metadata_no_surg <- metadata %>%
  filter(cd_resection != "yes" & diagnosis_full != "HC")

# Randomly selecting 6 samples of the 3 IBD subtypes with no surgery from
# timepoint 1. Adding treatment type column

set.seed(123) # Setting seed for reproducibility

metadata_no_surg <- metadata_no_surg %>%
  filter(timepoint == 1) %>%
  group_by(ibd_subtype) %>%
  slice_sample(n = 6) %>%
  mutate(treatment_type = "IBD_no_surg")
  
# Combining dataframes
metadata_final <- bind_rows(metadata_surg, metadata_no_surg, metadata_healthy)

# Adding short diagnosis description and removing resection classification
# from IBD subtype
metadata_final <- metadata_final %>%
  mutate(diagnosis = case_when(diagnosis_full %in% c("UC", "CD") ~ "IBD",
                               diagnosis_full == "HC" ~ "HC"))

metadata_final$ibd_subtype <- gsub("_(nr)$", "", metadata_final$ibd_subtype)
  
# Select necessary columns
metadata_final <- metadata_final %>%
  select("sample-id", "ibd_subtype", "diagnosis", "treatment_type") %>%
  mutate(across(everything(), ~ gsub("HC", "Healthy_control", .)))

# Save metadata as a .tsv file
write.table(metadata_final, file="data/IBD_halfvarson_metadata_filtered.tsv", sep="\t", quote = T, row.names = F, col.names = T)
