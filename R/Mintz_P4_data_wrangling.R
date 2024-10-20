library(tidyverse)
library(dplyr)

metadata <- read_delim("fmt_p4_metadata.txt", delim="\t")

#remove CDI 
samples_to_remove <- c("CDI-PRE", "CDI-3MoPST", "CDI-1WkPST")
filtered_data <- metadata[!metadata$Group %in% samples_to_remove, ]
print(filtered_data)

metadata_p4_final_2 <- filtered_data %>%
  rename(treatment_type = Group)
print(metadata_p4_final_2)

#Filter names to match Halfverson paper##
metadata_p4_final_3<- metadata_p4_final_2 %>%
  mutate(treatment_type = ifelse(treatment_type == "Donor", "Healthy_control", treatment_type))
print(metadata_p4_final_3)
metadata_p4_final_4<- metadata_p4_final_3 %>%
  mutate(treatment_type = ifelse(treatment_type == "CDI+UC-3MoPST", "IBD_surg_long", treatment_type))
metadata_p4_final_4<- metadata_p4_final_4 %>%
mutate(treatment_type = ifelse(treatment_type == "UC-3MoPST", "IBD_surg_long", treatment_type))
metadata_p4_final_4<- metadata_p4_final_4 %>%
  mutate(treatment_type = ifelse(treatment_type == "CDI+UC-PRE", "IBD_no_surg", treatment_type))
metadata_p4_final_4<- metadata_p4_final_4 %>%
  mutate(treatment_type = ifelse(treatment_type == "UC-PRE", "IBD_no_surg", treatment_type))
metadata_p4_final_4<- metadata_p4_final_4 %>%
  mutate(treatment_type = ifelse(treatment_type == "UC-1WkPST", "IBD_surg_short", treatment_type))
metadata_p4_final_4<- metadata_p4_final_4 %>%
  mutate(treatment_type = ifelse(treatment_type == "CDI+UC-1WkPST", "IBD_surg_short", treatment_type))
print(metadata_p4_final_4)

# Save metadata as a .tsv file
write.table(metadata_p4_final_4, file="Mintz_P4_filtered.tsv", sep="\t", quote = T, row.names = F, col.names = T)
