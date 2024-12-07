library(tidyverse)
library(phyloseq)
library(indicspecies)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

set.seed(1)

#### Load unrarefied data ####
load("ibd_final.RData")

# glom to Species
ibd_species <- tax_glom(ibd_final, "Species", NArm = FALSE)
ibd_species_RA <- transform_sample_counts(ibd_species, fun=function(x) x/sum(x))

#ISA for Species
isa_ibd_species <- multipatt(t(otu_table(ibd_species_RA)), cluster = sample_data(ibd_species_RA)$`treatment_type`)
summary(isa_ibd_species)

taxtable <- tax_table(ibd_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Create indicator species table for species
isa_species_results <- isa_ibd_species$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable, by = "ASV") %>%
  filter(p.value < 0.05) %>%
  filter(!is.na(Species)) %>%
  filter(!str_detect(Species, "uncultured|human|unidentified")) %>%
  filter(rowSums(across(everything(), ~ . == 1)) <= 1)

View(isa_species_results)

#### Selecting for ASVs corresponding to single treatment groups ####

treatment_columns <- c("s.Healthy_control", "s.IBD_fmt_short",
                           "s.IBD_fmt_long", 
                           "s.IBD_surg_long", "s.IBD_surg_short", 
                           "s.IBD_no_fmt", "s.IBD_no_surg")

#Adding a treatment type column
isa_species_results_table <- isa_species_results %>%
  mutate(`Treatment Type` = apply(
    select(., all_of(treatment_columns)), 1, function(row) {
      treatments <- names(row)[which(row == 1)]
      if (length(treatments) > 0) {
        paste(treatments, collapse = ", ")
      } else {
        NA_character_
      }
    }
  ))
View(isa_species_results_table)

#### Table Formatting ####
#Formatted species table
formatted_species_table <- isa_species_results_table %>%
  select(
    `Treatment Type`,
    Phylum, 
    Genus, 
    Species, 
    `Observed Indicator Value` = stat,  
    `P-value` = p.value
  ) %>%
  filter(!is.na(`Treatment Type`)) %>%  # Remove rows with NA in Treatment Type
  mutate(
    Phylum = str_remove(Phylum, "^p__"),    # Remove 'p__' prefix from Phylum
    Genus = str_remove(Genus, "^g__"),      # Remove 'g__' prefix from Genus
    Species = str_remove(Species, "^s__"), # Remove 's__' prefix from Species
    `Treatment Type` = str_remove(`Treatment Type`, "^s\\.IBD_"),
    `Treatment Type` = str_replace(`Treatment Type`, "no_fmt", "No FMT"),
    `Treatment Type` = str_replace(`Treatment Type`, "no_surg", "No Surgery"),
    `Treatment Type` = str_replace(`Treatment Type`, "fmt_long", "FMT Long"),
    `Treatment Type` = str_replace(`Treatment Type`, "fmt_short", "FMT Short"),
    `Treatment Type` = str_replace(`Treatment Type`, "surg_short", "Surgery Short"),
    `Treatment Type` = str_replace(`Treatment Type`, "surg_long", "Surgery Long"),
    `P-value` = round(`P-value`, 3),        # Round p.value to three decimal points
    `Observed Indicator Value` = round(`Observed Indicator Value`, 3)
  ) %>%
  arrange(`Treatment Type`)  # Sort by Treatment Type in ascending order
view(formatted_species_table)

write.csv(formatted_species_table, "indicator_species_results.csv", row.names = FALSE)


