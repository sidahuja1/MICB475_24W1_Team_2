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
