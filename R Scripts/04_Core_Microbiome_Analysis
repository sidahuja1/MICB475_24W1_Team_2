Biocmanager::install("microbiome")
install.packages("ggVennDiagram")
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("ibd_final.RData")

#### "core" microbiome ####
#### Convert to relative abundance ####
ibd_RA <- transform_sample_counts(ibd_final, fun=function(x) x/sum(x))
ibd_RA

# Filter dataset by treatment_type
IBD_fmt_short <- subset_samples(ibd_RA, `treatment_type`=="IBD_fmt_short") 
IBD_fmt_long <- subset_samples(ibd_RA, `treatment_type`=="IBD_fmt_long")
IBD_no_fmt <- subset_samples(ibd_RA, `treatment_type`=="IBD_no_fmt")
IBD_surg_long <- subset_samples(ibd_RA, `treatment_type`=="IBD_surg_long")
IBD_surg_short <- subset_samples(ibd_RA, `treatment_type`=="IBD_surg_short")
IBD_no_surg <- subset_samples(ibd_RA, `treatment_type`=="IBD_no_surg")
Healthy_control <- subset_samples(ibd_RA, `treatment_type`=="Healthy_control")


#### Create ASV for each treatment type ####
#### Create ASV and list for each treatment type ####
####detection = 0.001 and prevalence = 0.1####
IBD_fmt_short_list <- core_members(IBD_fmt_short, detection=0.001, prevalence = 0.10)
IBD_fmt_long_list <- core_members(IBD_fmt_long, detection=0.001, prevalence = 0.10)
IBD_no_fmt_list <- core_members(IBD_no_fmt, detection=0.001, prevalence = 0.10)
IBD_surg_short_list <- core_members(IBD_surg_short, detection=0.001, prevalence = 0.10)
IBD_surg_long_list <- core_members(IBD_surg_long, detection=0.001, prevalence = 0.10)
IBD_no_surg_list <- core_members(IBD_no_surg, detection=0.001, prevalence = 0.10)
Healthy_control_list <- core_members(Healthy_control, detection=0.001, prevalence = 0.10)

####Create lists ####
IBD_fmt_list_full_Healthy_contol <- list(IBD_fmt_short = IBD_fmt_short_list, IBD_fmt_long = IBD_fmt_long_list, Healthy_control = Healthy_control_list )
IBD_surg_list_full_Healthy_contol <- list(IBD_surg_short = IBD_surg_short_list, IBD_surg_long = IBD_surg_long_list, Healthy_control = Healthy_control_list )

#### Create graphs used for manuscript and presentation ####
IBD_fmt_list_full_healthy_control_venn <- ggVennDiagram(x = IBD_fmt_list_full_Healthy_contol)
ggsave("venn_IBD_fmt_list_full_healthy_control_1.png", IBD_fmt_list_full_healthy_control_venn)

IBD_surg_list_full_healthy_control_venn <- ggVennDiagram(x = IBD_surg_list_full_Healthy_contol) 
ggsave("venn_IBD_surg_list_full_healthy_control_1.png", IBD_surg_list_full_healthy_control_venn)



