Biocmanager::install("microbiome")
install.packages("ggVennDiagram")

#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("ibd_final.RData")

#### "core" microbiome ####

# Convert to relative abundance
ibd_RA <- transform_sample_counts(ibd_final, fun=function(x) x/sum(x))

# Filter dataset by antibiotic use
IBD_fmt_short <- subset_samples(ibd_RA, `treatment_type`=="IBD_fmt_short")
IBD_fmt_long <- subset_samples(ibd_RA, `treatment_type`=="IBD_fmt_long")
IBD_no_fmt <- subset_samples(ibd_RA, `treatment_type`=="IBD_no_fmt")
IBD_surg_long <- subset_samples(ibd_RA, `treatment_type`=="IBD_surg_long")
IBD_surg_short <- subset_samples(ibd_RA, `treatment_type`=="IBD_surg_short")
IBD_no_surg <- subset_samples(ibd_RA, `treatment_type`=="IBD_no_surg")
Healthy_control <- subset_samples(ibd_RA, `treatment_type`=="Healthy_control")

# What ASVs are found in more than x% of samples in each treatment type category?
# trying changing the prevalence to see what happens
IBD_fmt_short_ASV <- core_members(IBD_fmt_short, detection=0.001, prevalence = 0.1)
IBD_fmt_long_ASV <- core_members(IBD_fmt_long , detection=0.001, prevalence = 0.1)
IBD_no_fmt_ASV <- core_members(IBD_no_fmt, detection=0.001, prevalence = 0.1)
IBD_surg_long_ASV <- core_members(IBD_surg_long, detection=0.001, prevalence = 0.1)
IBD_surg_short_ASV <- core_members(IBD_surg_short , detection=0.001, prevalence = 0.1)
IBD_no_surg_ASV <- core_members(IBD_no_surg , detection=0.001, prevalence = 0.1)
Healthy_control_ASV <- core_members(Healthy_control , detection=0.001, prevalence = 0.1)
# What are these ASVs? you can code it in two different ways to see the same things
prune_taxa(IBD_fmt_short_ASV,ibd_final) %>%
  tax_table()
prune_taxa(IBD_fmt_long_ASV,ibd_final) %>%
  tax_table()
prune_taxa(IBD_no_fmt_ASV,ibd_final) %>%
  tax_table()
prune_taxa(IBD_surg_long_ASV,ibd_final) %>%
  tax_table()
prune_taxa(IBD_surg_short_ASV,ibd_final) %>%
  tax_table()
prune_taxa(IBD_no_surg_ASV,ibd_final) %>%
  tax_table()
prune_taxa(Healthy_control_ASV,ibd_final) %>%
  tax_table()

# plot ASVs' relative abundance
IBD_fmt_short_ASV_plot <- prune_taxa(IBD_fmt_short_ASV,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`treatment_type`, scales ="free")
IBD_fmt_long_ASV_plot <- prune_taxa(IBD_fmt_long_ASV,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`treatment_type`, scales ="free")
IBD_no_fmt_ASV_plot <- prune_taxa(IBD_no_fmt_ASV,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`treatment_type`, scales ="free")
IBD_surg_long_ASV_plot <- prune_taxa(IBD_surg_long_ASV,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`treatment_type`, scales ="free")
IBD_no_surg_ASV_plot <- prune_taxa(IBD_no_surg_ASV,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`treatment_type`, scales ="free")
IBD_surg_short_ASV_plot <- prune_taxa(IBD_surg_short_ASV,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`treatment_type`, scales ="free")
Healthy_control_ASV_plot <- prune_taxa(Healthy_control_ASV,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`treatment_type`, scales ="free")
# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
IBD_fmt_short_list <- core_members(IBD_fmt_short, detection=0.001, prevalence = 0.10)
IBD_fmt_long_list <- core_members(IBD_fmt_long, detection=0.001, prevalence = 0.10)
IBD_no_fmt_list <- core_members(IBD_no_fmt, detection=0.001, prevalence = 0.10)
IBD_surg_short_list <- core_members(IBD_surg_short, detection=0.001, prevalence = 0.10)
IBD_surg_long_list <- core_members(IBD_surg_long, detection=0.001, prevalence = 0.10)
IBD_no_surg_list <- core_members(IBD_no_surg, detection=0.001, prevalence = 0.10)
Healthy_control_list <- core_members(Healthy_control, detection=0.001, prevalence = 0.10)



IBD_list_full <- list(IBD_fmt_short = IBD_fmt_short_list, IBD_fmt_long = IBD_fmt_long_list, IBD_no_fmt = IBD_no_fmt_list,IBD_surg_short = IBD_surg_short_list, IBD_surg_long = IBD_surg_long_list, IBD_no_surg = IBD_no_surg_list, Healthy_control = Healthy_control_list )
IBD_fmt_list_full <- list(IBD_fmt_short = IBD_fmt_short_list, IBD_fmt_long = IBD_fmt_long_list )
IBD_fmt_list_full_Healthy_contol <- list(IBD_fmt_short = IBD_fmt_short_list, IBD_fmt_long = IBD_fmt_long_list, Healthy_control = Healthy_control_list )
IBD_fmt_list_full_no_fmt <- list(IBD_fmt_short = IBD_fmt_short_list, IBD_fmt_long = IBD_fmt_long_list, IBD_no_fmt = IBD_no_fmt_list )
IBD_surg_list_full <- list( IBD_surg_short = IBD_surg_short_list, IBD_surg_long = IBD_surg_long_list)
IBD_surg_list_full_Healthy_contol <- list(IBD_surg_short = IBD_surg_short_list, IBD_surg_long = IBD_surg_long_list, Healthy_control = Healthy_control_list )
IBD_surg_FMT_list_full <- list(IBD_fmt_short = IBD_fmt_short_list, IBD_fmt_long = IBD_fmt_long_list,IBD_surg_short = IBD_surg_short_list, IBD_surg_long = IBD_surg_long_list)
IBD_surg_list_full_no_surg <- list(IBD_surg_short = IBD_surg_short_list, IBD_surg_long = IBD_surg_long_list, IBD_no_surg = IBD_no_surg_list )
IBD_no_surg_no_fmt_healthy_control <- list(IBD_no_surg = IBD_no_surg_list, IBD_no_fmt = IBD_no_fmt_list, Healthy_control = Healthy_control_list  )

#graphs used for manuscript and presentation 
IBD_fmt_list_full_no_fmt_venn <- ggVennDiagram(x = IBD_fmt_list_full_no_fmt)
ggsave("venn_IBD_fmt_list_full_no_fmt_1.png", IBD_fmt_list_full_no_fmt_venn)

IBD_surg_list_full_no_surg_venn <- ggVennDiagram(x = IBD_surg_list_full_no_surg)
ggsave("venn_IBD_surg_list_full_no_surg_1.png", IBD_surg_list_full_no_surg_venn)

#extra graphs 
IBD_list_full_venn <- ggVennDiagram(x = IBD_list_full)
ggsave("venn_IBD_list_full_1.png", IBD_list_full_venn)

IBD_fmt_list_full_venn <- ggVennDiagram(x = IBD_fmt_list_full)
ggsave("venn_IBD_fmt_list_full_1.png", IBD_fmt_list_full_venn)

IBD_fmt_list_full_healthy_control_venn <- ggVennDiagram(x = IBD_fmt_list_full_Healthy_contol)
ggsave("venn_IBD_fmt_list_full_healthy_control_1.png", IBD_fmt_list_full_healthy_control_venn)


IBD_surg_list_full_venn <- ggVennDiagram(x = IBD_surg_list_full)
ggsave("venn_IBD_surg_list_full_1.png", IBD_surg_list_full_venn)

IBD_surg_list_full_healthy_control_venn <- ggVennDiagram(x = IBD_surg_list_full_Healthy_contol)
ggsave("venn_IBD_surg_list_full_healthy_control_1.png", IBD_surg_list_full_healthy_control_venn)

IBD_surg_FMT_list_full_venn <- ggVennDiagram(x = IBD_surg_FMT_list_full)
ggsave("venn_IBD_surg_FMT_list_full_1.png", IBD_surg_FMT_list_full_venn)

IBD_no_surg_no_fmt_healthy_control_venn <- ggVennDiagram(x = IBD_no_surg_no_fmt_healthy_control)
ggsave("venn_IBD_no_surg_no_fmt_healthy_control_1.png", IBD_no_surg_no_fmt_healthy_control_venn)

