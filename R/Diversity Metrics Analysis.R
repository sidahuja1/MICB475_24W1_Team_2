library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(paletteer)
library(ggsignif)
library(patchwork)
library(microViz)
library(vegan)
library(FSA)
library(ggthemes)

# Load in phyloseq data
load("ibd_rare.RData")
load("ibd_final.RData")

# Set Seed for reproducbility
set.seed(1)

### ALPHA DIVERSITY ###

# Generating alpha diversity metrics and merging with sample data
alphadiv <- estimate_richness(ibd_rare)
samp_dat <- sample_data(ibd_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)


## Observed 

# Running Kruskal-Wallis multiple comparison test to see if any group differs
# in observed diversity
kruskal.test(Observed ~ `treatment_type`, data = samp_dat_wdiv)

# Running Dunn's test to see for pairwise comparisons
dunnTest(Observed ~ treatment_type,data = samp_dat_wdiv, method = "bonferroni")

# Renaming groups for figure axis
samp_dat_wdiv <- samp_dat_wdiv %>%
  mutate(
    treatment_formatted = case_when(
      grepl("Healthy_control", treatment_type, ignore.case = TRUE) ~ "Healthy Control",
      grepl("IBD_no_fmt", treatment_type, ignore.case = TRUE) ~ "No FMT IBD",
      grepl("IBD_fmt_short", treatment_type, ignore.case = TRUE) ~ "FMT Short",
      grepl("IBD_fmt_long", treatment_type, ignore.case = TRUE) ~ "FMT Long",
      grepl("IBD_no_surg", treatment_type, ignore.case = TRUE) ~ "No Surgery IBD",
      grepl("IBD_surg_short", treatment_type, ignore.case = TRUE) ~ "Surgery Short",
      grepl("IBD_surg_long", treatment_type, ignore.case = TRUE) ~ "Surgery Long"))


# Setting factor order for ordering of box plot
samp_dat_wdiv$treatment_formatted <- factor(samp_dat_wdiv$treatment_formatted, 
                                       levels = c("Healthy Control", 
                                                  "No FMT IBD", "FMT Short", 
                                                  "FMT Long", "No Surgery IBD",
                                                  "Surgery Short", "Surgery Long"))

# Categorizing groups for coloring boxes in figure
samp_dat_wdiv <- samp_dat_wdiv %>%
  mutate(
    Group = case_when(
      grepl("Healthy", treatment_type, ignore.case = TRUE) ~ "Healthy",
      grepl("fmt", treatment_type, ignore.case = TRUE) ~ "FMT",
      grepl("surg", treatment_type, ignore.case = TRUE) ~ "Surgery"))

# Setting factor order for ordering of legend
samp_dat_wdiv <- samp_dat_wdiv %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Healthy", "FMT", "Surgery")))

# ALPHA DIVERSITY --> Observed plot
obs_plot <- ggplot(samp_dat_wdiv, aes(x=`treatment_formatted`, y=Observed)) +
  geom_boxplot(aes(fill = Group), alpha = 0.5) +
  xlab("Treatment") +
  geom_signif(comparisons = list(c("Healthy Control","No FMT IBD"), 
                                 c("Healthy Control", "Surgery Long"),
                                 c("Healthy Control", "Surgery Short")),
              y_position = c(480, 530, 505),
              tip_length = 0.01,
              vjust = 0.4,
              annotations = c("***","****","****")) +
  theme_classic()
obs_plot


## Shannon 

# Running Kruskal-Wallis multiple comparison test to see if any group differs
# in Shannon diversity
kruskal_obs <- kruskal.test(Shannon ~ `treatment_type`, data = samp_dat_wdiv)

# Running Dunn's test to see for pairwise comparisons
dunnTest(Shannon ~ treatment_type,data = samp_dat_wdiv, method = "bonferroni")

# ALPHA DIVERSITY --> Shannon plot
shan_plot <- ggplot(samp_dat_wdiv, aes(x=treatment_formatted, y=Shannon)) +
  geom_boxplot(aes(fill = Group), alpha = 0.5) +
  xlab("Treatment") +
  theme_classic() +
  geom_signif(comparisons = list(c("Healthy Control", "Surgery Long"),
                                 c("Healthy Control", "Surgery Short")),
              y_position = c(5.3, 5.0),
              tip_length = 0.01,
              vjust = 0.4,
              annotations = c("***","***"))
shan_plot


## Faith's phylogenetic diversity 

# Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(ibd_rare)), phy_tree(ibd_rare),
                 include.root=F) 

# Add PD to metadata table
sample_data(ibd_rare)$PD <- phylo_dist$PD

samp_dat <- sample_data(ibd_rare)
samp_dat <- data.frame(samp_dat)

# Running Kruskal-Wallis multiple comparison test to see if any group differs
# in Shannon diversity
kruskal_obs <- kruskal.test(PD ~ `treatment_type`, data = samp_dat)

# Running Dunn's test to see for pairwise comparisons
dunnTest(PD ~ treatment_type, data = samp_dat, method = "bonferroni")

# Renaming groups for figure axis
samp_dat <- samp_dat %>%
  mutate(
    treatment_formatted = case_when(
      grepl("Healthy_control", treatment_type, ignore.case = TRUE) ~ "Healthy Control",
      grepl("IBD_no_fmt", treatment_type, ignore.case = TRUE) ~ "No FMT IBD",
      grepl("IBD_fmt_short", treatment_type, ignore.case = TRUE) ~ "FMT Short",
      grepl("IBD_fmt_long", treatment_type, ignore.case = TRUE) ~ "FMT Long",
      grepl("IBD_no_surg", treatment_type, ignore.case = TRUE) ~ "No Surgery IBD",
      grepl("IBD_surg_short", treatment_type, ignore.case = TRUE) ~ "Surgery Short",
      grepl("IBD_surg_long", treatment_type, ignore.case = TRUE) ~ "Surgery Long"))


# Setting factor order for ordering of box plot
samp_dat$treatment_formatted <- factor(samp_dat$treatment_formatted, 
                                            levels = c("Healthy Control", 
                                                       "No FMT IBD", "FMT Short", 
                                                       "FMT Long", "No Surgery IBD",
                                                       "Surgery Short", "Surgery Long"))

# Categorizing groups for coloring boxes in figure
samp_dat <- samp_dat %>%
  mutate(
    Group = case_when(
      grepl("Healthy", treatment_type, ignore.case = TRUE) ~ "Healthy",
      grepl("fmt", treatment_type, ignore.case = TRUE) ~ "FMT",
      grepl("surg", treatment_type, ignore.case = TRUE) ~ "Surgery"))

# Setting factor order for ordering of legend
samp_dat <- samp_dat %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Healthy", "FMT", "Surgery")))

# ALPHA DIVERSITY --> PD plot
plot.pd <- ggplot(samp_dat, aes(treatment_formatted, PD)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5) +
  xlab("Treatment") +
  theme_classic() +
  ylab("Phylogenetic Diversity") +
  geom_signif(comparisons = list(c("Healthy Control", "No FMT IBD"),
                                 c("Healthy Control", "Surgery Long"),
                                 c("Healthy Control", "Surgery Short")),
              y_position = c(30, 33, 31.5),
              tip_length = 0.01,
              vjust = 0.4,
              annotations = c("**","***", "***"))
plot.pd


### BETA DIVERSITY ###


sample_data(ibd_rare) <- samp_dat

wuf_dm <- UniFrac(ibd_rare, weighted = T)
pcoa_wuf <- ordinate(ibd_rare, method="PCoA", distance='unifrac')

adonis2(wuf_dm ~ Group, data=samp_dat)

gg_pcoa <- plot_ordination(ibd_rare, pcoa_wuf, color = "treatment_formatted") +
  labs(col = "Treatment") +
  theme_few()+
  stat_ellipse(aes(color = Group), linetype = 2, type = 'norm')+
  geom_point(size = 3, alpha = 0.5) +
  scale_color_paletteer_d("ggthemes::Classic_10") +
  annotate("text", x=0.35, y=-0.5, label= "P < 0.01", size = 5)
gg_pcoa
ggsave('gg_pcoa.png', width = 6.5, height = 5, device='png', dpi=300)
#### Taxonomy bar plots #####

pseq <- ibd_rare %>%
  tax_fix() %>%
  phyloseq_validate()

taxa_names_df <- as.data.frame(tax_table(pseq))
taxa_names_df$Phylum <- gsub("^p__", "", taxa_names_df$Phylum)
tax_table(pseq) <- as.matrix(taxa_names_df)

taxa <- pseq %>%
  phyloseq::merge_samples(group = "treatment_formatted") %>%
  comp_barplot(tax_level = "Phylum", n_taxa = 12, bar_width = 0.8, 
               sample_order = rev(c("Healthy Control", 
                                "No FMT IBD", 
                                "FMT Short", 
                                "FMT Long", 
                                "No Surgery IBD",
                                "Surgery Short", 
                                "Surgery Long"))) +
  labs(x = 'Treatment', y = 'Average Abundance') +
  coord_flip()
taxa
ggsave('taxa.png', width = 8, height = 4, device='png', dpi=300)


# multi <- (obs_plot + shan_plot) /plot_spacer() + (plot.pd + taxa)+
#   plot_annotation(tag_levels = 'A') +
#   plot_layout(heights = c(1, 0.08, 1), widths = c(0.5,0.5)) &
#   theme(plot.tag = element_text(face = 'bold', size = 15)) +
#   theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
#         axis.text=element_text(size=9),
#         axis.title=element_text(size=11.5))
# multi

multi <- (obs_plot + plot_spacer() + shan_plot + plot_spacer() + plot.pd)+
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect', ncol = (5), widths = c(1,0.1, 1, 0.1, 1))&
  theme(plot.tag = element_text(face = 'bold', size = 15),
        axis.text.x = element_text(angle = 40, vjust =1, hjust=1),
        axis.text=element_text(size=9),
        axis.title=element_text(size=11.5))
multi

ggsave('multi.png', width = 10, height = 4.5, device='png', dpi=300)
