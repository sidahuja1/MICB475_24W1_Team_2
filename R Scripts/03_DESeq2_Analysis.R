#Install the DESeq2 package
BiocManager::install("DESeq2")

#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)


#### Load data ####
load("~/Desktop/MICB475_Project/MICB_475_Final_Project_DESeq2/ibd_final.RData")


## For Reproducibility
set.seed(1)


## Adding pseudo-count of 1 to reads, converting data to a DESeq2 object & performing differential abundance analysis based on treatment type ## 
ibd_plus1 <- transform_sample_counts(ibd_final, function(x) x+1)
ibd_deseq <- phyloseq_to_deseq2(ibd_plus1, ~`treatment_type`)
DESEQ_ibd <- DESeq(ibd_deseq)



##### IBD_surg_long vs IBD_no_surg #####
ibd_res_surg_long_vs_no_surg <- results(DESEQ_ibd, tidy=TRUE, 
                                        contrast = c("treatment_type", "IBD_surg_long", "IBD_no_surg"))

surglong_nosurg_sigASVs <- ibd_res_surg_long_vs_no_surg %>% 
  filter(padj < 0.01 & abs(log2FoldChange) > 2.5) %>%
  dplyr::rename(ASV = row)
View(surglong_nosurg_sigASVs)
surglong_nosurg_sigASVs_vec <- surglong_nosurg_sigASVs %>% pull(ASV)

surglong_nosurg_DESeq <- prune_taxa(surglong_nosurg_sigASVs_vec, ibd_final)
surglong_nosurg_sigASVs <- tax_table(surglong_nosurg_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(surglong_nosurg_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(surglong_nosurg_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Count the number of increased, decreased, and total features
feature_counts <- surglong_nosurg_sigASVs %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Increased", "Decreased")) %>%
  group_by(Direction) %>%
  summarise(Count = n()) %>%
  bind_rows(data.frame(Direction = "Total", Count = nrow(surglong_nosurg_sigASVs)))

print(feature_counts)

# Save the plot
ggsave("surglong_nosurg_barplot.png", width = 16, height = 6, dpi = 300)


#### IBD_surg_short vs IBD_no_surg ####
ibd_res_surg_short_vs_no_surg <- results(DESEQ_ibd, tidy=TRUE, 
                                         contrast = c("treatment_type", "IBD_surg_short", "IBD_no_surg"))

surgshort_nosurg_sigASVs <- ibd_res_surg_short_vs_no_surg %>% 
  filter(padj < 0.01 & abs(log2FoldChange) > 2.5) %>%
  dplyr::rename(ASV = row)
View(surgshort_nosurg_sigASVs)
surgshort_nosurg_sigASVs_vec <- surgshort_nosurg_sigASVs %>% pull(ASV)

surgshort_nosurg_DESeq <- prune_taxa(surgshort_nosurg_sigASVs_vec, ibd_final)
surgshort_nosurg_sigASVs <- tax_table(surgshort_nosurg_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(surgshort_nosurg_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(surgshort_nosurg_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Count the number of increased, decreased, and total features
feature_counts <- surgshort_nosurg_sigASVs %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Increased", "Decreased")) %>%
  group_by(Direction) %>%
  summarise(Count = n()) %>%
  bind_rows(data.frame(Direction = "Total", Count = nrow(surgshort_nosurg_sigASVs)))

print(feature_counts)


# Save the plot
ggsave("surgshort_nosurg_barplot.png", width = 16, height = 6, dpi = 300)


#### IBD_surg_long vs IBD_surg_short ####
ibd_res_surg_long_vs_surg_short <- results(DESEQ_ibd, tidy=TRUE, 
                                           contrast = c("treatment_type", "IBD_surg_long", "IBD_surg_short"))

surglong_surgshort_sigASVs <- ibd_res_surg_long_vs_surg_short %>% 
  filter(padj < 0.01 & abs(log2FoldChange) > 2.5) %>%
  dplyr::rename(ASV = row)
View(surglong_surgshort_sigASVs)
surglong_surgshort_sigASVs_vec <- surglong_surgshort_sigASVs %>% pull(ASV)

surglong_surgshort_DESeq <- prune_taxa(surglong_surgshort_sigASVs_vec, ibd_final)
surglong_surgshort_sigASVs <- tax_table(surglong_surgshort_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(surglong_surgshort_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(surglong_surgshort_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Count the number of increased, decreased, and total features
feature_counts <- surglong_surgshort_sigASVs %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Increased", "Decreased")) %>%
  group_by(Direction) %>%
  summarise(Count = n()) %>%
  bind_rows(data.frame(Direction = "Total", Count = nrow(surglong_surgshort_sigASVs)))

print(feature_counts)

# Save the plot
ggsave("surglong_surgshort_barplot.png", width = 16, height = 6, dpi = 300)



#### IBD_fmt_long vs IBD_no_fmt ####
ibd_res_fmt_long_vs_no_fmt <- results(DESEQ_ibd, tidy=TRUE, 
                                      contrast = c("treatment_type", "IBD_fmt_long", "IBD_no_fmt"))

fmtlong_nofmt_sigASVs <- ibd_res_fmt_long_vs_no_fmt %>% 
  filter(padj < 0.01 & abs(log2FoldChange) > 2.5) %>%
  dplyr::rename(ASV = row)
View(fmtlong_nofmt_sigASVs)
fmtlong_nofmt_sigASVs_vec <- fmtlong_nofmt_sigASVs %>% pull(ASV)

fmtlong_nofmt_DESeq <- prune_taxa(fmtlong_nofmt_sigASVs_vec, ibd_final)
fmtlong_nofmt_sigASVs <- tax_table(fmtlong_nofmt_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(fmtlong_nofmt_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(fmtlong_nofmt_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Count the number of increased, decreased, and total features
feature_counts <- fmtlong_nofmt_sigASVs %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Increased", "Decreased")) %>%
  group_by(Direction) %>%
  summarise(Count = n()) %>%
  bind_rows(data.frame(Direction = "Total", Count = nrow(fmtlong_nofmt_sigASVs)))

print(feature_counts)

# Save the plot
ggsave("fmtlong_nofmt_barplot.png", width = 16, height = 6, dpi = 300)



#### IBD_fmt_short vs IBD_no_fmt ####
ibd_res_fmt_short_vs_no_fmt <- results(DESEQ_ibd, tidy=TRUE, 
                                       contrast = c("treatment_type", "IBD_fmt_short", "IBD_no_fmt"))

fmtshort_nofmt_sigASVs <- ibd_res_fmt_short_vs_no_fmt %>% 
  filter(padj < 0.01 & abs(log2FoldChange) > 2.5) %>%
  dplyr::rename(ASV = row)
View(fmtshort_nofmt_sigASVs)
fmtshort_nofmt_sigASVs_vec <- fmtshort_nofmt_sigASVs %>% pull(ASV)

fmtshort_nofmt_DESeq <- prune_taxa(fmtshort_nofmt_sigASVs_vec, ibd_final)
fmtshort_nofmt_sigASVs <- tax_table(fmtshort_nofmt_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(fmtshort_nofmt_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(fmtshort_nofmt_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Count the number of increased, decreased, and total features
feature_counts <- fmtshort_nofmt_sigASVs %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Increased", "Decreased")) %>%
  group_by(Direction) %>%
  summarise(Count = n()) %>%
  bind_rows(data.frame(Direction = "Total", Count = nrow(fmtshort_nofmt_sigASVs)))

print(feature_counts)

# Save the plot
ggsave("fmtshort_nofmt_barplot.png", width = 16, height = 6, dpi = 300)


#### IBD_fmt_long vs IBD_fmt_short ####
ibd_res_fmt_short_vs_fmt_long <- results(DESEQ_ibd, tidy=TRUE, 
                                         contrast = c("treatment_type", "IBD_fmt_long", "IBD_fmt_short"))

fmtshort_fmtlong_sigASVs <- ibd_res_fmt_short_vs_fmt_long %>% 
  filter(padj < 0.01 & abs(log2FoldChange) > 2.5) %>%
  dplyr::rename(ASV = row)
View(fmtshort_fmtlong_sigASVs)
fmtshort_fmtlong_sigASVs_vec <- fmtshort_fmtlong_sigASVs %>% pull(ASV)

fmtshort_fmtlong_DESeq <- prune_taxa(fmtshort_fmtlong_sigASVs_vec, ibd_final)
fmtshort_fmtlong_sigASVs <- tax_table(fmtshort_fmtlong_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(fmtshort_fmtlong_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(fmtshort_fmtlong_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Count the number of increased, decreased, and total features
feature_counts <- fmtshort_fmtlong_sigASVs %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Increased", "Decreased")) %>%
  group_by(Direction) %>%
  summarise(Count = n()) %>%
  bind_rows(data.frame(Direction = "Total", Count = nrow(fmtshort_fmtlong_sigASVs)))

print(feature_counts)

# Save the plot
ggsave("fmtshort_fmtlong_barplot.png", width = 16, height = 6, dpi = 300)
