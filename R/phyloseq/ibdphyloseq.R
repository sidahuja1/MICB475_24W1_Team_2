library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

meta <- read_delim("merged-metadata.tsv", delim="\t", skip=0)
meta <- meta[-1, ]

otu <- read_delim(file = "merge-export/table_export/feature-table.txt", delim="\t", skip=1)

tax <- read_delim(file = "merge-export/taxonomy_export/taxonomy.tsv", delim="\t")

phylotree <- read.tree("merge-export/tree_export/tree.nwk")

# Format OTU table
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

# Format sample metadata
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)
class(SAMP)

# Format taxonomy
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

# Create phyloseq object
ibd <- phyloseq(OTU, SAMP, TAX, phylotree)

otu_table(ibd)
sample_data(ibd)
tax_table(ibd)
phy_tree(ibd)

# Remove non-bacterial sequences
ibd_filt <- subset_taxa(ibd,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Remove ASVs that have less than 5 counts total
ibd_filt_nolow <- filter_taxa(ibd_filt, function(x) sum(x)>5, prune = TRUE)

# Remove samples with less than 100 reads
ibd_final <- prune_samples(sample_sums(ibd_filt_nolow)>100, ibd_filt_nolow)

# Rarefy samples
rarecurve(t(as.data.frame(otu_table(ibd_final))), cex=0.1)
ibd_rare <- rarefy_even_depth(ibd_final, rngseed = 1, sample.size = 44459)

# Saving
save(ibd_final, file="ibd_final.RData")
save(ibd_rare, file="ibd_rare.RData")
